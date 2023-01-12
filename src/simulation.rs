use std::fs;
use std::io::Write;
use std::ops::Deref;
use std::sync::mpsc;
use std::sync::{Arc, Mutex};
use std::thread::{self, JoinHandle};

use rand::seq::IteratorRandom;
use rand::{Rng, thread_rng};

use crate::config::{self, SimulationConfig, GridConfig, PTConfig};
use crate::grid::{Grid,History, Sweep};


pub struct PTEnviroment{
    pub config: PTConfig,
    pub pt_ids: Vec<equalTGridIds>,
}

impl PTEnviroment {
    pub fn new(config:PTConfig)-> Self{
        return PTEnviroment {
            config,
            pt_ids: Vec::new()}
        }
    pub fn delete_id(&mut self, id:u32){
        for element in self.pt_ids.iter_mut(){
            element.delete_id(id)
        }
        self.pt_ids.retain(|elem|elem.ids.len() > 0);
    }

}

pub struct equalTGridIds{
    pub T:f64,
    pub was_switched: Vec<bool>,
    pub current_pt_acceptance_prob: f64,
    pub ids: Vec<u32>
}     

impl equalTGridIds{
    pub fn delete_id(&mut self, id:u32){
        self.ids.retain(|&existing_id|existing_id != id);
    }
    pub fn update_current_acceptance_prob(&mut self){
        let sum = self.was_switched.iter().filter(|&&elem|elem==true).count();
        self.current_pt_acceptance_prob = sum as f64 / self.was_switched.len() as f64;
    }
}


pub struct Simulation{
    pub config: SimulationConfig,
    pub default_grid_config: GridConfig,
    grids: Vec<Grid>,
    grids_to_delete: Vec<u32>,
    pub current_grid_ids: Vec<u32>,
    pub running: bool,
    pub pt_enviroment: PTEnviroment,
    pub overlapp_histo: Vec<u32>,
    pub linked_overlapp_histo: Vec<u32>,
}

impl Simulation{
    pub fn new()->Self{
        let global_config = config::load();
        let config = global_config.simulation_config.clone();
        let default_grid_config = global_config.grid_config.clone();
        let mut grids = Vec::new();
        let mut current_grid_ids : Vec<u32> = Vec::new();
        for grid_id in 0..config.num_grids{
            grids.push(Grid::new(default_grid_config.clone(),grid_id).unwrap());
            current_grid_ids.push(grid_id);
        }
        let mut grids_to_delete: Vec<u32> = Vec::new();
        let running = false;
        let mut pt_enviroment = PTEnviroment::new(global_config.pt_config.clone());

        let num_bins = (2.0/config.histo_width+1.0).floor() as usize;
        let overlapp_histo = vec![0;num_bins];
        let linked_overlapp_histo = vec![0;num_bins];
        
        let mut sim =  Simulation{
            config,
            default_grid_config,
            grids,
            current_grid_ids,
            grids_to_delete,
            running,
            pt_enviroment,
            overlapp_histo,
            linked_overlapp_histo
        };
        sim.init_dir();
        return sim
    }

    #[cfg(not(target_arch = "wasm32"))]
    fn init_dir(&self){
        fs::remove_dir_all("./results/");
        fs::create_dir("./results/");
    }

    #[cfg(target_arch = "wasm32")]
    fn init_dir(&self){
        return;
    }    

    pub fn simulation_step(&mut self){
        self.thread_step();
        self.pt_exchange();
        self.update_overlapp_histogramm(0.5);
    }

    pub fn custom_run(&mut self){
        let grid_id1 = self.new_grid().unwrap();
        self.pt_init(grid_id1).unwrap();
        self.running = true;
    }

    pub fn update_overlapp_histogramm(&mut self, target_T: f64){
        if self.running == false {return}
        match self.pt_enviroment.pt_ids.iter().find(|&elem|elem.T == target_T){
            None => return,
            Some(someGrids) =>{
                let ids = someGrids.ids.clone();
                if someGrids.ids.len()<2{return}
                let first_grid_id = someGrids.ids.first().unwrap();
                let mut other_grid_ids = someGrids.ids.iter().skip(1);
                while let Some(other_grid_id) = other_grid_ids.next(){
                    let overlapp = self.grid(first_grid_id.clone()).unwrap().calc_overlap(
                        self.grid(other_grid_id.clone()).unwrap()).unwrap();
                    let capacity = self.grid(first_grid_id.clone()).unwrap().capacity as f64;
                    let normed_orverlapp = overlapp/capacity;
                    let bin = ((normed_orverlapp+1.0)/self.config.histo_width).floor() as usize;
                    self.overlapp_histo[bin]+=1;

                    let linked_overlapp = self.grid(first_grid_id.clone()).unwrap().calc_linked_overlap(
                        self.grid(other_grid_id.clone()).unwrap()).unwrap();
                    
                    let num_links = 2.0*self.grid(first_grid_id.clone()).unwrap().dimensions().len() as f64;
                    let normed_linked_orverlapp = linked_overlapp/(num_links*capacity);
                    let bin = ((normed_linked_orverlapp+1.0)/self.config.histo_width).floor() as usize;
                    self.linked_overlapp_histo[bin]+=1;
                }
                
            },
        };
    
    }

    pub fn pt_init(&mut self, original_grid_id: u32)->Result<(),String>{
        let K = self.pt_enviroment.config.num_T_steps;
        let Q = self.pt_enviroment.config.num_grids_equal_T;
        let T_start = self.pt_enviroment.config.T_start;
        let beta_start = 1.0/T_start;
        let T_end = self.pt_enviroment.config.T_end;
        let beta_end= 1.0 / T_end;
        //let delta_T: f64 = (T_end-T_start)/(K as f64 - 1.0);
        let delta_beta: f64 = (beta_end-beta_start)/(K as f64 - 1.0);

        let original_grid = self.grid(original_grid_id).unwrap();
        //let mut grid_config = self.default_grid_config.clone();

        let mut current_beta:f64 = beta_start;
        for k in 0..K{
            let current_T:f64 = 1.0/current_beta;
            let mut equal_T_Ids = equalTGridIds{
                T: current_T,
                was_switched: Vec::new(),
                current_pt_acceptance_prob: 0.0,
                ids: Vec::new()};       
            for q in 0..Q{
                let mut pt_id: u32;
                if k==0 && q == 0{
                    pt_id = original_grid_id;
                }else {
                    pt_id = self.clone_grid(original_grid_id)?;
                    self.grid_mut(pt_id).unwrap().set_T(current_T);
                }
                equal_T_Ids.ids.push(pt_id);
            }
            self.pt_enviroment.pt_ids.push(equal_T_Ids);
            current_beta += delta_beta;
        }
        Ok(())
    }

    fn pt_acceptance_probability(
        &mut self, T_lower:f64, T_higher:f64, energy_lower:f64, energy_higher:f64)
        ->f64{
        let mut acceptance_prob = f64::min(
            1.0,
            ((1.0/T_lower-1.0/T_higher)*(energy_lower-energy_higher)).exp());
        return acceptance_prob
    }

    pub fn pt_exchange(&mut self){
        if !self.running{
            return
        }

        if self.pt_enviroment.pt_ids.len() == 0{
            return
        }

        let K = self.pt_enviroment.config.num_T_steps;
        let Q = self.pt_enviroment.config.num_grids_equal_T;
        let T_start = self.pt_enviroment.config.T_start;
        let beta_start= 1.0/T_start;
        let T_end = self.pt_enviroment.config.T_end;
        let beta_end = 1.0 / T_end;
        //let delta_T: f64 = (T_end-T_start)/(K as f64 - 1.0);
        let delta_beta: f64 = (beta_end-beta_start)/(K as f64 - 1.0);

        let mut rng = thread_rng();

        let mut was_switched:Vec<bool> = Vec::new();

        
        let mut current_beta:f64 = beta_start;
        for k in 0..K-1{
            let lower_T:f64 = 1.0/(current_beta+delta_beta);
            let higher_T = 1.0/current_beta;
            /* let mut lower_equalTGrids = self.pt_enviroment.pt_ids.iter_mut().find(|e|e.T == current_T).unwrap();
            let mut higher_equalTGrids = self.pt_enviroment.pt_ids.iter_mut().find(|e|e.T == current_T+delta_T).unwrap();
             */
            let mut lower_ids = self.pt_enviroment.pt_ids.iter().find(|e|e.T == lower_T).unwrap().ids.clone();
            let mut higher_ids = self.pt_enviroment.pt_ids.iter().find(|e|e.T == higher_T).unwrap().ids.clone();
            //self.pt_enviroment.pt_ids.retain(|e|(e.T != current_T) && (e.T != current_T+delta_T));



            for _ in 0..lower_ids.len(){
                let lower_id = lower_ids.iter().choose(&mut rng).unwrap().clone();
                let higher_id = higher_ids.iter().choose(&mut rng).unwrap().clone();

                let mut lower_energy = self.grid(lower_id).unwrap().calc_energy();
                let mut higher_energy = self.grid(higher_id).unwrap().calc_energy();

                let acceptance_prop = self.pt_acceptance_probability(
                    lower_T, higher_T, lower_energy, higher_energy);

                if rand::thread_rng().gen_bool(acceptance_prop){
                    self.grid_mut(lower_id).unwrap().set_T(higher_T);
                    self.grid_mut(higher_id).unwrap().set_T(lower_T);
                    lower_ids.retain(|&e|e != lower_id);
                    lower_ids.push(higher_id);
                    higher_ids.retain(|&e|e != higher_id);
                    higher_ids.push(lower_id);
                    was_switched.push(true);
                }else{
                    was_switched.push(false);
                }

            }

            self.pt_enviroment.pt_ids.iter_mut().find(|e|e.T == lower_T)
                .unwrap().was_switched.append(&mut was_switched.clone());
            self.pt_enviroment.pt_ids.iter_mut().find(|e|e.T == higher_T)
                .unwrap().was_switched.append(&mut was_switched);
            
            self.pt_enviroment.pt_ids.iter_mut().find(|e|e.T == lower_T)
                .unwrap().update_current_acceptance_prob();
            self.pt_enviroment.pt_ids.iter_mut().find(|e|e.T == higher_T)
                .unwrap().update_current_acceptance_prob();


            self.pt_enviroment.pt_ids.iter_mut().find(|e|e.T == lower_T)
                .unwrap().ids = lower_ids;
            self.pt_enviroment.pt_ids.iter_mut().find(|e|e.T == higher_T)
                .unwrap().ids = higher_ids;


            //self.pt_enviroment.pt_ids.push(equalTGridIds { T: current_T, was_switched: vec![true,true], ids: lower_ids });
            //self.pt_enviroment.pt_ids.push(equalTGridIds { T: current_T+delta_T, was_switched: vec![true,true], ids: higher_ids });
            current_beta += delta_beta;


        }

    }

    pub fn new_grid(&mut self)->Result<u32,String>{
        let new_id = self.find_free_grid_id()?;
        let grid = Grid::new(self.default_grid_config.clone(), new_id)?;
        self.grids.push(grid);
        self.current_grid_ids.push(new_id);
        Ok(new_id)
    }

    pub fn clone_grid(&mut self,grid_id: u32)->Result<u32,String>{
        let mut orig_grid = match self.grid(grid_id){
            Some(grid) => grid,
            None => return Err("Grid id not found. Can't clone.".to_owned()),
        };
        let mut cloned_grid = orig_grid.clone();
        let new_id = self.find_free_grid_id()?;
        cloned_grid.set_id(new_id);
        cloned_grid.init_output_file();
        self.grids.push(cloned_grid);
        self.current_grid_ids.push(new_id);
        return Ok(new_id)
    }

    fn find_free_grid_id(&self)->Result<u32,String>{
        for new_id in 0..u32::MAX{
            match self.current_grid_ids.iter().find(|&&id|id==new_id){
                Some(_) => continue,
                None => return Ok(new_id),
            }
        }
        Err("Can't find free grid id.".to_owned())
    }

    pub fn delete_grid(&mut self, grid_id:u32){
        self.grids.retain(|grid| grid.id() != grid_id);
        self.current_grid_ids.retain(|&existing_id|existing_id != grid_id);
        self.pt_enviroment.delete_id(grid_id);
    }

    pub fn queue_grid_deletion(&mut self, grid_id:u32){
        self.grids_to_delete.push(grid_id);
    }

    pub fn queue_all_grids_deletion(&mut self){
        for grid_id in self.current_grid_ids.clone().iter(){
            self.queue_grid_deletion(*grid_id)
        }
    }

    pub fn delete_queded_grids(&mut self){
        while let Some(id) = self.grids_to_delete.pop(){
            self.delete_grid(id)
        };
    }

    pub fn current_grid_ids(&self)->&Vec<u32>{return &self.current_grid_ids}

    pub fn num_grids(&self)->u32{return self.grids.len() as u32}

    pub fn run(&mut self){
        for _ in 0..self.config.num_steps{
            self.thread_step();
        }
    }

    #[cfg(not(target_arch = "wasm32"))]
    pub fn thread_step(&mut self){
        if !self.running{
            return
        }
        //function to compute multiple metropolis steps for multiple grids using multithreading


        //communication channel to parent
        let (tx,rx) = mpsc::channel();

        //Loop over every existing grid
        for _ in 0..self.grids.len(){
            //pop grid to gain ownership
            let mut thread_grid = self.grids.pop().unwrap();
            //clone channel sender
            let thread_tx = tx.clone();
            //amount of metropolis steps to perform before thread closes
            let thread_steps = self.config.thread_steps.clone();

            thread::spawn(move ||{
                //run N metropolis steps in thread
                for _ in 0..thread_steps{
                    thread_grid.metropolis_step();
                }
                //calculate and save temperature and energy of grid
                thread_grid.update_history();
                //send results
                thread_tx.send(thread_grid).unwrap();
            });
        }
        //original Sender needs to be dropped. Every grid copied from this one, so tx isn't used.
        //if it isn't dropped we get a deadlock because rx waits for tx
        drop(tx);

        //gather updated grids
        for grid in rx{
            self.grids.push(grid);
        };


    }

    #[cfg(target_arch = "wasm32")]
    pub fn thread_step(&mut self){
        if !self.running{
            return
        }
        //function to compute multiple metropolis steps for multiple grids using multithreading

        //amount of metropolis steps to perform before thread closes
        let thread_steps = self.config.thread_steps.clone();
        for grid in self.grids.iter_mut(){
            for _ in 0..thread_steps{
                grid.metropolis_step();
            }
            grid.update_history();
        }
    }

    pub fn grid(&self, grid_id: u32)-> Option<&Grid>{
        for grid in self.grids.iter(){
            if grid_id == grid.id(){
                return Some(grid)
            }
        }
        return None
    }

    pub fn grid_mut(&mut self, grid_id: u32)-> Option<&mut Grid>{
        for grid in self.grids.iter_mut(){
            if grid_id == grid.id(){
                return Some(grid)
            }
        }
        return None
    }    


    pub fn history_by_id(&self, grid_id:u32)->Option<&History>{
        return Some(self.grid(grid_id)?.history());

    }

}