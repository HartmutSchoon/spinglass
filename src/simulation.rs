use std::collections::VecDeque;
use std::fs;
use std::sync::mpsc;
use std::sync::{Arc, Mutex};
use std::thread::{self, JoinHandle, current};
use std::path;
use std::time;

use rand::seq::IteratorRandom;
use rand::seq::SliceRandom;
use rand::{Rng, thread_rng};

use crate::config::{self, SimulationConfig, GridConfig, PTConfig};
use crate::grid::{Grid,History, Sweep};


#[derive(Clone)]
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

#[derive(Clone)]
pub struct equalTGridIds{
    pub T:f64,
    pub was_switched: VecDeque<bool>,
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

    pub fn add_single_switch(&mut self, switch: bool){
        if self.was_switched.len() >= 1000 {
            self.was_switched.pop_front();
        }
        self.was_switched.push_back(switch)
    }

    pub fn add_vec_switch(&mut self, switch_vec: Vec<bool>){
        for &elem in switch_vec.iter(){
            self.add_single_switch(elem);
        }
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
    pub fn new(global_config: config::Config)->Self{

        let config = global_config.simulation_config.clone();
        let default_grid_config = global_config.grid_config.clone();
        let grids = Vec::new();
        let current_grid_ids : Vec<u32> = Vec::new();
        let grids_to_delete: Vec<u32> = Vec::new();
        let running = false;
        let pt_enviroment = PTEnviroment::new(global_config.pt_config.clone());

        let num_bins = (2.0/config.histo_width+1.0).floor() as usize;
        let overlapp_histo = vec![0;num_bins];
        let linked_overlapp_histo = vec![0;num_bins];
        
        let sim =  Simulation{
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
        fs::create_dir_all(&self.default_grid_config.save_path)
            .expect(&format!("Can't create save path: {}", self.default_grid_config.save_path));


        //Look every second for 10 seconds if directory exists. If not panic!
/*         let mut loop_counter = 0;
        loop{
            match path::Path::new(&self.default_grid_config.save_path).exists(){
                true => break,
                false =>{
                    if loop_counter >= 10 {panic!("Can't create directory {}",self.default_grid_config.save_path)};
                    loop_counter += 1;
                    thread::sleep(time::Duration::from_secs(1));
                }

            }
        } */

        fs::copy("./config.toml", self.default_grid_config.save_path.clone() + "/config.toml") 
            .expect(&format!("Can't copy config.toml to save path: {}", self.default_grid_config.save_path));
        // let config_string = toml::to_string(&self.config).unwrap();
        // fs::File::create("./results/config.toml").unwrap().write_all(config_string.as_bytes());
    }

    #[cfg(target_arch = "wasm32")]
    fn init_dir(&self){
        return;
    }    

    pub fn simulation_step(&mut self){
        //self.thread_step();
        self.not_threaded_step();
        self.pt_exchange();
        self.update_all_grid_histories();
        self.update_overlapp_histogramm(0.5);
    }

    pub fn custom_run(&mut self){
        let grid_id1 = self.new_grid().unwrap();
        self.pt_init(grid_id1).unwrap();
        self.update_all_grid_histories();
        self.running = true;
    }

    pub fn update_all_grid_histories(&mut self){
        if self.running == false {return}
        for grid_id in self.current_grid_ids.clone().iter(){
            // we need to check if grid_id is in a pt_enviroment
            let grid_ref = self.grid(*grid_id).unwrap();
            let T = grid_ref.T();
            let capacity = grid_ref.capacity();
            let mut current_pt_ids= self.pt_enviroment.pt_ids.clone();
            //this one now only stores the one equalTIds where its T is equal to this grid T
            current_pt_ids.retain(|elem|elem.T == T);
            //now call first() to see if there is a element
            match current_pt_ids.first(){
                Some(equalTGriddIds) =>{
                    //We found a equalT struct with this T. Now we need to find our grid_id
                    let mut other_ids = equalTGriddIds.ids.clone();
                    let len = other_ids.len();
                    other_ids.retain(|&elem|elem != *grid_id);
                    //if other_ids didn't shrink in size the current grid_id wasn't in a pt_enviroment
                    //therefore there can't be a linked overlapp
                    //it would be possible to create a equalTGrid with only this grid in it
                    //Then it also wouldn't be possible to compute the linked overlapp
                    // then the size of other_ids would be 0
                    if len == other_ids.len() || other_ids.len() == 0{
                        self.grid_mut(*grid_id).unwrap().update_history(f64::NAN, f64::NAN)
                    }
                    
                    //sample some other id at the same temperature
                    let other_id = other_ids.choose(&mut rand::thread_rng()).unwrap();

                    //calc linked overlapp
                    let linked_overlapp = self.grid(*grid_id).unwrap().calc_linked_overlap(
                        self.grid(*other_id).unwrap()).unwrap();
    
                    //norm it by num of links
                    //z_half is the half of the number of links each particle has. This corresponds to
                    //the number of unique links per particle
                    //see :Monte Carlo simulations of spin glasses at low temperatures
                    //Helmut G. Katzgraber, Matteo Palassini, and A. P. Young
                    let z_half = self.grid(*grid_id).unwrap().dimensions().len() as f64;
                    let normed_linked_overlapp = linked_overlapp/(z_half*capacity as f64);

                    let overlap = self.grid(*grid_id).unwrap().calc_overlap(
                        self.grid(*other_id).unwrap()).unwrap();
                    let normed_overlap = overlap/capacity as f64;

                    //let grid pdate history with this overlapp
                    self.grid_mut(*grid_id).unwrap().update_history(normed_overlap,normed_linked_overlapp);

                },
                None => {
                    //we did not find a equalTGridIDs with matching T -> grid isnt in Pt enviroment
                    //->no overlapp can be calculated
                    self.grid_mut(*grid_id).unwrap().update_history(f64::NAN,f64::NAN);
                },
            }
        }
    }


    pub fn update_overlapp_histogramm(&mut self, target_T: f64){
        /* This function calculates the linked overlap for a target temperature
        between the first and every other grid existing at this temperature.
        It adds a normed linked overlap to self.linked_overlapp_histo, which ist 
        later used to draw a overlapp histogramm. */
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
        
        let Q = self.pt_enviroment.config.num_grids_equal_T;
        let T_start = self.pt_enviroment.config.T_start;
        //let beta_start = 1.0/T_start;
        let T_end = self.pt_enviroment.config.T_end;
        //let beta_end= 1.0 / T_end;
        //let delta_beta: f64 = (beta_end-beta_start)/(K as f64 - 1.0);


        let original_grid = self.grid(original_grid_id).unwrap();

        let mut current_T = T_start;
        //let mut current_beta:f64 = beta_start;
        loop{
        //for k in 0..K{
            //let current_T:f64 = 1.0/current_beta;
            if current_T> T_end{
                break
            }

            let delta_T = self.pt_lin_fun(current_T);
            //let delta_T = self.pt_logistic_fun(current_T);
            //let delta_T = (T_end-T_start)/(K as f64 -1.0);

            let mut equal_T_Ids = equalTGridIds{
                T: current_T,
                was_switched: VecDeque::new(),
                current_pt_acceptance_prob: 0.0,
                ids: Vec::new()};       
            for q in 0..Q{
                let mut pt_id: u32;
                if current_T == T_start && q == 0{
                    pt_id = original_grid_id;
                }else {
                    pt_id = self.clone_grid(original_grid_id)?;
                    self.grid_mut(pt_id).unwrap().set_T(current_T);
                    self.grid_mut(pt_id).unwrap().shuffle_spins();
                }
                equal_T_Ids.ids.push(pt_id);
            }
            self.pt_enviroment.pt_ids.push(equal_T_Ids);

            current_T += delta_T;
            //current_beta += delta_beta;
        }
        Ok(())
    }

    pub fn pt_exchange(&mut self){
        if !self.running{
            return
        }

        if self.pt_enviroment.pt_ids.len() == 0{
            return
        }

        let Q = self.pt_enviroment.config.num_grids_equal_T;
        let T_start = self.pt_enviroment.config.T_start;
        //let beta_start = 1.0/T_start;
        let T_end = self.pt_enviroment.config.T_end;
        //let beta_end= 1.0 / T_end;
        //let delta_beta: f64 = (beta_end-beta_start)/(K as f64 - 1.0);

        let mut rng = thread_rng();

        let mut current_T = T_start;
        //let mut current_beta:f64 = beta_start;
        loop{
        //for k in 0..K-1{
            let mut was_switched:Vec<bool> = Vec::new();

            let delta_T = self.pt_lin_fun(current_T);
            //let delta_T = self.pt_logistic_fun(current_T);
            //let delta_T = (T_end-T_start)/(K as f64 -1.0);
            let lower_T:f64 = current_T;
            let higher_T = current_T + delta_T;
            if higher_T > T_end{
                break
            }

            /* let lower_T:f64 = 1.0/current_beta;
            let higher_T = 1.0/(current_beta+delta_beta); */

            let mut lower_ids = self.pt_enviroment.pt_ids.iter().find(|e|e.T == lower_T).unwrap().ids.clone();
            let mut higher_ids = self.pt_enviroment.pt_ids.iter().find(|e|e.T == higher_T).unwrap().ids.clone();

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
                .unwrap().add_vec_switch(was_switched.clone());
            self.pt_enviroment.pt_ids.iter_mut().find(|e|e.T == higher_T)
                .unwrap().add_vec_switch(was_switched);
            
            self.pt_enviroment.pt_ids.iter_mut().find(|e|e.T == lower_T)
                .unwrap().update_current_acceptance_prob();
            self.pt_enviroment.pt_ids.iter_mut().find(|e|e.T == higher_T)
                .unwrap().update_current_acceptance_prob();


            self.pt_enviroment.pt_ids.iter_mut().find(|e|e.T == lower_T)
                .unwrap().ids = lower_ids;
            self.pt_enviroment.pt_ids.iter_mut().find(|e|e.T == higher_T)
                .unwrap().ids = higher_ids;


            //current_beta += delta_beta;
            current_T += delta_T;
        }

    }

    fn pt_acceptance_probability(
        &mut self, T_lower:f64, T_higher:f64, energy_lower:f64, energy_higher:f64)
        ->f64{
        let mut acceptance_prob = f64::min(
            1.0,
            ((1.0/T_lower-1.0/T_higher)*(energy_lower-energy_higher)).exp());
        return acceptance_prob
    }

    fn pt_lin_fun(& self, T:f64) -> f64{
        //let L = self.pt_enviroment.config.logic_L;
        let m = self.pt_enviroment.config.linear_m;
        let dT0 = self.pt_enviroment.config.linear_dT0;
        let dT = m*T+dT0;
        return dT
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
            let steps_per_sweep = self.config.steps_per_sweep.clone();

            thread::spawn(move ||{
                //run N metropolis steps in thread
                for _ in 0..steps_per_sweep{
                    thread_grid.metropolis_step();
                }
                //calculate and save temperature and energy of grid
                //thread_grid.update_history();
                //send results
                thread_tx.send(thread_grid).unwrap();
            });
        }
        //original sender needs to be dropped. Every grid copied from this one, so tx isn't used.
        //if it isn't dropped we get a deadlock because rx waits for tx
        drop(tx);

        //gather updated grids
        for grid in rx{
            self.grids.push(grid);
        };


    }

    pub fn not_threaded_step(&mut self){
        for grid in self.grids.iter_mut(){
            for _ in 0..self.config.steps_per_sweep{
                grid.metropolis_step();
            }
        }
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



