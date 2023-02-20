use crate::{config::*, particles::*};
use log::{error, warn};
use rand::{self, Rng};
use rand_distr::Distribution;

use std::fs;
use std::thread;
use std::time;
use std::path;

#[derive(Clone)]
pub struct History {
    pub run: VecDeque<u32>,
    pub T: VecDeque<f64>,
    pub energy: VecDeque<f64>,
    pub magnetization: VecDeque<i32>,
    pub linked_overlapp: VecDeque<f64>,
    /*These two are implementations of the equlibrium condition of this paper:
    Monte Carlo simulations of spin glasses at low temperatures
    Helmut G. Katzgraber, Matteo Palassini, and A. P. Young*/
    pub katz_energy: VecDeque<f64>,
    pub av_linked_overlapp: VecDeque<f64>,
    pub av_katz_energy: VecDeque<f64>,
    current_size: usize,
    capacity: usize,
}
impl History {
    pub fn new(capacity: usize) -> Self {
        return History {
            run: VecDeque::new(),
            T: VecDeque::new(),
            energy: VecDeque::new(),
            magnetization: VecDeque::new(),
            linked_overlapp: VecDeque::new(),
            katz_energy: VecDeque::new(),
            av_linked_overlapp: VecDeque::new(),
            av_katz_energy: VecDeque::new(),
            current_size: 0,
            capacity,

        };
    }
    pub fn add(&mut self,
        run: u32,
        T: f64,
        energy: f64,
        magnetization: i32,
        linked_overlapp: f64,
        katz_energy: f64,
        av_linked_overlapp: f64,
        av_katz_energy: f64) {
        if self.capacity == 0 {
            return;
        }
        if self.current_size == self.capacity {
            self.remove_oldest()
        }

        self.run.push_back(run);
        self.T.push_back(T);
        self.energy.push_back(energy);
        self.magnetization.push_back(magnetization);
        self.linked_overlapp.push_back(linked_overlapp);
        self.katz_energy.push_back(katz_energy);
        self.av_linked_overlapp.push_back(av_linked_overlapp);
        self.av_katz_energy.push_back(av_katz_energy);

        self.current_size += 1;
    }
    pub fn remove_oldest(&mut self) {
        if self.current_size == 0 {
            return;
        }

        self.run.pop_front();
        self.T.pop_front();
        self.energy.pop_front();
        self.magnetization.pop_front();
        self.linked_overlapp.pop_front();
        self.katz_energy.pop_front();
        self.av_linked_overlapp.pop_front();
        self.av_katz_energy.pop_front();

        self.current_size -= 1;
    }

    pub fn current_size(&self) -> usize {
        return self.current_size;
    }
    pub fn capacity(&self) -> usize {
        return self.capacity;
    }
}

#[derive(Clone)]
pub struct Sweep {
    T_start: f64,
    T_end: f64,
    num_steps: u64,
    delta_T: f64,
}

impl Sweep {
    pub fn new(T_start: f64, T_end: f64, num_steps: u64) -> Self {
        let delta_T = (T_end - T_start) / (num_steps as f64);
        return Sweep {
            T_start,
            T_end,
            num_steps,
            delta_T,
        };
    }
    pub fn update_T(&mut self) {
        self.delta_T = (self.T_end - self.T_start) / (self.num_steps as f64);
    }
    pub fn set_T_start(&mut self, T_start: f64) {
        self.T_start = T_start;
        self.update_T()
    }
    pub fn set_T_end(&mut self, T_end: f64) {
        self.T_end = T_end;
        self.update_T()
    }
    pub fn set_num_steps(&mut self, num_steps: u64) {
        self.num_steps = num_steps;
        self.update_T()
    }

    pub fn T_start(&self) -> f64 {
        return self.T_start;
    }
    pub fn T_end(&self) -> f64 {
        return self.T_end;
    }
    pub fn num_steps(&self) -> u64 {
        return self.num_steps;
    }
}

impl Default for Sweep {
    fn default() -> Self {
        return Sweep::new(5.0, 0.0, 10000000);
    }
}

#[derive(Clone)]
pub struct GridPos {
    id: u32,
    pos: Vec<u16>,
    particle: Option<Particle>,
}

impl GridPos {
    pub fn id(&self) -> u32 {
        return self.id;
    }
    pub fn pos(&self) -> &Vec<u16> {
        return &self.pos;
    }
    pub fn particle(&self) -> Option<&Particle> {
        return self.particle.as_ref();
    }
    pub fn particle_as_mut(&mut self) -> Option<&mut Particle> {
        return self.particle.as_mut();
    }
}

//#[derive(Clone)]
pub struct Grid {
    id: u32,
    run: u32,
    pub cloned_from: Option<u32>,
    pub config: GridConfig,
    dimensions: Vec<u16>,
    pub capacity: u32,
    num_particles: u32,
    T: f64,
    external_field: f64,
    //TODO: warum den hier eigl nicht als hashmap? macht mehr sinn wenn ich indexe raus nehmen will
    grid_positions: Vec<GridPos>,
    history: History,
    sweep: Option<Sweep>,
    output_file: Option<fs::File>,
    //rng_thread: ThreadRng,
    //rng_grid_idx: Uniform<usize>
}

impl Grid {
    pub fn new(config: GridConfig, id: u32) -> Result<Self, String> {
        let run = 0;
        let cloned_from = None;
        // pre compute fields
        let dimensions = config.grid_dimensions.clone();
        let capacity: u32 = dimensions.iter().map(|&e| u32::from(e)).product();
        let num_particles = config.num_particles;
        let T = config.T;
        //let mut rng_thread = rand::thread_rng();
        //let rng_grid_idx = Uniform::from(
        //    0..usize::try_from(capacity).unwrap());
        let mut external_field = config.external_field;
        let mut grid_positions: Vec<GridPos> = Vec::new();
        let mut energy_history: Vec<f64> = Vec::new();
        let mut temp_history: Vec<f64> = Vec::new();
        let mut history = History::new(config.history_capacity);
        let mut sweep = None;
        let output_file = None;
        // create grid
        let mut grid = Grid {
            //TODO: warum unterscheiden zwischen sachen in Config und in Grid struct? Warum nicht z.B. alles in Config?
            id,
            run,
            cloned_from,
            config,
            dimensions,
            capacity,
            num_particles,
            grid_positions,
            T,
            external_field,
            history,
            sweep,
            output_file,
            //rng_thread,
            //rng_grid_idx,
        };
        grid.create_grid_positions()
            .init_particles()
            .link_all_particles()
            .init_output_file();
        //grid.init_particles(grid.config.num_particles).link_all_particles();
        //grid.update_history(0.0);
        return Ok(grid);
    }

    #[cfg(not(target_arch = "wasm32"))]
    pub fn init_output_file(&mut self) -> &mut Self {
        let path = format!("{}/grid_{}.tsv",self.config.save_path, self.id);  

        fs::remove_file(&path);
/*         let mut loop_counter = 0;
        loop{
            match path::Path::new(&path).exists(){
                false => break,
                true=>{
                    if loop_counter >= 10 {panic!("Cant' create file handle for grid! Path: {}", &path)};
                    loop_counter += 1;
                    thread::sleep(time::Duration::from_secs(1));
                }
            };
        }; */
        
        self.output_file = Some(
            fs::OpenOptions::new()
                .create_new(true)
                .append(true)
                .open(path)
                .expect("Unable to open file"),
        );
        let header = String::from("run\tT\tenergy\tmagnetization\tlinked_overlapp\n");

        self.output_file
            .as_mut()
            .unwrap()
            .write_all(header.as_bytes());

        return self;
    }

    #[cfg(target_arch = "wasm32")]
    pub fn init_output_file(&self) {
        return;
    }

    pub fn write_to_output_file(&mut self,
        run: u32,
        T: f64,
        energy: f64,
        magnetization: i32,
        linked_overlapp:f64) {
        if let Some(file) = self.output_file.as_mut() {
            let output = String::from(format!("{}\t{}\t{}\t{}\t{}\n",
                 run, T, energy, magnetization,linked_overlapp));
            file.write_all(output.as_bytes());
        }
    }

    fn create_grid_positions(&mut self) -> &mut Self {
        //position
        let mut curr_pos: Vec<u16> = vec![0; self.dimensions.len()];
        let mut old_pos = curr_pos.clone();

        //push last GridPos by hand
        self.grid_positions.push(GridPos {
            id: 0,
            pos: curr_pos,
            particle: None,
        });

        //Loop over all possible positions, without last one
        for id in 1..self.capacity {
            //Update Position to new one
            curr_pos = self.next_pos(&old_pos).unwrap();
            old_pos = curr_pos.clone();
            //create empty GridPos and push to grid
            self.grid_positions.push(GridPos {
                id,
                pos: curr_pos,
                particle: None,
            });
        }
        return self;
    }

    fn init_particles(&mut self) -> &mut Self {
        match self.config.particle_position_method {
            ParticlePositionMethod::RandomParticlePositions => self.random_init(),
            ParticlePositionMethod::FillGrid => self.grid_fill_init(),
        };
        return self;
    }

    fn random_init(&mut self) -> &mut Self {
        let num_particles = self.config.num_particles;
        if num_particles >= self.capacity {
            warn!(target: "grid","Tried to randomly fill more particles to grid than it has capacity.
                Filled grid with one particle at each position");
            return self.grid_fill_init();
        }

        const MAX_ITER: u32 = 10000;
        let mut iter = 0;
        let mut particle_counter = 0;
        while particle_counter < num_particles {
            if iter >= MAX_ITER {
                //TODO: Handle this one graceful!
                panic!("In random particle init: Can't find location for new particle!");
            }

            let rnd_idx = rand::thread_rng().gen_range(0..self.capacity) as usize;
            let particle = &mut self.grid_positions[rnd_idx].particle;
            if particle.is_none() {
                let spin = match self.config.spin_up_init {
                    true => Spin::Up,
                    //Distribution<Spin> is implemented, therefore rng can sample Spins
                    false => rand::thread_rng().gen(),
                };
                *particle = Some(Particle::new(rnd_idx.try_into().unwrap(), spin));
                particle_counter += 1;
                iter = 0;
            }
            iter += 1;
        }
        return self;
    }
    fn grid_fill_init(&mut self) -> &mut Self {
        for grid_pos in self.grid_positions.iter_mut() {
            let spin = match self.config.spin_up_init {
                true => Spin::Up,
                //Distribution<Spin> is implemented, therefore rng can sample Spins
                false => rand::thread_rng().gen(),
            };
            grid_pos.particle = Some(Particle::new(grid_pos.id, spin));
        }
        self.num_particles = self.capacity;
        return self;
    }

    //TODO: irgendwie nicht so intuitiv. Das ding füllt explizit "von ob nach
    //unten". Sollte die funktion umbenennen und weitere grid fill methoden
    //einbauen
    /* fn init_particles2(&mut self, num:u32)->Result<&mut Self,String>{
        //function to set "num" particles on the grid by filling one dimension
        //after another.
        //eg: (idx=0):[0,0] -> (idx=1):[1,0] -> (idx=2):[0,1] -> (idx=3):[1,1]
        //for 2x2 grid and num=4.
        //
        //Input:    num                       amount of particles to set on grid
        //Output:   Result Ok(self)           init was succesful, return
        //                                    self in Ok()
        //          Result Err(String)        Error occured, return error
        //                                    message in Err()

        //check if "num" particles fit on grid
        if u32::try_from(num).unwrap() > self.capacity{
            return Err(String::from(format!("Can't fit {} particles on grid with capacity of {}!", num, self.capacity)));
        }

        //position of new particle
        let mut old_pos: Vec<u16> = vec![0;self.dimensions.len()];
        //position of new particle
        let mut curr_pos: Vec<u16> = vec![0;self.dimensions.len()];

        let mut spin = match self.config.spin_up_init{
            true => Spin::Up,
            false => rand::thread_rng().gen(),
        };

        //Init first particle
        self.grid_positions[0]= Some(
            Particle::new(
                0,
                spin,
                curr_pos));

        //Loop through all positions
        for idx in 1..self.capacity{
            curr_pos = self.next_pos(&old_pos)?;
            old_pos = curr_pos.clone();

            let conv_idx:usize = idx.try_into().unwrap();
            if idx < self.config.num_particles{
                spin = match self.config.spin_up_init{
                    true => Spin::Up,
                    false => rand::thread_rng().gen(),
                };
                let conv_idx:usize = idx.try_into().unwrap();
                self.grid_positions[conv_idx] = Some(
                    Particle::new(conv_idx,spin,curr_pos));
            }else{
                self.grid_positions[conv_idx] = None;
            }

        }
        return Ok(self);

    } */

    pub fn link_all_particles(&mut self) -> &mut Self {
        //function to "link" all existing particles on grid to their neighbours
        // by adding the indices of the neighbour positions to the particles
        // neighbour field.

        let mut neighbour_indices: Vec<usize>;

        for id in 0..self.capacity {
            let grid_pos = &self.grid_positions[id as usize];
            match grid_pos.particle() {
                None => continue,
                Some(_) => {
                    let neighbour_ids = self.calc_filled_neighbour_indices(id);
                    for neighbour_id in neighbour_ids {
                        self.link_particles(id, neighbour_id);
                    }
                }
            }
        }
        return self;
    }

    pub fn link_particles(&mut self, grid_id1: u32, grid_id2: u32) -> Result<&mut Self, String> {
        //function which takes two grid positions and tries to link both particles to each other.
        //Output: Ok(Self)      both particles have been linked. It is garuanteed that both links are
        //                      the same. If one particle was already linked to the other the link will
        //                      will be copied.
        //        Err(String)   functions returnes Err() if one of the grid positions wasn't filled with
        //                      a particle
        if self.grid_positions[grid_id1 as usize].particle.is_none() {
            return Err("No particle found at position 1".to_owned());
        }

        if self.grid_positions[grid_id2 as usize].particle.is_none() {
            return Err("No particle found at position 2.".to_owned());
        }

        let link_from1_to2 = self.grid_positions[grid_id1 as usize]
            .particle()
            .unwrap()
            .find_link(grid_id2);
        let link_from2_to1 = self.grid_positions[grid_id2 as usize]
            .particle()
            .unwrap()
            .find_link(grid_id1);

        //if both particles have links to each other
        if link_from1_to2.is_some() && link_from2_to1.is_some() {
            //If links aren't the same
            if link_from1_to2.unwrap() != link_from2_to1.unwrap() {
                error!(target:"grid","Found particles with mismatched couplings!");
            }
            return Ok(self);
        }

        //if particle 1 has link to particle 2 but particle 2 hasn't to particle 1
        if link_from1_to2.is_some() && link_from2_to1.is_none() {
            //link particle 2 to 1 with coupling from link from particle 1 to 2
            self.grid_positions[grid_id2 as usize]
                .particle_as_mut()
                .unwrap()
                .link_neighbour(Link::new_from(grid_id1, &link_from1_to2.unwrap()));
            return Ok(self);
        }

        //if particle 2 has link to particle 1 but particle 1 hasn't to particle 1
        if link_from2_to1.is_some() && link_from1_to2.is_none() {
            //link particle 1 to 2 with coupling from link from particle 2 to 1
            self.grid_positions[grid_id1 as usize]
                .particle_as_mut()
                .unwrap()
                .link_neighbour(Link::new_from(grid_id2, &link_from2_to1.unwrap()));
            return Ok(self);
        }

        // At this point both particles don't have links to each other

        //check if both limits from config are the same. If so don't use rng
        //and set coupling to constant
        /* let coupling = match self.config.coupling_limits[0] == self.config.coupling_limits[1] {
            true => self.config.coupling_limits[0],
            false => rand::thread_rng()
                .gen_range(self.config.coupling_limits[0]..self.config.coupling_limits[1]),
        }; */

        let mean = self.config.coupling_mean;
        let std_dev = self.config.coupling_variance.sqrt();
        
        let coupling = rand_distr::Normal::new(mean, std_dev).unwrap().sample(&mut rand::thread_rng());

        //create link
        let link_from1_to2 = Link::new(grid_id2, coupling);
        //copy link from first link
        let link_from2_to1 = Link::new_from(grid_id1, &link_from1_to2);
        //set links
        self.grid_positions[grid_id1 as usize]
            .particle_as_mut()
            .unwrap()
            .link_neighbour(link_from1_to2);
        self.grid_positions[grid_id2 as usize]
            .particle_as_mut()
            .unwrap()
            .link_neighbour(link_from2_to1);
        return Ok(self);
    }

    pub fn calc_overlap(&self, other: &Self) -> Result<f64,String> {
        if self.dimensions.iter().zip(other.dimensions.iter()).filter(|&(a,b)| a!=b).count() > 0{
            return Err("Grid dimensions mismatch!".to_owned())
        }

        let mut overlap = 0.0;
        for own_position in self.grid_positions.iter() {
            if own_position.particle.is_none() {
                continue;
            }

            let other_position = other.find_grid_position(own_position.id()).unwrap();
            if other_position.particle.is_none() {
                continue;
            };

            let own_spin: f64 = own_position.particle.as_ref().unwrap().spin().into();
            let other_spin: f64 = other_position.particle.as_ref().unwrap().spin().into();

            overlap += own_spin * other_spin;
        }
        return Ok(overlap);
    }

    pub fn calc_linked_overlap(&self, other: &Self) -> Result<f64,String> {
        if self.dimensions.iter().zip(other.dimensions.iter()).filter(|&(a,b)| a!=b).count() > 0{
            return Err("Grid dimensions mismatch!".to_owned())
        }

        let mut linked_overlap = 0.0;
        for own_position in self.grid_positions.iter() {
            if own_position.particle.is_none() {
                continue;
            }

            let other_position = other.find_grid_position(own_position.id()).unwrap();
            if other_position.particle.is_none() {
                continue;
            };

            let own_spin: f64 = own_position.particle.as_ref().unwrap().spin().into();
            let other_spin: f64 = other_position.particle.as_ref().unwrap().spin().into();

            let neighbours = own_position.particle.as_ref().unwrap().neighbours();
            for neighbour in neighbours.iter(){
                let position_id =  neighbour.target_id() as usize;
                let own_neighbour_spin=f64::from(self.grid_positions[position_id].particle.as_ref().unwrap().spin());
                let other_neighbour_spin = match other.grid_positions[position_id].particle.as_ref(){
                    Some(particle) => f64::from(particle.spin()),
                    None => 0.0,
                };
                linked_overlap += 0.5 * own_spin * own_neighbour_spin * other_spin * other_neighbour_spin;
            } 
        }
        return Ok(linked_overlap);
    }

    pub fn calc_energy(&self) -> f64 {
        //calculate energy/amount of particles of whole grid.
        //Returns: Energy

        //spin of selected particle
        let mut spin: &Spin;
        //neighbours of this particle
        let mut neighbour_links: &Vec<Link>;
        //idx and spin of those neighbours
        let mut neighbour_id: u32;
        let mut neighbour_spin: &Spin;
        //resulting energy
        let mut energy = 0.0;

        //loop through every position
        for position in self.grid_positions.iter() {
            //if there is particle at this position
            match position.particle() {
                None => continue,
                Some(particle) => {
                    //get refs to its spin and neighbours
                    spin = particle.spin();
                    neighbour_links = particle.neighbours();
                }
            }
            //loop through neighbour positions
            for neighbour_link in neighbour_links {
                neighbour_id = neighbour_link.target_id();
                //if there is particle at neighbour position
                match self.grid_positions[neighbour_id as usize].particle.as_ref() {
                    None => continue,
                    //get ref to its spin
                    Some(neighbour) => neighbour_spin = neighbour.spin(),
                }
                //Add energy contribution of this neighbour
                //0.5 because every link gets counted twice
                energy -=
                    0.5 * neighbour_link.coupling() * f64::from(spin) * f64::from(neighbour_spin);
            }
            //Add external field to energy
            energy -= 0.5 * self.external_field * f64::from(spin);
        }
        return energy;
    }

    pub fn calc_magnetization(&self) -> i32 {
        let mut spin_sum:i32 = 0;
        for position in self.grid_positions.iter() {
            if let Some(particle) = position.particle.as_ref() {
                spin_sum += i32::from(particle.spin());
            }
        }
        return spin_sum;
    }

    pub fn calc_av_linked_overlapp(&self)->f64{
        //this function computes the average over to "last half" of all data value of the linked overlapp.
        let start_idx = ((self.history.current_size as f64)/2.0) as usize;
        let end_idx = self.history.current_size;
        let mut sum = 0.0;
        for idx in start_idx..end_idx{
            sum += self.history.linked_overlapp[idx];
        }
        return sum/(end_idx-start_idx) as f64;
    }

    pub fn calc_av_energy(&self)->f64{
        //this function computes the average over to "last half" of all data value of the linked overlapp.
        let start_idx = ((self.history.current_size as f64)/2.0) as usize;
        let end_idx = self.history.current_size;
        let mut sum = 0.0;
        for idx in start_idx..end_idx{
            sum += self.history.energy[idx];
        }
        return sum/(end_idx-start_idx) as f64;
    }

    pub fn calc_av_katz(&self)->f64{
        //this function computes the average over to "last half" of all data value of the linked overlapp.
        let start_idx = ((self.history.current_size as f64)/2.0) as usize;
        let end_idx = self.history.current_size;
        let mut sum = 0.0;
        for idx in start_idx..end_idx{
            sum += self.history.katz_energy[idx];
        }
        return sum/(end_idx-start_idx) as f64;
    }

    pub fn metropolis_step(&mut self) -> Option<usize> {
        //perform metropolis_step on grid

        //energy difference
        let mut dH: f64 = 0.0;
        //spin of selected particle
        let mut spin: &Spin;
        //links to neighbours of selected particle
        let mut neighbour_links: &Vec<Link>;
        //idx and spin of a neighbour
        let mut neighbour_id: u32;
        let mut neighbour_spin: &Spin;

        //update run
        self.run += 1;

        //select random grid position
        //let grid_idx = self.rng_grid_idx.sample(&mut self.rng_thread);
        let grid_id = rand::thread_rng().gen_range(0..usize::try_from(self.capacity).unwrap());

        //is there some particle at this position?
        match self.grid_positions[grid_id].particle.as_ref() {
            None => return None,
            Some(particle) => {
                //get its spin and neighbours
                spin = particle.spin();
                neighbour_links = particle.neighbours();
            }
        };

        //Loop through neighbour positions
        for neighbour_link in neighbour_links {
            neighbour_id = neighbour_link.target_id();
            // is there some particle at neighbour position?
            match self.grid_positions[neighbour_id as usize].particle.as_ref() {
                None => continue,
                //get its spin
                Some(neighbour) => neighbour_spin = neighbour.spin(),
            }
            //Add contribution of this neighbour to energy difference
            dH += 2.0 * neighbour_link.coupling() * f64::from(spin) * f64::from(neighbour_spin);
        }
        //add external field to energie difference
        dH += 2.0 * self.external_field * f64::from(spin);

        // calculate metropolis acceptance probability
        let metro_p = f64::min(1.0, (-dH / self.T).exp());

        //Change Temperature if there should be a sweep
        if let Some(sweep) = &self.sweep {
            //check if sweep should end
            let sign = sweep.delta_T.signum();
            // the case where T grows towards T_end
            if (sign > 0.0) && (self.T >= sweep.T_end) ||
            //The case where T shrinks towards T_end
            (sign < 0.0) && (self.T <= sweep.T_end)
            {
                self.T = sweep.T_end;
                self.sweep = None;
            } else {
                self.T += sweep.delta_T;
            }
        };

        //random spin flip with metro_p and return
        match rand::thread_rng().gen_bool(metro_p) {
            false => return None,
            true => {
                self.grid_positions[grid_id]
                    .particle
                    .as_mut()
                    .unwrap()
                    .flip_spin();
                return Some(grid_id);
            }
        }
    }

    fn calc_filled_neighbour_indices(&self, grid_id: u32) -> Vec<u32> {
        // function to calculate all indeces from all neighbours of the given
        //particle
        //Input:    particle_idx   Particle whose neighbour indeces should be
        //                         calulated
        //Output:                   Vector of indeces

        let neighbour_positions = self.calc_neighbour_positions(grid_id);
        let mut neighbour_ids: Vec<u32> = Vec::new();
        for neighbour_pos in neighbour_positions {
            let neighbour_id = self.pos_to_id(&neighbour_pos).unwrap();
            match self.grid_positions[neighbour_id as usize].particle {
                Some(_) => neighbour_ids.push(neighbour_id),
                None => continue,
            }
        }
        return neighbour_ids;
    }

    fn calc_neighbour_positions(&self, grid_id: u32) -> Vec<Vec<u16>> {
        //function for calculating and returning the neighbouring positions to
        //a given position of a particle.
        //if position is on the edge of the grid, neighbours returned are on opposite side
        // eg: 3x3 grid, pos [0,0]: [[2,0],[1,0],[0,2],[0,1]]
        //
        //Input: particle   particle on grid for which the neighbours should be calculated
        //Output:           Vector of vectors containing neighbour positions

        let pos = self.grid_positions[grid_id as usize].pos();

        //Storage for neighbours
        let mut neighbour_positions = Vec::new();
        for axis in 0..self.dimensions.len() {
            //for every axis of position "pos" neighbours are at pos[axis]+/-1

            //store current posistion
            let mut smaller_neighbour = pos.clone();
            let mut bigger_neighbour = pos.clone();

            // try substracting 1 without leaving the grid, otherwise go to opposite side
            if smaller_neighbour[axis] <= 0 {
                smaller_neighbour[axis] = self.dimensions[axis] - 1;
            } else {
                smaller_neighbour[axis] -= 1;
            }

            // try adding 1 without leaving the grid, otherwise go to opposite side
            if bigger_neighbour[axis] >= self.dimensions[axis] - 1 {
                bigger_neighbour[axis] = 0;
            } else {
                bigger_neighbour[axis] += 1;
            }

            // push neighbours to solution Vector
            neighbour_positions.push(smaller_neighbour);
            neighbour_positions.push(bigger_neighbour);
        }
        return neighbour_positions;
    }

    fn is_valid_grid_position(&self, pos: &Vec<u16>) -> Result<(), String> {
        //function to check if given position pos is a valid position on grid
        //if dimensions of input position and grid don't match
        //
        //Input: pos    position to check

        if pos.len() != self.dimensions.len() {
            return Err("position dimensions don't match grid dimensions".to_owned());
        }
        //if one element of input position is larger than grid limits
        for axis in 0..self.dimensions.len() {
            if pos[axis] >= self.dimensions[axis] {
                return Err("Position can't fit into grid".to_owned());
            }
        }
        return Ok(());
    }

    fn idx_to_pos(&self, idx: u32) -> Result<Vec<u16>, String> {
        //function to convert from an index of a particle vector and given grid dimensions to the position of a particle.
        //particles are initialized by filling one dimension after another:
        //(idx=0):[0,0] -> (idx=1):[1,0] -> (idx=2):[0,1] -> (idx=3):[1,1] for 2x2 grid.
        //
        //Input: idx    index of particle vector

        //check if index is to large
        if idx >= self.capacity {
            return Err(String::from("Index to large for grid dimension"));
        }

        //init position vector with matching dimensions
        let mut pos: Vec<u16> = vec![0; self.dimensions.len()];

        //Try casting the index into "pos" type and store it as first entry
        match idx.try_into() {
            Ok(val) => pos[0] = val,
            Err(_) => {
                return Err(String::from("Index to large. Can't cast to u32."));
            }
        }

        // Main Loop: Loop through every dimension. If this dimension is to "full" (according to the grid size), push multiples of
        // current dimension size limit into next dimension
        // e.g 5x5 Grid, idx = 11: pos=[11,0]->pos=[11%5,11/5]=[1,2]
        for idx in 0..pos.len() - 1 {
            pos[idx + 1] = pos[idx] / self.dimensions[idx];
            pos[idx] = pos[idx] % self.dimensions[idx];
        }
        return Ok(pos);
    }

    fn pos_to_id(&self, pos: &Vec<u16>) -> Result<u32, String> {
        //function to convert from a position and given grid dimensions to the index of a particle.
        //particles are initialized by filling one dimension after another:
        //(idx=0):[0,0] -> (idx=1):[1,0] -> (idx=2):[0,1] -> (idx=3):[1,1] for 2x2 grid.
        //
        //Input: pos    position to convert

        //check if given position is valid
        self.is_valid_grid_position(pos)?;

        // Resulting index
        let mut id: u32 = 0;
        //Loop through every dimension
        for axis in 0..self.dimensions.len() {
            // For every axis we need to multiply the dimension limits of the previous dimensions
            //For example we have [4,3,2] grid and pos of [_,_,2]: we had to assign 2x(4x3) full
            //2D matrizes of particles to get to a position containing [_,_,2]

            let mut mult = 1;
            for prev_axis in 0..axis {
                mult *= self.dimensions[prev_axis];
            }
            id += u32::from(pos[axis]) * u32::from(mult);
        }

        return Ok(id);
    }

    fn next_pos(&self, curr_pos: &Vec<u16>) -> Result<Vec<u16>, String> {
        //returnes the "next" position to the current position or an error if current position is last position on grid.
        //For ex in given 2x2 grid with positions:
        //  {0,0} {0,1}
        //  {0,1} {1,1}
        //
        // the function returns:
        // current position: [0,0] Returns: Ok([1,0])
        // current position: [1,0] Returns: Ok([0,1])
        // current position: [0,1] Returns: Ok([1,1])
        // current position: [1,1] Returns: Err(msg)
        //
        // Input: curr_pos      current position
        // Output:  "next" position

        //check if given position is valid
        self.is_valid_grid_position(curr_pos)?;

        let axis = 0;
        return next_pos_recursion(curr_pos, &self.dimensions, axis);

        fn next_pos_recursion(
            curr_pos: &Vec<u16>,
            grid_dimensions: &Vec<u16>,
            axis: usize,
        ) -> Result<Vec<u16>, String> {
            // function body. Finds the next position by recursion.
            // function recursion is moved to own function to be able to set axis=0 default value.
            //
            //Input: grid_dimensions  stores size of grid
            //       axis              used for recursion, corresponding to a current dimension beeing looked at

            //stopping condition: if trying to set next position in dimension not available
            if axis >= grid_dimensions.len() {
                return Err(String::from("No \"next\" position available on grid."));
            }

            let mut next_pos = curr_pos.clone();

            //If there is room in this dimension add 1
            if next_pos[axis] + 1 < grid_dimensions[axis] {
                next_pos[axis] += 1;
                return Ok(next_pos);

            //Else try using next dimension
            } else {
                next_pos[axis] = 0;
                match next_pos_recursion(&next_pos, grid_dimensions, axis + 1) {
                    Ok(pos) => return Ok(pos),
                    Err(msg) => Err(msg),
                }
            }
        }
    }

    pub fn find_grid_position(&self, target_id: u32) -> Result<&GridPos, String> {
        if self.grid_positions[target_id as usize].id() == target_id {
            return Ok(&self.grid_positions[target_id as usize]);
        }
        let matching_positions: Vec<&GridPos> = self
            .grid_positions
            .iter()
            .filter(|&elem| elem.id() == target_id)
            .collect();

        if matching_positions.len() == 0 {
            return Err(
                format!("Did not find a grid postion with target_id {}", target_id).to_owned(),
            );
        }
        if matching_positions.len() >= 1 {
            return Err("Found multiple grid positions with same id!".to_owned());
        }
        return Ok(matching_positions.first().unwrap());
    }

    pub fn shuffle_spins(&mut self){
        //function to randomly set all spins in a grid without changing the links
        // between the spins
        for grid_pos in self. grid_positions.iter_mut(){
            match grid_pos.particle.as_mut(){
                Some(particle) => {
                    particle.spin = rand::thread_rng().gen()
                }
                None => (),
            }
        }
    }

    pub fn grid_positions(&self) -> &Vec<GridPos> {
        return &self.grid_positions;
    }
    pub fn id(&self) -> u32 {
        return self.id;
    }
    pub fn T(&self) -> f64 {
        return self.T;
    }
    pub fn set_T(&mut self, T: f64) {
        self.T = T
    }
    pub fn external_field(&self) -> f64 {
        return self.external_field;
    }
    pub fn set_external_field(&mut self, external_field: f64) {
        self.external_field = external_field
    }
    pub fn set_id(&mut self, id: u32) {
        self.id = id
    }
    pub fn update_history(&mut self, linked_overlapp:f64) {
        let run = self.run;
        let T = self.T;
        let energy = self.calc_energy();
        let magnetization = self.calc_magnetization();
        let av_energy = self.calc_av_energy();
        let katz_energy = 
            1.0-self.T*energy.abs()/(self.dimensions.len() as f64 * self.config.coupling_variance*self.capacity as f64);
        let av_linked_overlapp = self.calc_av_linked_overlapp();
        let av_katz_energy = self.calc_av_katz();

        self.history.add(
            run, T, energy, magnetization, linked_overlapp, katz_energy, av_linked_overlapp, av_katz_energy);
        
        if run%self.config.skip_save==0{
            self.write_to_output_file(run, T, energy, magnetization,linked_overlapp);
        }
    }
    pub fn history(&self) -> &History {
        return &self.history;
    }
    pub fn dimensions(&self) -> &Vec<u16> {
        return &self.dimensions;
    }
    pub fn capacity(&self) -> u32 {
        return self.capacity;
    }
    pub fn set_sweep(&mut self, sweep: Sweep) {
        self.T = sweep.T_start;
        self.sweep = Some(sweep)
    }
    pub fn sweep(&self) -> Option<&Sweep> {
        return self.sweep.as_ref();
    }
}

use std::{cmp::min, collections::VecDeque, fmt, fs::File, io::Write};
impl fmt::Display for Grid {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} Grid {} \n", "-".repeat(20), "-".repeat(20));
        let mut output = String::new();
        output.push_str("Dimension: ");

        for axis in self.dimensions.iter() {
            output.push_str(&axis.to_string());
            output.push('x');
        }
        output.pop();
        write!(f, "{} \n", output);
        write!(f, "Capacity: {} particles \n", self.capacity.to_string());
        write!(f, "Temperature: {}K \n", self.T.to_string());

        for (grid_idx, grid_pos) in self.grid_positions.iter().enumerate() {
            match grid_pos.particle() {
                None => write!(f, "No particle at grid index {}\n", grid_idx),
                Some(particle) => write!(f, "{}\n", particle.to_string()),
            };
        }
        write!(f, "")
    }
}

impl Clone for Grid {
    fn clone(&self) -> Self {
        let output_file: Option<fs::File> = None;
        let grid = Grid {
            id: self.id,
            run: self.run,
            cloned_from: Some(self.id),
            config: self.config.clone(),
            dimensions: self.dimensions.clone(),
            capacity: self.capacity,
            num_particles: self.num_particles,
            T: self.T,
            external_field: self.external_field,
            grid_positions: self.grid_positions.clone(),
            history: self.history.clone(),
            sweep: self.sweep.clone(),
            output_file,
        };
        return grid;
    }
}

/* #[cfg(test)]
mod tests{
    use crate::{*,particles::*};

    #[test]
    fn calc_neighbour_indices(){
        //create grid, calc neighbour idx check for correctness
        //test config produces 5x5 grid
        let config: config::Config = config::load();
        let my_grid = grid::Grid::new(config.grid_config.clone(),0).unwrap();
        let neighbour_idx =
            my_grid.calc_filled_neighbour_indices(0);

        assert_eq!(neighbour_idx[0],4);
        assert_eq!(neighbour_idx[1],1);
        assert_eq!(neighbour_idx[2],20);
        assert_eq!(neighbour_idx[3],5);
    }

    #[test]
    fn neighbour_positions(){
        let config: config::Config = config::load();
        let my_grid = grid::Grid::new(config.grid_config.clone(),0).unwrap();
        let neighbour_positions =
            my_grid.calc_neighbour_positions(0);

        assert_eq!(neighbour_positions[0], vec![4,0]);
        assert_eq!(neighbour_positions[1], vec![1,0]);
        assert_eq!(neighbour_positions[2], vec![0,4]);
        assert_eq!(neighbour_positions[3], vec![0,1]);
        }

    #[test]
    fn idx_to_pos(){
        let config: config::Config = config::load();
        let my_grid = grid::Grid::new(config.grid_config.clone(),0).unwrap();

        let idx = 7 as u32;
        let pos = my_grid.idx_to_pos(idx).unwrap();

        assert_eq!(pos, vec![2,1]);

        let to_large_idx = 100 as u32;
        let _err = my_grid.idx_to_pos(to_large_idx).unwrap_err();

        }

    #[test]
    fn pos_to_idx(){
        let config: config::Config = config::load();
        let my_grid = grid::Grid::new(config.grid_config.clone(),0).unwrap();

        let valid_pos = vec![2,1];
        let idx = my_grid.pos_to_id(&valid_pos).unwrap();

        assert_eq!(idx,7);

        let invalid_pos = vec![10,1,3];
        let _err = my_grid.pos_to_id(&invalid_pos).unwrap_err();
    }

    #[test]
    fn next_pos(){
        let config: config::Config = config::load();
        let my_grid = grid::Grid::new(config.grid_config.clone(),0).unwrap();

        let mut pos = vec![3,0];
        pos = my_grid.next_pos(&pos).unwrap();
        assert_eq!(pos, vec![4,0]);
        pos = my_grid.next_pos(&pos).unwrap();
        assert_eq!(pos, vec![0,1]);
    }
} */
