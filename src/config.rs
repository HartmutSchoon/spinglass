use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::prelude::*;

use crate::grid::Sweep;

#[derive(Debug,Deserialize,Serialize,Clone,PartialEq,Copy)]
pub enum ParticlePositionMethod{
    RandomParticlePositions,
    FillGrid,
}

#[derive(Debug,Deserialize,Serialize,Clone)]
pub struct Config{
    pub grid_config: GridConfig,
    pub pt_config: PTConfig,
    pub simulation_config:SimulationConfig
}

impl Config{

    pub fn load()-> Config{
        let mut config_str = String::new();
        File::open("./config.toml").unwrap().read_to_string(&mut config_str);
        let config:Config = toml::from_str(&config_str).unwrap();
    return config
}

}

impl Default for Config{
    fn default()->Self{
        return Config{
            grid_config: GridConfig::default(),
            pt_config: PTConfig::default(),
            simulation_config:SimulationConfig::default(),
        };
    }
}




#[derive(Debug,Deserialize,Serialize,Clone)]
pub struct GridConfig{
    pub grid_dimensions:Vec<u16>,
    pub particle_position_method: ParticlePositionMethod,
    pub num_particles:u32,
    pub coupling_mean: f64,
    pub coupling_variance: f64,
    pub T:f64,
    pub spin_up_init: bool,
    pub external_field:f64,
    pub history_capacity:usize,  
    
    //path to results folder where grid data is stored
    pub save_path: String,
    //Amount of skipped datapoints while saving. E.x. skip_save=100
    //results in only saving every 100 Datapoints
    pub skip_save: u32,
}

impl Default for GridConfig{
    fn default() -> Self{
        return GridConfig{
            grid_dimensions:vec![100,100],
            particle_position_method: ParticlePositionMethod::FillGrid,
            num_particles:1000,
            T:0.5,
            spin_up_init: false,
            external_field:0.0,
            coupling_mean: 1.0,
            coupling_variance: 0.0,
            history_capacity:10_000,
            save_path: String::from("./results"),
            skip_save: 1_000_000,
        }
    }
}

#[derive(Debug,Deserialize,Serialize,Clone)]
pub struct SimulationConfig{
    pub num_grids:u32,
    pub steps_per_sweep:u32,
    pub num_sweeps:u32, 
    pub histo_width: f64,

}

impl Default for SimulationConfig{
    fn default() -> Self{
        return SimulationConfig{
            num_grids:0,
            steps_per_sweep:1000,
            num_sweeps:10,
            histo_width: 0.1,
            
        }
    }
}

#[derive(Debug,Deserialize,Serialize,Clone)]
pub struct PTConfig{
    pub num_grids_equal_T: u32,
    pub T_start: f64,
    pub T_end: f64,
    //deltaT(T) is is produced with linear spacing-> deltaT grows linear
    //delta(T) = linear_m*T+linear_dT0
    pub linear_m: f64,
    pub linear_dT0:f64,
}

impl Default for PTConfig{
    fn default() -> Self{
        return PTConfig{
            num_grids_equal_T: 2,
            T_start: 0.5,
            T_end: 0.8,
            linear_m: 0.035,
            linear_dT0: 0.03,
        }
    }
}
