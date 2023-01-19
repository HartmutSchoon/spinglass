use serde::Deserialize;
use std::fs::File;
use std::io::prelude::*;

use crate::grid::Sweep;

#[derive(Deserialize,Clone,PartialEq,Copy)]
pub enum ParticlePositionMethod{
    RandomParticlePositions,
    FillGrid,
}

#[derive(Deserialize,Clone)]
pub struct Config{
    pub grid_config: GridConfig,
    pub pt_config: PTConfig,
    pub simulation_config:SimulationConfig,
    pub app_config:AppConfig,
}

impl Default for Config{
    fn default()->Self{
        return Config{
            grid_config: GridConfig::default(),
            pt_config: PTConfig::default(),
            simulation_config:SimulationConfig::default(),
            app_config:AppConfig::default(),
        };
    }
}


#[derive(Deserialize,Clone)]
pub struct PTConfig{
    pub num_T_steps: u32,
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
            num_T_steps: 10,
            num_grids_equal_T: 4,
            T_start: 0.5,
            T_end: 0.8,
            linear_m: 0.035,
            linear_dT0: 0.03,
        }
    }
}

#[derive(Deserialize,Clone)]
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
}

impl Default for GridConfig{
    fn default() -> Self{
        return GridConfig{
            grid_dimensions:vec![100,100],
            particle_position_method: ParticlePositionMethod::FillGrid,
            num_particles:1000,
            coupling_mean: 1.0,
            coupling_variance: 0.0,
            T:0.5,
            spin_up_init: false,
            external_field:0.0,
            history_capacity:10_000,
        }
    }
}

#[derive(Deserialize,Clone)]
pub struct SimulationConfig{
    pub num_grids:u32,
    pub thread_steps:u32,
    pub num_steps:u32,
    pub histo_width: f64,
}

impl Default for SimulationConfig{
    fn default() -> Self{
        return SimulationConfig{
            num_grids:0,
            thread_steps:1000,
            num_steps:100000,
            histo_width: 0.1,
        }
    }
}

#[derive(Deserialize,Clone)]
pub struct AppConfig{
}

impl Default for AppConfig{
    fn default()->Self{
        return AppConfig{};
    }
}

pub fn load()-> Config{
    let config = Config::default();
    return config
}