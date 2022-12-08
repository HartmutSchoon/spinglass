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


#[derive(Deserialize,Clone)]
pub struct PTConfig{
    pub num_T_steps: u32,
    pub num_grids_equal_T: u32,
    pub T_start: f64,
    pub T_end: f64,
}
#[derive(Deserialize,Clone)]
pub struct GridConfig{
    pub grid_dimensions:Vec<u16>,
    pub particle_position_method: ParticlePositionMethod,
    pub num_particles:u32,
    pub coupling_limits: [f64;2],
    pub T:f64,
    pub spin_up_init: bool,
    pub external_field:f64,
    pub history_capacity:usize,
}

#[derive(Deserialize,Clone)]
pub struct SimulationConfig{
    pub num_grids:u32,
    pub thread_steps:u32,
    pub num_steps:u32,
}

#[derive(Deserialize,Clone)]
pub struct AppConfig{
}

pub fn load()-> Config{
    let mut config_str = String::new();
    File::open("config.toml").unwrap().read_to_string(&mut config_str);
    let config:Config = toml::from_str(&config_str).unwrap();
    return config
}