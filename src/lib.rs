#![warn(clippy::all, rust_2018_idioms)]

mod config;
pub use config::*;

mod grid;
pub use grid::*;

mod particles;
pub use particles::*;

mod simulation;
pub use simulation::*;

mod app;
pub use app::App;
