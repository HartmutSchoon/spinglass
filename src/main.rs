#![warn(clippy::all, rust_2018_idioms)]
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]
// hide console window on Windows in release

use core::num;
use std::{env,fs};

// When compiling natively:
#[cfg(not(target_arch = "wasm32"))]
fn main(){
    //Read user input

    match env::args().find(|elem|elem == "no_ui"){
        Some(_) => run_without_ui(),
        None => run_with_ui(),
    }

}

fn run_with_ui(){
    
    // Log to stdout (if you run with `RUST_LOG=debug`).
    tracing_subscriber::fmt::init();

    let native_options = eframe::NativeOptions::default();
    eframe::run_native(
        "eframe template",
        native_options,
        Box::new(|cc| Box::new(spinglass::App::new(cc))),
    );
}

fn run_without_ui(){
    use spinglass::{Simulation,Config};    
    use std::time::{Duration, Instant};


    let mut config = Config::load(); 
    if let Some(p) = env::args().find(|elem|elem.contains("path=")) {
        let path = p.replace("path=", "");
        config.grid_config.save_path = path.clone();
        fs::remove_dir_all(path.clone());
        fs::create_dir(path.clone());
    };
    let num_sweeps = config.simulation_config.num_sweeps;

    let mut sim = Simulation::new(config.clone());
    sim.custom_run();

    let start = Instant::now();
    let mut duration:Duration;
    
    for run in 0..num_sweeps{
        sim.simulation_step();
        duration = start.elapsed();
        println!("Completed run {}/{}, ETC: {:.2}s",
            run+1,
            num_sweeps,
            (num_sweeps-run) as f64 * duration.as_secs_f64()/run as f64 );
    }
}

// when compiling to web using trunk.
#[cfg(target_arch = "wasm32")]
fn main() {
    // Make sure panics are logged using `console.error`.
    console_error_panic_hook::set_once();

    // Redirect tracing to console.log and friends:
    tracing_wasm::set_as_global_default();

    let web_options = eframe::WebOptions::default();

    wasm_bindgen_futures::spawn_local(async {
        eframe::start_web(
            "the_canvas_id", // hardcode it
            web_options,
            Box::new(|cc| Box::new(spinglass::App::new(cc))),
        )
        .await
        .expect("failed to start eframe");
    });
}

