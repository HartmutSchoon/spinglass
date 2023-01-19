#![warn(clippy::all, rust_2018_idioms)]
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]
// hide console window on Windows in release

use core::num;

// When compiling natively:
#[cfg(not(target_arch = "wasm32"))]
fn main() {
    //run_without_ui();
    run_with_ui();
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
    use spinglass::Simulation;
    use std::time::{Duration, Instant};

    let num_sweeps = 100_000;
    let mut sim = Simulation::new();
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

/*TODO: 
    -Die couplings müssen nach Gauß initiert werden
    -Multithreading in WASM
    -Ab und zu stürzt der ab bei PT wenn man ein grid schließt

*/