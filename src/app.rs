
use eframe::egui::{self,*,style::*,plot::*};
use eframe::epaint::RectShape;


use crate::config::{self, Config};
use crate::{particles, grid};
use crate::simulation::Simulation;

pub struct ParticleRect{
    pos: [f64;4],
    color: [f32;4],
}

struct AppState{
    pub sweep: grid::Sweep,
    pub show_overlapp: bool,
    pub show_pt_acceptance: bool,
}
impl AppState{
    pub fn new(sweep: grid::Sweep)->Self{
        return Self {
            sweep,
            show_overlapp:false,
            show_pt_acceptance: false};
    }

}
impl Default for AppState {
    fn default() -> Self {
        return Self { 
            sweep: Default::default(),
            show_overlapp: false,
            show_pt_acceptance:false,};
    }
}

//TODO: was passiert hier eigl bei 3D grids? wird nur vorne geplottet?
pub struct App {
    pub sim: Simulation,
    pub config: Config, //TODO: der ist glaube ich unnötig, benutze gerade nur simConfig!
    app_state: AppState,
}

impl App{
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        //load config and split into smaller parts
        let config = config::Config::default(); 
        let mut sim = Simulation::new(config.clone());
        let app_state = AppState::default();

        return App {
            sim,
            config,
            app_state,
        }   
    }

    fn calc_rectangles_from_grid(& self, canvas_size:Rect, grid: &grid::Grid)
    ->Mesh{

        let mut mesh = Mesh::default();

        let grid_x_dim = f32::from(grid.config.grid_dimensions[0]);
        let grid_y_dim = f32::from(grid.config.grid_dimensions[1]);

        let rectangle_width:f32 = (canvas_size.width() / grid_x_dim).ceil();
        let rectangle_height:f32 = (canvas_size.height() / grid_y_dim).ceil();

        let mut rect_buff: Vec<Shape> = Vec::new();
        //Das ist falsch so! der soll nicht alle grid_positions zeichnen sondern nur die in den ersten zwei dims
        for grid_pos in grid.grid_positions(){
            let color:Color32;
            match grid_pos.particle(){
                None => color=Color32::GRAY,
                Some(particle) =>{
                    color = match particle.spin(){
                        particles::Spin::Up => Color32::WHITE,
                        particles::Spin::Down => Color32::BLACK,
                    };
                }
            }

            let pos = grid_pos.pos();
            let min = (Vec2{
                    x: f32::from(pos[0])*rectangle_width,
                    y: f32::from(pos[1])*rectangle_height}
                +canvas_size.min.to_vec2())
                .to_pos2();
            

            let max = (Vec2{
                x: f32::from(pos[0]+1)*rectangle_width,
                y: f32::from(pos[1]+1)*rectangle_height}
                +canvas_size.max.to_vec2())
                .to_pos2();

           /*  rect_buff.push(Shape::Rect(RectShape{
                rect: Rect{min,max},
                rounding: Rounding::none(),
                fill: color,
                stroke: Stroke::default()})) */
            mesh.add_colored_rect(
                Rect {min,max},
                color
            );
        }
        return mesh;           
    }
}

impl eframe::App for App{
    fn update(&mut self, ctx: &egui::Context, frame: &mut eframe::Frame) {
        SidePanel::right("Menue").show(ctx, |ui|{
            Frame::group(ui.style())
            .show(ui,|ui|{
                ui.label("Create new Grid");

                let mut grid_dimensions = self.sim.default_grid_config.grid_dimensions.clone();
                let mut num_dims = grid_dimensions.len();
                ui.add(DragValue::new(&mut num_dims).speed(0.1)
                        .prefix("#Dimensions: ").clamp_range(1..=10));
                grid_dimensions.resize(num_dims, 0);

                for axis in 0..num_dims{
                    let mut num_particles = grid_dimensions[axis]; 
                    ui.add(DragValue::new(&mut num_particles).speed(1)
                        .prefix(format!("[{}]: ",axis)).clamp_range(1..=u16::MAX));
                        grid_dimensions[axis] = num_particles;
                }
                self.sim.default_grid_config.grid_dimensions = grid_dimensions.clone();
                let capacity: f64 = grid_dimensions.iter().map(|&e|e as f64).product();
                let scaled_capacity = capacity/1000.0;
                ui.label(format!("Grid capacity: {}K",scaled_capacity));

                let mut position_method = self.sim.default_grid_config.particle_position_method.clone();
                ui.radio_value(&mut position_method,
                  config::ParticlePositionMethod::FillGrid,
                         "Fill grid");
                ui.radio_value(&mut position_method,
                config::ParticlePositionMethod::RandomParticlePositions,
                        "Random Positions");
                self.sim.default_grid_config.particle_position_method = position_method;
                if position_method == config::ParticlePositionMethod::RandomParticlePositions{
                    let mut num_particles = self.sim.default_grid_config.num_particles;
                    let capacity = grid_dimensions.iter().map(|&e|u32::from(e)).product();
                    ui.add(DragValue::new(&mut num_particles).speed(1)
                        .prefix("#Particles: ").clamp_range(0..=capacity));
                    self.sim.default_grid_config.num_particles = num_particles;
                }
                let mut spin_up_init = self.sim.default_grid_config.spin_up_init;
                ui.add(Checkbox::new(&mut spin_up_init, "All Spins Up "));
                self.sim.default_grid_config.spin_up_init = spin_up_init;

                let mut T:f32 = self.sim.default_grid_config.T as f32;
                ui.add(DragValue::new(&mut T).speed(0.01).prefix("Temperature: "));
                self.sim.default_grid_config.T = T as f64;

                let mut external_field:f32 = self.sim.default_grid_config.external_field as f32;
                ui.add(DragValue::new(&mut external_field).speed(0.01).prefix("External field: "));
                self.sim.default_grid_config.external_field = external_field as f64;

                let mut coupling_mean:f32 = self.sim.default_grid_config.coupling_mean as f32;
                ui.add(DragValue::new(&mut coupling_mean).speed(0.01).prefix("Coupling mean: "));
                self.sim.default_grid_config.coupling_mean = coupling_mean as f64;
                
                let mut coupling_variance:f32 = self.sim.default_grid_config.coupling_variance as f32;
                ui.add(DragValue::new(&mut coupling_variance).speed(0.01).prefix("Coupling variance: "));
                self.sim.default_grid_config.coupling_variance = coupling_variance as f64;

                if ui.add(egui::Button::new("Create!")).clicked(){
                    self.sim.new_grid();
                };
            });
            Frame::group(ui.style())
            .show(ui,|ui|{
                ui.label("Simulation Parameters");
                ui.label("#Steps Grids run in Multithreading: ");
                let mut steps_per_sweep = self.sim.config.steps_per_sweep;
                let mut slider_steps = steps_per_sweep as f64 / 1000.0;
                ui.add(DragValue::new(&mut slider_steps).speed(100)
                    .prefix("#Threaded Steps: ").suffix("K"));
                steps_per_sweep = (slider_steps * 1000.0) as u32;
                self.sim.config.steps_per_sweep = steps_per_sweep;
                if ui.add(egui::Button::new("Run Simulation")).clicked(){
                    self.sim.running=true;
                };
                if ui.add(egui::Button::new("Pause Simulation")).clicked(){
                    self.sim.running=false;
                };
                if ui.add(egui::Button::new("Delete all grids")).clicked(){
                    self.sim.queue_all_grids_deletion();
                };

                if ui.add(egui::Button::new("Start Custom Run")).clicked(){
                    self.sim.custom_run();
                }
            });
            Frame::group(ui.style())
            .show(ui,|ui|{
                ui.label("Additional Windows");
                if ui.add(egui::Button::new("Show Overlapp")).clicked(){
                    self.app_state.show_overlapp = !self.app_state.show_overlapp;
                }
                if ui.add(egui::Button::new("Show PT Acceptance Ration")).clicked(){
                    self.app_state.show_pt_acceptance = !self.app_state.show_pt_acceptance;
                }

            });
        });

        CentralPanel::default()
        .show(ctx, |ui|{
            let current_grid_ids = self.sim.current_grid_ids().clone();
            for grid_id in current_grid_ids{
                Window::new(format!("Grid {}",grid_id)).show(ctx, |ui| {
                    ui.with_layout(Layout::left_to_right(Align::TOP), |ui|{
                        Frame::canvas(ui.style())
                        .fill(Color32::GRAY)
                        //.outer_margin(Margin{left:10.0, right: 10.0, top: 10.0, bottom: 10.0})
                        .show(ui, |ui|{
                            ui.ctx().request_repaint();
                            let desired_canvas_size = ui.available_size()*Vec2{x:0.5, y:0.5};
                            let (response, painter) = 
                                ui.allocate_painter(desired_canvas_size, Sense::hover());
                            let canvas_size = response.rect;
                            
                            let grid = self.sim.grid(grid_id).unwrap();
        
                            let mut mesh = self.calc_rectangles_from_grid(
                                canvas_size,
                                grid);
                            painter.add(mesh);
                        });//Frame
        
                        Grid::new(format!("settings grid {}",grid_id)).show(ui, |ui|{
                            egui::Frame::group(ui.style())
                            .show(ui, |ui| {
                                ui.vertical(|ui|{
                                    
                                    let mut output_string = String::new();
                                    let mut dimensions = self.sim.grid(grid_id).unwrap().dimensions();
                                    output_string.push('[');
                                    for axis in dimensions{
                                        output_string.push_str(&axis.to_string());
                                        output_string.push(',');
                                    }
                                    output_string.pop();
                                    output_string.push(']');
                                    ui.label(format!("Dimensions: {}",output_string));
                                    let capacity = self.sim.grid(grid_id).unwrap().capacity();
                                    let scaled_capacity = (capacity as f64)/1000.0;
                                    ui.label(format!("Grid capacity: {}K",scaled_capacity));

                                
                                    if let Some(id) = self.sim.grid(grid_id).unwrap().cloned_from{
                                        ui.label(format!("Clone of Grid {}",id));
                                    }

                                    let coupling_mean = self.sim.grid(grid_id).unwrap().config.coupling_mean;
                                    let coupling_variance = self.sim.grid(grid_id).unwrap().config.coupling_variance;
                                    ui.label(format!("Coupling mean: {:.2}", coupling_mean));
                                    ui.label(format!("Coupling variance: {:.2}", coupling_variance));
            
                                    let mut T = self.sim.grid(grid_id).unwrap().T();
                                    ui.add(DragValue::new(&mut T).speed(0.01).prefix("Temperature: "));
                                    self.sim.grid_mut(grid_id).unwrap().set_T(T);
            
                                    let mut external_field = self.sim.grid(grid_id).unwrap().external_field();
                                    ui.add(DragValue::new(&mut external_field).speed(0.01).prefix("External Field: "));
                                    self.sim.grid_mut(grid_id).unwrap().set_external_field(external_field); 
            
                                    
                                    ui.horizontal(|ui| {
                                        if ui.add(egui::Button::new("Clone Grid")).clicked(){
                                            self.sim.clone_grid(grid_id);
                                        };
            
                                        if ui.add(egui::Button::new("Delete Grid")).clicked(){
                                            self.sim.queue_grid_deletion(grid_id);
                                        };
                                    });
                                });
                            });
            
                            egui::Frame::group(ui.style())
                            .show(ui, |ui| {
                                ui.vertical(|ui|{
                                    ui.label("Linear Temperature Sweep");
                                    let mut T_start = self.app_state.sweep.T_start();
                                    ui.add(DragValue::new(&mut T_start).speed(0.01).prefix("Start: "));
                                    self.app_state.sweep.set_T_start(T_start);
            
                                    let mut T_end = self.app_state.sweep.T_end();
                                    ui.add(DragValue::new(&mut T_end).speed(0.01).prefix("End: "));
                                    self.app_state.sweep.set_T_end(T_end);
            
                                    let mut num_steps = self.app_state.sweep.num_steps();
                                    let mut slider_steps:f64 = num_steps as f64 / 1000000.0;
                                    ui.add(DragValue::new(&mut slider_steps).speed(0.1)
                                        .prefix("#Steps: ").suffix("M"));
                                    num_steps = (slider_steps * 1000000.0)as u64;
                                    self.app_state.sweep.set_num_steps(num_steps);
                                    
            
                                    if ui.add(egui::Button::new("Run!")).clicked(){
                                        self.sim.grid_mut(grid_id).unwrap().set_sweep(
                                            self.app_state.sweep.clone()
                                        );
                                    };
            
            
                                });
                            });

                            ui.end_row();
                            egui::Frame::group(ui.style())
                            .show(ui, |ui| {
                                ui.vertical(|ui|{
                                    ui.label("Parallel Tempering");

                                    ui.add(DragValue::new(&mut self.sim.pt_enviroment.config.num_grids_equal_T)
                                        .speed(0.1)
                                        .prefix("#Grids per temperature: "));

                                    ui.add(DragValue::new(&mut self.sim.pt_enviroment.config.T_start)
                                        .speed(0.1)
                                        .prefix("Lowest temperature: "));
                                    
                                    ui.add(DragValue::new(&mut self.sim.pt_enviroment.config.T_end)
                                        .speed(0.1)
                                        .prefix("highest temperature: "));
                                    

                                    if ui.add(Button::new("Initialize!")).clicked(){
                                        self.sim.pt_init(grid_id);
                                    }

                                });
                            });
                        });
        
                    });
                    
            
                    Plot::new(format!("Plot {}",grid_id))
                    .include_x(0.0)
                    //.include_x(100.0)
                    .legend(Legend::default())
                    .show(ui, |plot_ui|{
                        let history = self.sim.history_by_id(grid_id).unwrap();
                        let T:PlotPoints = (0..history.current_size()).map(|i|{
                            let x = i as f64;
                            let T= history.T[i];
                            [x,T]
                            }).collect();
                        let t_line=Line::new(T).name("Temperature");
                        plot_ui.line(t_line);
                        let E:PlotPoints = (0..history.current_size()).map(|i|{
                            let x = i as f64;
                            let energy= history.energy[i]/self.sim.grid(grid_id).unwrap().capacity as f64;
                            [x,energy]
                        }).collect();
                        let e_line = Line::new(E).name("Energy");
                        plot_ui.line(e_line);

                        let magnetization:PlotPoints = (0..history.current_size()).map(|i|{
                            let x = i as f64;
                            let mag= history.magnetization[i] as f64/self.sim.grid(grid_id).unwrap().capacity as f64;
                            [x,mag]
                        }).collect();
                        let mag_line = Line::new(magnetization).name("Magnetization");
                        plot_ui.line(mag_line);

                        let linked_overlapp:PlotPoints = (0..history.current_size()).map(|i|{
                            let x = i as f64;
                            let overlapp = history.linked_overlapp[i] as f64;
                            [x,overlapp]
                        }).collect();
                        let overlapp_line = Line::new(linked_overlapp).name("Linked Overlapp");
                        plot_ui.line(overlapp_line);

                        let katz_energy:PlotPoints = (0..history.current_size()).map(|i|{
                            let x = i as f64;
                            let katz_en = history.katz_energy[i] as f64;
                            [x,katz_en]
                        }).collect();
                        let ke_line = Line::new(katz_energy).name("1-2T|U|/zJ^2");
                        plot_ui.line(ke_line);

                        let av_linked_overlapp:PlotPoints = (0..history.current_size()).map(|i|{
                            let x = i as f64;
                            let overlapp = history.av_linked_overlapp[i] as f64;
                            [x,overlapp]
                        }).collect();
                        let overlapp_line = Line::new(av_linked_overlapp).name("<ql>");
                        plot_ui.line(overlapp_line);

                        let av_katz_energy:PlotPoints = (0..history.current_size()).map(|i|{
                            let x = i as f64;
                            let av_katz_en = history.av_katz_energy[i] as f64;
                            [x,av_katz_en]
                        }).collect();
                        let av_ke_line = Line::new(av_katz_energy).name("<1-2T|U|/zJ^2>");
                        plot_ui.line(av_ke_line); 
                        
                    });
                
                });
            }//Grid Loop
            if self.app_state.show_overlapp {
                Window::new("Overlapp").show(ctx, |ui| {
                    Plot::new("Overlapp_Plot")
                    .legend(Legend::default())
                    .show(ui, |plot_ui|{
                        let overlapp:PlotPoints = (0..self.sim.overlapp_histo.len()).map(|i|{
                            let x = i as f64*self.sim.config.histo_width-1.0;
                            let overlapp= (self.sim.overlapp_histo[i] as f64/self.sim.overlapp_histo.len() as f64);
                            [x,overlapp]
                        }).collect(); 
                        let overlapps_line = Line::new(overlapp).name("Overlapp");
                        plot_ui.line(overlapps_line);
                        let linked_overlapp:PlotPoints = (0..self.sim.linked_overlapp_histo.len()).map(|i|{
                            let x = i as f64*self.sim.config.histo_width-1.0;
                            let overlapp= (self.sim.linked_overlapp_histo[i] as f64/self.sim.overlapp_histo.len() as f64);
                            [x,overlapp]
                        }).collect(); 
                        let linked_overlapps_line = Line::new(linked_overlapp).name("Linked Overlapp");
                        plot_ui.line(linked_overlapps_line);
                    });

                });
            };
            if self.app_state.show_pt_acceptance {
                Window::new("Parallel Tempering").show(ctx, |ui| {
                    for equalTGridIDs in self.sim.pt_enviroment.pt_ids.iter(){
                        ui.label(format!(
                            "T: {:.2}, Acceptance Prob: {:.2} %",
                            equalTGridIDs.T,
                            equalTGridIDs.current_pt_acceptance_prob*100.0));
                    }
                });
            }
            self.sim.simulation_step();
        });//CentralPanel
        self.sim.delete_queded_grids();
/*         if self.sim.current_grid_ids.len() == 0{
            
        } */
    }

}




    
