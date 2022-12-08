
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
}
impl AppState{
    pub fn new(sweep: grid::Sweep)->Self{
        return Self {sweep};
    }

}
impl Default for AppState {
    fn default() -> Self {
        Self { sweep: Default::default() }
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
        let config = config::load(); 
        let mut sim = Simulation::new();
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

                let mut min_coupling:f32 = self.sim.default_grid_config.coupling_limits[0] as f32;
                ui.add(DragValue::new(&mut min_coupling).speed(0.01).prefix("Min coupling: "));
                self.sim.default_grid_config.coupling_limits[0] = min_coupling as f64;

                let mut max_coupling:f32 = self.sim.default_grid_config.coupling_limits[1] as f32;
                ui.add(DragValue::new(&mut max_coupling).
                        speed(0.01).prefix("Max coupling:: ").clamp_range(min_coupling..=f32::INFINITY));
                self.sim.default_grid_config.coupling_limits[1] = max_coupling as f64;



                if ui.add(egui::Button::new("Create!")).clicked(){
                    self.sim.new_grid();
                };
            });
            Frame::group(ui.style())
            .show(ui,|ui|{
                ui.label("Simulation Parameters");
                ui.label("#Steps Grids run in Multithreading: ");
                let mut thread_steps = self.sim.config.thread_steps;
                ui.add(DragValue::new(&mut thread_steps).speed(100)
                    .prefix("#Threaded Steps: "));
                self.sim.config.thread_steps = thread_steps;
                if ui.add(egui::Button::new("Run Simulation")).clicked(){
                    self.sim.running=true;
                };
                if ui.add(egui::Button::new("Pause Simulation")).clicked(){
                    self.sim.running=false;
                };
                if ui.add(egui::Button::new("Delete all grids")).clicked(){
                    self.sim.queue_all_grids_deletion();
                };

            });

            ui.label("Nur für sichtbarmachung von PT");
            if ui.add(egui::Button::new("Initialize PT:")).clicked(){
                self.sim.pt_init();
            };
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

                    });
            
                    Plot::new(format!("Plot {}",grid_id))
                    .include_x(0.0)
                    .include_x(100.0)
                    .legend(Legend::default())
                    .show(ui, |plot_ui|{
                        let history = self.sim.history_by_id(grid_id).unwrap();
                        let T:PlotPoints = (0..history.current_size()).map(|i|{
                            let x = (i as f64 * self.sim.config.thread_steps as f64);
                            let T= history.T();
                            [x,T[i]]
                            }).collect();
                        let t_line=Line::new(T);
                        let E:PlotPoints = (0..history.current_size()).map(|i|{
                            let x = (i as f64 * self.sim.config.thread_steps as f64);
                            let energy= history.energy();
                            [x,energy[i]]
                        }).collect();
                        let e_line = Line::new(E);
                        plot_ui.line(t_line);
                        plot_ui.line(e_line);
                    });
                
                });//Window
            }//Grid Loop
            if self.sim.pt_enviroment.grids.len() > 0{
                let max_grids_same_T = 
                    self.sim.pt_enviroment.pt_ids.iter().map(|e|e.ids.len()).max().unwrap() as f32;
                let num_diff_T = self.sim.pt_enviroment.pt_ids.len() as f32;
                Window::new("Parallel Tempering").show(ctx, |win_ui| {
                    let available_height = win_ui.available_height()-10.0;
                    let available_width = win_ui.available_width()-10.0;
                    //ui.vertical(|ui|{
                        for  equal_T_grids in self.sim.pt_enviroment.pt_ids.iter(){
                            win_ui.with_layout(Layout::left_to_right(Align::TOP), |row_ui|{                               
                                //ui.label(format!("T: {}",equal_T_grids.T.to_string()));
                                for &grid_id in equal_T_grids.ids.iter(){
                                    Frame::group(row_ui.style())
                                    //.outer_margin(Margin{left: 10.0, right: 10.0, top: 10.0, bottom: 10.0})
                                    .show(row_ui, |group_ui|{
                                        let group_size = Vec2{
                                            x: available_width/max_grids_same_T-0.0,
                                            y: available_height/num_diff_T-30.0,
                                        };
                                        group_ui.allocate_exact_size(group_size,Sense::hover());
                                        
                                        Frame::canvas(group_ui.style())
                                        .fill(Color32::GRAY)
                                        //.outer_margin(Margin{left:10.0, right: 10.0, top: 10.0, bottom: 10.0})
                                        .show(group_ui, |frame_ui|{
        
                                            //ui.ctx().request_repaint();
                                            /* let desired_canvas_size = Vec2 { 
                                                x: available_width/max_grids_same_T-10.0,
                                                y: available_height/num_diff_T-10.0}; */
                                            let desired_canvas_size = group_size*Vec2 { x: 0.5, y: 0.5 };
                                            let (response, painter) = 
                                                frame_ui.allocate_painter(desired_canvas_size, Sense::hover());
                                            let canvas_size = response.rect; 

                                            /* let desired_canvas_size = size*Vec2 { x: 0.5, y: 0.5 };
                                            let (canvas_size, _) = 
                                                ui.allocate_exact_size(desired_canvas_size, Sense::hover());

                                              */
                                            let grid = self.sim.pt_enviroment.grid(grid_id).unwrap();
        
                                            let mut mesh = self.calc_rectangles_from_grid(
                                                canvas_size,
                                                grid);
                                            painter.add(mesh);
        
                                        });

                                    });
                                    
                                }
                            });
                        }

                    //});
                });     
            }
            
            self.sim.simulation_step();
        });//CentralPanel
        self.sim.delete_queded_grids();
    }
}




    
