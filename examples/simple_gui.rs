//! Simple GUI demonstration of the tectonic plate simulation
//! 
//! This example shows:
//! - Real-time plate movement visualization
//! - Plate statistics (type, area, age)
//! - Interactive simulation controls
//! 
//! To run this example:
//! ```bash
//! cargo run --example simple_gui --features gui
//! ```

use eframe::egui;
use tectonics::{TectonicWorld, SphericalPoint, PlateType, BoundaryType, WorldConfig};
use std::f64::consts::PI;

fn main() -> Result<(), eframe::Error> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1000.0, 700.0]),
        ..Default::default()
    };
    
    eframe::run_native(
        "Tectonic Plate Simulation",
        options,
        Box::new(|_cc| Ok(Box::new(TectonicApp::new()))),
    )
}

#[derive(Clone, Copy, PartialEq)]
enum ProjectionType {
    Equirectangular,
    Orthographic,
}

struct TectonicApp {
    world: TectonicWorld,
    time: f64,
    dt: f64,
    is_running: bool,
    step_count: u32,
    plate_colors: Vec<egui::Color32>,
    selected_plate: Option<usize>,
    
    // Feature toggles
    show_voronoi_plates: bool,
    show_boundary_evolution: bool,
    show_plate_motion: bool,
    enable_boundary_noise: bool,
    show_grid: bool,
    
    // Projection settings
    projection: ProjectionType,
    view_lon: f64,  // Longitude of the center of the orthographic view
    view_lat: f64,  // Latitude of the center of the orthographic view
}

impl TectonicApp {
    fn new() -> Self {
        let world = TectonicWorld::with_config(42, WorldConfig::default());
        
        // Generate colors for plates
        let plate_colors = (0..world.num_plates())
            .map(|i| {
                let hue = (i as f32 * 360.0 / world.num_plates() as f32) / 360.0;
                egui::Color32::from_rgb(
                    (hue * 255.0) as u8,
                    ((1.0 - hue) * 255.0) as u8,
                    ((hue * 0.5 + 0.5) * 255.0) as u8,
                )
            })
            .collect();
        
        Self {
            world,
            time: 0.0,
            dt: 0.1,
            is_running: false,
            step_count: 0,
            plate_colors,
            selected_plate: None,
            
            // Default feature toggles
            show_voronoi_plates: true,
            show_boundary_evolution: true,
            show_plate_motion: true,
            enable_boundary_noise: true,
            show_grid: true,
            
            // Default projection settings
            projection: ProjectionType::Orthographic,
            view_lon: 0.0,
            view_lat: 0.0,
        }
    }
    
    fn step_simulation(&mut self) {
        // Only step if plate motion is enabled
        if self.show_plate_motion {
            self.world.tick(self.dt, self.enable_boundary_noise);
            self.time += self.dt;
            self.step_count += 1;
        }
    }
    
    fn reset_simulation(&mut self) {
        self.world = TectonicWorld::new(6, 42);
        self.time = 0.0;
        self.step_count = 0;
        self.selected_plate = None;
    }
    
    fn spherical_to_screen(&self, point: SphericalPoint, screen_size: egui::Vec2) -> Option<egui::Pos2> {
        match self.projection {
            ProjectionType::Equirectangular => {
                // Normalize longitude to [0, 2π] range
                let lon = ((point.lon + PI) % (2.0 * PI) + (2.0 * PI)) % (2.0 * PI);
                let lat = point.lat;
                
                // Map to screen coordinates
                let x = (lon / (2.0 * PI)) * screen_size.x as f64;
                let y = ((PI/2.0 - lat) / PI) * screen_size.y as f64;
                
                // Clamp to screen bounds
                let x = x.max(0.0).min(screen_size.x as f64 - 1.0);
                let y = y.max(0.0).min(screen_size.y as f64 - 1.0);
                
                Some(egui::pos2(x as f32, y as f32))
            }
            ProjectionType::Orthographic => {
                // Orthographic projection: viewing a sphere from infinite distance
                // Center the view on (view_lat, view_lon)
                
                // Rotate point relative to view center
                let cos_c = self.view_lat.sin() * point.lat.sin() + 
                            self.view_lat.cos() * point.lat.cos() * (point.lon - self.view_lon).cos();
                
                // Don't draw points on the back of the sphere
                if cos_c < 0.0 {
                    return None;
                }
                
                // Project to 2D
                let x = point.lat.cos() * (point.lon - self.view_lon).sin();
                let y = self.view_lat.cos() * point.lat.sin() - 
                       self.view_lat.sin() * point.lat.cos() * (point.lon - self.view_lon).cos();
                
                // Scale to screen (using the smaller dimension to keep aspect ratio)
                let scale = screen_size.x.min(screen_size.y) as f64 * 0.45; // 0.45 to leave margin
                let center_x = screen_size.x as f64 * 0.5;
                let center_y = screen_size.y as f64 * 0.5;
                
                let screen_x = center_x + x * scale;
                let screen_y = center_y - y * scale; // Negative because screen y increases downward
                
                // Check if point is within screen bounds
                if screen_x >= 0.0 && screen_x < screen_size.x as f64 &&
                   screen_y >= 0.0 && screen_y < screen_size.y as f64 {
                    Some(egui::pos2(screen_x as f32, screen_y as f32))
                } else {
                    None
                }
            }
        }
    }
}

impl eframe::App for TectonicApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Auto-step if running
        if self.is_running {
            self.step_simulation();
            ctx.request_repaint();
        }
        
        // Right panel - Controls and stats
        egui::SidePanel::right("controls").resizable(true).default_width(300.0).show(ctx, |ui| {
            // Simulation Controls
            ui.heading("Simulation Controls");
            ui.horizontal(|ui| {
                if ui.button("Reset").clicked() {
                    self.reset_simulation();
                }
                
                if ui.button("Step").clicked() {
                    self.step_simulation();
                }
                
                if ui.button(if self.is_running { "Pause" } else { "Run" }).clicked() {
                    self.is_running = !self.is_running;
                }
            });
            
            ui.add_space(10.0);
            
            ui.horizontal(|ui| {
                ui.label("Time step:");
                ui.add(egui::Slider::new(&mut self.dt, 0.01..=1.0).text("dt"));
            });
            
            ui.add_space(10.0);
            ui.label(format!("Time: {:.2}", self.time));
            ui.label(format!("Steps: {}", self.step_count));
            
            if self.world.area_is_conserved() {
                ui.colored_label(egui::Color32::GREEN, "✅ Area conservation maintained");
            } else {
                ui.colored_label(egui::Color32::RED, "❌ Area conservation violated!");
            }
            
            ui.add_space(20.0);
            
            // Projection Controls
            ui.heading("Projection");
            ui.horizontal(|ui| {
                ui.radio_value(&mut self.projection, ProjectionType::Orthographic, "Orthographic (Globe)");
                ui.radio_value(&mut self.projection, ProjectionType::Equirectangular, "Equirectangular (Flat)");
            });
            
            if matches!(self.projection, ProjectionType::Orthographic) {
                ui.horizontal(|ui| {
                    ui.label("View Center:");
                    if ui.button("← Lon").clicked() { self.view_lon -= 0.3; }
                    if ui.button("Lon →").clicked() { self.view_lon += 0.3; }
                    if ui.button("↑ Lat").clicked() { self.view_lat = (self.view_lat + 0.3).min(PI/2.0); }
                    if ui.button("Lat ↓").clicked() { self.view_lat = (self.view_lat - 0.3).max(-PI/2.0); }
                });
                ui.label(format!("  Lat: {:.1}°, Lon: {:.1}°", 
                    self.view_lat.to_degrees(), self.view_lon.to_degrees()));
            }
            
            ui.add_space(20.0);
            
            // Feature Toggles
            ui.heading("Feature Toggles");
            ui.checkbox(&mut self.show_voronoi_plates, "Show Voronoi Plates");
            ui.checkbox(&mut self.show_boundary_evolution, "Enable Boundary Evolution");
            ui.checkbox(&mut self.show_plate_motion, "Enable Plate Motion");
            ui.checkbox(&mut self.enable_boundary_noise, "Enable Boundary Noise");
            ui.checkbox(&mut self.show_grid, "Show Coordinate Grid");
            
            ui.add_space(20.0);
            
            // Plate Statistics
            ui.heading("Plate Statistics");
            egui::ScrollArea::vertical().show(ui, |ui| {
                for (i, plate) in self.world.plates().iter().enumerate() {
                    let color = self.plate_colors[i];
                    let plate_type = match &plate.plate_type {
                        PlateType::Oceanic { density } => format!("Oceanic (ρ={:.2})", density),
                        PlateType::Continental { density } => format!("Continental (ρ={:.2})", density),
                    };
                    
                    ui.horizontal(|ui| {
                        // Color indicator
                        ui.colored_label(color, "●");
                        
                        // Plate info
                        ui.vertical(|ui| {
                            ui.label(format!("Plate {}: {}", i, plate_type));
                            ui.label(format!("Position: ({:.1}°, {:.1}°)", 
                                plate.seed.lat.to_degrees(), 
                                plate.seed.lon.to_degrees()));
                            ui.label(format!("Motion: ω={:.3} rad/time", plate.motion.angular_velocity));
                            ui.label(format!("Age: {:.1} time units", plate.age));
                        });
                        
                        // Select button
                        if ui.button("Select").clicked() {
                            self.selected_plate = Some(i);
                        }
                    });
                    
                    ui.add_space(5.0);
                }
            });
            
            ui.add_space(20.0);
            
            // Boundary Statistics
            ui.heading("Boundary Statistics");
            
            let mut boundary_counts = std::collections::HashMap::new();
            let mut total_complexity = 0.0;
            let mut total_age = 0.0;
            
            for boundary in self.world.boundaries() {
                *boundary_counts.entry(boundary.boundary_type.description()).or_insert(0) += 1;
                total_complexity += boundary.complexity();
                total_age += boundary.get_age();
            }
            
            ui.label(format!("Total boundaries: {}", self.world.boundaries().len()));
            ui.label(format!("Average complexity: {:.2}", 
                total_complexity / self.world.boundaries().len() as f64));
            ui.label(format!("Average age: {:.2}", 
                total_age / self.world.boundaries().len() as f64));
            
            // Debug: Show first few boundary geometries
            ui.add_space(10.0);
            ui.label("Debug - First 3 boundaries:");
            for (i, boundary) in self.world.boundaries().iter().take(3).enumerate() {
                ui.label(format!("  Boundary {}: {} points", i, boundary.geometry.len()));
                if !boundary.geometry.is_empty() {
                    let first_point = boundary.geometry[0];
                    ui.label(format!("    First point: ({:.1}°, {:.1}°)", 
                        first_point.lat.to_degrees(), first_point.lon.to_degrees()));
                }
            }
            
            ui.add_space(10.0);
            ui.label("Boundary types:");
            for (boundary_type, count) in boundary_counts {
                ui.label(format!("  {}: {}", boundary_type, count));
            }
        });

        // Main visualization area - fills remaining space
        egui::CentralPanel::default().show(ctx, |ui| {
            // Create a frame for the visualization that fills all available space
            let (rect, response) = ui.allocate_exact_size(
                ui.available_size(),
                egui::Sense::click(),
            );
            
            // Draw the visualization background
            ui.painter().rect_filled(rect, 0.0, egui::Color32::from_gray(20));
            
            // Create a painter that clips to the visualization area
            let painter = ui.painter().with_clip_rect(rect);
            
            // Draw coordinate grid (if enabled)
            if self.show_grid {
                let grid_color = egui::Color32::from_gray(40);
                let grid_stroke = egui::Stroke::new(1.0, grid_color);
                
                match self.projection {
                    ProjectionType::Equirectangular => {
                        // Draw longitude lines (vertical)
                        for i in 0..8 {
                            let x = (i as f32 / 7.0) * rect.width();
                            painter.add(egui::Shape::line_segment(
                                [egui::pos2(x, rect.top()), egui::pos2(x, rect.bottom())],
                                grid_stroke,
                            ));
                        }
                        
                        // Draw latitude lines (horizontal)
                        for i in 0..6 {
                            let y = (i as f32 / 5.0) * rect.height();
                            painter.add(egui::Shape::line_segment(
                                [egui::pos2(rect.left(), y), egui::pos2(rect.right(), y)],
                                grid_stroke,
                            ));
                        }
                    }
                    ProjectionType::Orthographic => {
                        // Draw sphere outline
                        let scale = rect.size().x.min(rect.size().y) * 0.45;
                        let center = egui::pos2(
                            rect.left() + rect.width() * 0.5,
                            rect.top() + rect.height() * 0.5
                        );
                        painter.circle_stroke(
                            center,
                            scale,
                            egui::Stroke::new(2.0, egui::Color32::from_gray(60))
                        );
                    }
                }
            }
            
            // Draw boundaries (if boundary evolution is enabled)
            if self.show_boundary_evolution {
                for boundary in self.world.boundaries() {
                    if boundary.geometry.len() < 2 {
                        continue;
                    }
                    
                    let color = match &boundary.boundary_type {
                        BoundaryType::Divergent { .. } => egui::Color32::from_rgb(0, 255, 0), // Green
                        BoundaryType::Convergent { .. } => egui::Color32::from_rgb(255, 0, 0), // Red
                        BoundaryType::Transform { .. } => egui::Color32::from_rgb(255, 255, 0), // Yellow
                    };
                    
                    // Draw boundary line segments
                    for i in 0..boundary.geometry.len() - 1 {
                        let point1 = boundary.geometry[i];
                        let point2 = boundary.geometry[i + 1];
                        
                        if let (Some(screen_pos1), Some(screen_pos2)) = (
                            self.spherical_to_screen(point1, rect.size()),
                            self.spherical_to_screen(point2, rect.size())
                        ) {
                            painter.add(egui::Shape::line_segment(
                                [screen_pos1, screen_pos2],
                                egui::Stroke::new(3.0, color),
                            ));
                        }
                    }
                }
            }
            
            // Draw plates (if enabled)
            if self.show_voronoi_plates {
                for plate_id in 0..self.world.num_plates() {
                    if let Some(plate) = self.world.get_plate(plate_id) {
                        if let Some(pos) = self.spherical_to_screen(plate.seed, rect.size()) {
                            let color = self.plate_colors[plate_id];
                            
                            // Draw plate center
                            let radius = if self.selected_plate == Some(plate_id) { 8.0 } else { 6.0 };
                            painter.circle_filled(pos, radius, color);
                            
                            // Draw plate ID
                            painter.text(
                                pos + egui::vec2(10.0, -10.0),
                                egui::Align2::LEFT_TOP,
                                format!("{}", plate_id),
                                egui::FontId::proportional(12.0),
                                egui::Color32::WHITE,
                            );
                        }
                    }
                }
            }
            
            
            // Handle clicks on the visualization
            if response.clicked() {
                self.selected_plate = None;
            }
        });
    }
}
