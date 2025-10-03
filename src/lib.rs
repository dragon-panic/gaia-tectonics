//! # Tectonics - Spherical Tectonic Plate Simulation
//! 
//! A Voronoi-based tectonic plate simulation library designed for realistic geological 
//! modeling on spherical surfaces. This library provides:
//! 
//! - **Spherical Voronoi tessellation** for plate boundary generation
//! - **Physics-based plate motion** using Euler pole rotation
//! - **Dynamic boundary evolution** with fractal detail
//! - **Realistic plate interactions** (convergent, divergent, transform)
//! - **Deterministic simulation** for reproducible results
//! 
//! ## Features
//! 
//! - 🌍 **Spherical geometry**: All calculations work on a unit sphere
//! - 🔄 **Plate dynamics**: Plates rotate according to Euler pole mechanics
//! - 🏔️ **Boundary classification**: Automatic detection of convergent/divergent/transform boundaries
//! - 📈 **Evolving complexity**: Boundaries gain fractal detail over geological time
//! - 🎯 **Efficient queries**: Fast plate ownership and boundary lookups
//! - 🔬 **Physics validation**: Area conservation and motion accuracy
//! 
//! ## Quick Start
//! 
//! Add this to your `Cargo.toml`:
//! 
//! ```toml
//! [dependencies]
//! tectonics = "0.1"
//! ```
//! 
//! ## Basic Usage
//! 
//! ```rust
//! use tectonics::{TectonicWorld, SphericalPoint};
//! 
//! // Create a world with 12 plates using seed 42 for deterministic generation
//! let mut world = TectonicWorld::new(12, 42);
//! 
//! // Query which plate owns a point
//! let point = SphericalPoint::new(0.0, 0.0); // Equator, Prime Meridian
//! let plate_id = world.plate_at(point);
//! println!("Point is on plate {}", plate_id);
//! 
//! // Get information about a plate
//! if let Some(plate) = world.get_plate(plate_id) {
//!     println!("Plate type: {:?}", plate.plate_type);
//!     println!("Plate age: {} years", plate.age);
//! }
//! 
//! // Advance the simulation by 1 million years
//! world.tick(1e6, true);
//! 
//! // Examine boundaries
//! for boundary in world.boundaries() {
//!     println!("Boundary between plates {} and {}", 
//!              boundary.plate_a, boundary.plate_b);
//!     println!("Type: {:?}", boundary.boundary_type);
//!     println!("Length: {:.4} radians", boundary.length());
//! }
//! ```
//! 
//! ## Advanced Configuration
//! 
//! ```rust
//! use tectonics::{TectonicWorld, WorldConfig};
//! 
//! // Create a custom world with 80% continental plates
//! let config = WorldConfig::new()
//!     .with_number_of_plates(15)
//!     .with_continental_fraction(0.8);
//! 
//! let world = TectonicWorld::with_config(42, config);
//! ```
//! 
//! ## Boundary Types
//! 
//! The library automatically classifies boundaries based on plate motions:
//! 
//! - **Divergent**: Plates moving apart (spreading centers, mid-ocean ridges)
//! - **Convergent**: Plates moving together (subduction zones, mountain building)
//! - **Transform**: Plates sliding past each other (strike-slip faults)
//! 
//! ## Architecture
//! 
//! This library focuses purely on tectonic plate dynamics. For terrain generation
//! and heightmap creation, consider using the companion `tectonics-terrain` crate.
//! 
//! ## Modules
//! 
//! - [`geometry`]: Core geometric types and spherical math
//! - [`world`]: Main simulation container and world-level operations
//! - [`boundary`]: Plate boundary types and evolution
//! - [`constants`]: Physical constants and simulation parameters
//! 
//! ## Examples
//! 
//! See the `examples/` directory for complete examples:
//! 
//! - `basic_usage.rs` - Simple command-line example
//! - `simple_gui.rs` - Interactive visualization (requires `--features gui`)
//! 
//! ## References
//! 
//! This library implements concepts from:
//! - Spherical Voronoi tessellation
//! - Euler pole rotation mechanics
//! - Plate tectonic theory

pub mod geometry;
pub mod world;
pub mod boundary;
pub mod constants;

// Re-export core types for convenience
pub use geometry::*;
pub use world::*;
pub use boundary::*;
pub use constants::*;

