//! Basic usage example for the tectonics library
//! 
//! This example demonstrates Phase 1 functionality:
//! - Creating a tectonic world with multiple plates
//! - Querying which plate owns a given point
//! - Accessing plate properties

use tectonics::{TectonicWorld, SphericalPoint, WorldConfig};

fn main() {
    println!("🌍 Tectonic Plate Simulation - Phase 4 Demo");
    println!("===========================================");
    
    // Create a world with 8 tectonic plates using default config (random 20-40% continental)
    let mut world = TectonicWorld::new(8, 42);
    
    println!("Created world with {} plates and {} boundaries:", 
             world.num_plates(), world.boundaries().len());
    
    // Display initial state
    display_plates(&world);
    display_boundaries(&world);
    
    println!("\n🚀 Simulating plate motion over time:");
    
    // Simulate plate motion over multiple time steps
    let time_steps = 5;
    let dt = 0.1; // Time step
    
    for step in 0..time_steps {
        println!("\n--- Time Step {} ---", step);
        
        // Advance the simulation
        world.tick(dt, true);
        
        // Show how plates have moved
        println!("Plate positions after motion:");
        for plate in world.plates() {
            println!("  Plate {}: ({:.1}°, {:.1}°)", 
                     plate.id,
                     plate.seed.lat.to_degrees(),
                     plate.seed.lon.to_degrees());
        }
        
        // Check area conservation
        if world.area_is_conserved() {
            println!("✅ Area conservation maintained");
        } else {
            println!("❌ Area conservation violated!");
        }
        
        // Show boundary type changes
        let mut boundary_type_counts = std::collections::HashMap::new();
        for boundary in world.boundaries() {
            *boundary_type_counts.entry(boundary.boundary_type.description()).or_insert(0) += 1;
        }
        println!("Boundary types: {:?}", boundary_type_counts);
        
        // Show boundary evolution (Phase 4)
        if step == 0 {
            println!("Initial boundary complexity:");
            for (i, boundary) in world.boundaries().iter().take(3).enumerate() {
                println!("  Boundary {}: {} points, age: {:.1}, complexity: {:.2}", 
                         i, boundary.point_count(), boundary.get_age(), boundary.complexity());
            }
        } else if step == 4 {
            println!("Final boundary complexity:");
            for (i, boundary) in world.boundaries().iter().take(3).enumerate() {
                println!("  Boundary {}: {} points, age: {:.1}, complexity: {:.2}", 
                         i, boundary.point_count(), boundary.get_age(), boundary.complexity());
            }
        }
    }
    
    println!("\n🔧 Testing different continental fractions:");
    
    // Test with 30% continental plates
    let config_30 = WorldConfig::new().with_continental_fraction(0.3).with_number_of_plates(8);
    let world_30 = TectonicWorld::with_config(42, config_30);
    println!("\nWorld with 30% continental:");
    display_plates(&world_30);
    
    // Test with 70% continental plates
    let config_70 = WorldConfig::new().with_continental_fraction(0.7).with_number_of_plates(8);
    let world_70 = TectonicWorld::with_config(42, config_70);
    println!("\nWorld with 70% continental:");
    display_plates(&world_70);
    
    println!("\n🔍 Testing plate ownership queries:");
    
    // Test some specific points
    let test_points = vec![
        ("North Pole", SphericalPoint::new(std::f64::consts::PI/2.0, 0.0)),
        ("South Pole", SphericalPoint::new(-std::f64::consts::PI/2.0, 0.0)),
        ("Equator at 0°", SphericalPoint::new(0.0, 0.0)),
        ("Equator at 90°E", SphericalPoint::new(0.0, std::f64::consts::PI/2.0)),
        ("Equator at 180°", SphericalPoint::new(0.0, std::f64::consts::PI)),
    ];
    
    for (name, point) in test_points {
        let plate_id = world.plate_at(point);
        let plate = world.get_plate(plate_id).unwrap();
        
        println!("  {}: belongs to Plate {} ({})", 
                 name, 
                 plate_id,
                 match plate.plate_type {
                     tectonics::PlateType::Oceanic { .. } => "Oceanic",
                     tectonics::PlateType::Continental { .. } => "Continental",
                 });
    }
    
    println!("\n🔬 Testing boundary evolution over longer time:");
    
    // Create a new world for longer evolution test
    let mut evolution_world = TectonicWorld::new(6, 123);
    
    println!("Initial state:");
    for (i, boundary) in evolution_world.boundaries().iter().take(3).enumerate() {
        println!("  Boundary {}: {} points, complexity: {:.2}", 
                 i, boundary.point_count(), boundary.complexity());
    }
    
    // Simulate over longer time period
    for _ in 0..20 {
        evolution_world.tick(1.0, true);
    }
    
    println!("After 20 time steps:");
    for (i, boundary) in evolution_world.boundaries().iter().take(3).enumerate() {
        println!("  Boundary {}: {} points, age: {:.1}, complexity: {:.2}", 
                 i, boundary.point_count(), boundary.get_age(), boundary.complexity());
    }
    
    println!("\n✅ Phase 4 complete! Ready for Phase 5 (LOD System & Simplification)");
}

fn display_plates(world: &TectonicWorld) {
    let continental_count = world.plates().iter()
        .filter(|p| matches!(p.plate_type, tectonics::PlateType::Continental { .. }))
        .count();
    let oceanic_count = world.num_plates() - continental_count;
    
    println!("  {} continental, {} oceanic plates", continental_count, oceanic_count);
    
    for plate in world.plates() {
        let plate_type = match plate.plate_type {
            tectonics::PlateType::Oceanic { density } => format!("Oceanic (density: {:.2})", density),
            tectonics::PlateType::Continental { density } => format!("Continental (density: {:.2})", density),
        };
        
        println!("    Plate {}: {} at ({:.1}°, {:.1}°)", 
                 plate.id, 
                 plate_type,
                 plate.seed.lat.to_degrees(),
                 plate.seed.lon.to_degrees());
    }
}

fn display_boundaries(world: &TectonicWorld) {
    println!("\nBoundary network:");
    
    let mut boundary_counts = std::collections::HashMap::new();
    let mut boundary_types = std::collections::HashMap::new();
    
    for boundary in world.boundaries() {
        // Count boundaries per plate
        *boundary_counts.entry(boundary.plate_a).or_insert(0) += 1;
        *boundary_counts.entry(boundary.plate_b).or_insert(0) += 1;
        
        // Count boundary types
        *boundary_types.entry(boundary.boundary_type.description()).or_insert(0) += 1;
    }
    
    println!("  Boundary types: {:?}", boundary_types);
    
    println!("  Plate connectivity:");
    for plate in world.plates() {
        let count = boundary_counts.get(&plate.id).unwrap_or(&0);
        println!("    Plate {}: {} boundaries", plate.id, count);
    }
    
    // Show some example boundaries
    println!("  Example boundaries:");
    for (i, boundary) in world.boundaries().iter().take(3).enumerate() {
        let plate_a_type = match world.get_plate(boundary.plate_a).unwrap().plate_type {
            tectonics::PlateType::Oceanic { .. } => "Oceanic",
            tectonics::PlateType::Continental { .. } => "Continental",
        };
        let plate_b_type = match world.get_plate(boundary.plate_b).unwrap().plate_type {
            tectonics::PlateType::Oceanic { .. } => "Oceanic",
            tectonics::PlateType::Continental { .. } => "Continental",
        };
        
        println!("    Boundary {}: {} {} <-> {} {} ({})", 
                 i,
                 plate_a_type, boundary.plate_a,
                 plate_b_type, boundary.plate_b,
                 boundary.boundary_type.description());
    }
}
