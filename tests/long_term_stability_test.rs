//! Long-term stability tests for plate tessellation validity
//! 
//! These tests run the simulation for many steps to ensure that
//! the tessellation remains valid and boundaries stay properly connected

use tectonics::*;

#[test]
fn plates_remain_valid_over_thousand_steps() {
    let mut world = TectonicWorld::new(8, 42);
    
    for step in 0..1000 {
        world.tick(0.1, true);
        
        // Run full validation every 100 steps to catch issues
        if step % 100 == 0 || step < 10 {
            validate_tessellation(&world, step);
        }
    }
}

#[test]
fn plates_remain_valid_with_fast_motion() {
    let mut world = TectonicWorld::new(6, 42);
    
    // Use larger timesteps to stress-test the system
    for step in 0..500 {
        world.tick(0.5, true); // Larger dt
        
        if step % 50 == 0 {
            validate_tessellation(&world, step);
        }
    }
}

#[test]
fn many_plates_remain_valid() {
    let mut world = TectonicWorld::new(15, 42);
    
    for step in 0..500 {
        world.tick(0.1, true);
        
        if step % 100 == 0 {
            validate_tessellation(&world, step);
        }
    }
}

/// Comprehensive validation of tessellation properties
fn validate_tessellation(world: &TectonicWorld, step: usize) {
    println!("Validating tessellation at step {}", step);
    
    // 1. Check that all boundaries are continuous (no gaps)
    validate_boundary_continuity(world, step);
    
    // 2. Check that triple junctions are properly connected
    validate_triple_junctions(world, step);
    
    // 3. Check that each plate has a closed boundary network
    validate_plate_closure(world, step);
    
    // 4. Check that boundaries don't overlap
    validate_no_overlaps(world, step);
    
    // 5. Check that the tessellation covers the whole sphere
    validate_complete_coverage(world, step);
}

fn validate_boundary_continuity(world: &TectonicWorld, step: usize) {
    for (i, boundary) in world.boundaries().iter().enumerate() {
        if boundary.geometry.len() < 2 {
            continue;
        }
        
        let mut max_gap: f64 = 0.0;
        for j in 0..boundary.geometry.len() - 1 {
            let gap = distance(boundary.geometry[j], boundary.geometry[j + 1]);
            max_gap = max_gap.max(gap);
        }
        
        // No gap should exceed ~100km (0.015 radians)
        assert!(
            max_gap < 0.02,
            "Step {}: Boundary {} has gap of {} radians (plates {}-{})",
            step, i, max_gap, boundary.plate_a, boundary.plate_b
        );
    }
}

fn validate_triple_junctions(world: &TectonicWorld, step: usize) {
    // Collect all boundary endpoints
    let mut endpoints: Vec<(PlateId, PlateId, bool, SphericalPoint)> = Vec::new();
    
    for boundary in world.boundaries() {
        if boundary.geometry.is_empty() {
            continue;
        }
        endpoints.push((boundary.plate_a, boundary.plate_b, true, boundary.geometry[0]));
        endpoints.push((boundary.plate_a, boundary.plate_b, false, *boundary.geometry.last().unwrap()));
    }
    
    // Check that endpoints that should be connected are actually close
    for i in 0..endpoints.len() {
        let mut nearby_count = 0;
        let (plate_a1, plate_b1, _, pos1) = endpoints[i];
        
        for (j, &endpoint) in endpoints.iter().enumerate() {
            if i == j {
                continue;
            }
            
            let (plate_a2, plate_b2, _, pos2) = endpoint;
            let dist = distance(pos1, pos2);
            
            // If endpoints are very close, they should share a plate
            if dist < 0.001 { // Within ~6km
                nearby_count += 1;
                
                let shares_plate = plate_a1 == plate_a2 || plate_a1 == plate_b2 ||
                                  plate_b1 == plate_a2 || plate_b1 == plate_b2;
                
                assert!(
                    shares_plate,
                    "Step {}: Endpoints at ({:.6}, {:.6}) are close but don't share plates: ({},{}) vs ({},{})",
                    step, pos1.lat, pos1.lon, plate_a1, plate_b1, plate_a2, plate_b2
                );
            }
        }
        
        // Each endpoint should connect to at least 2 other boundaries (forming a triple junction)
        assert!(
            nearby_count >= 2,
            "Step {}: Endpoint at ({:.6}, {:.6}) only connects to {} other boundaries (expected >= 2)",
            step, pos1.lat, pos1.lon, nearby_count
        );
    }
}

fn validate_plate_closure(world: &TectonicWorld, step: usize) {
    // For each plate, verify it has a connected network of boundaries
    for plate in world.plates() {
        let plate_boundaries: Vec<_> = world.boundaries().iter()
            .filter(|b| b.involves_plate(plate.id))
            .collect();
        
        assert!(
            !plate_boundaries.is_empty(),
            "Step {}: Plate {} has no boundaries",
            step, plate.id
        );
        
        // Each plate should have at least 3 boundaries to form a closed region
        assert!(
            plate_boundaries.len() >= 3,
            "Step {}: Plate {} only has {} boundaries (need >= 3 for closure)",
            step, plate.id, plate_boundaries.len()
        );
    }
}

fn validate_no_overlaps(world: &TectonicWorld, step: usize) {
    let boundaries: Vec<_> = world.boundaries().iter().collect();
    
    // Check each pair of non-adjacent boundaries
    for i in 0..boundaries.len() {
        for j in (i + 1)..boundaries.len() {
            let b1 = boundaries[i];
            let b2 = boundaries[j];
            
            // Skip adjacent boundaries (they share plates)
            if b1.plate_a == b2.plate_a || b1.plate_a == b2.plate_b ||
               b1.plate_b == b2.plate_a || b1.plate_b == b2.plate_b {
                continue;
            }
            
            // Check if any points are too close (indicating overlap)
            for p1 in &b1.geometry {
                for p2 in &b2.geometry {
                    let dist = distance(*p1, *p2);
                    
                    // Non-adjacent boundaries should not have points within 10km
                    assert!(
                        dist > 0.0015,
                        "Step {}: Non-adjacent boundaries {} and {} have overlapping points (dist: {})",
                        step, i, j, dist
                    );
                }
            }
        }
    }
}

fn validate_complete_coverage(world: &TectonicWorld, step: usize) {
    // Sample random points on the sphere and verify each belongs to exactly one plate
    let sample_count = 100;
    let mut rng = rand::thread_rng();
    
    use rand::Rng;
    for _ in 0..sample_count {
        let lat = (rng.gen::<f64>() - 0.5) * std::f64::consts::PI;
        let lon = rng.gen::<f64>() * 2.0 * std::f64::consts::PI;
        let point = SphericalPoint::new(lat, lon);
        
        // This should succeed - every point should have a plate
        let plate_id = world.plate_at(point);
        
        // Verify this is a valid plate
        assert!(
            plate_id < world.plates().len(),
            "Step {}: Invalid plate ID {} for point ({:.6}, {:.6})",
            step, plate_id, lat, lon
        );
    }
}

