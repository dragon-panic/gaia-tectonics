//! Tests for boundary connectivity and integrity during evolution
//! 
//! These tests ensure that plate boundaries maintain proper connectivity
//! as the simulation evolves, preventing gaps and discontinuities.

use tectonics::*;

#[test]
fn boundaries_remain_continuous_after_evolution() {
    let mut world = TectonicWorld::new(6, 42);
    
    // Evolve the simulation for several steps
    for _ in 0..50 {
        world.tick(0.1, true);
        
        // Check that each boundary is continuous (no gaps between adjacent points)
        for boundary in world.boundaries() {
            if boundary.geometry.len() < 2 {
                continue;
            }
            
            let mut max_segment_length: f64 = 0.0;
            for i in 0..boundary.geometry.len() - 1 {
                let seg_length = distance(boundary.geometry[i], boundary.geometry[i + 1]);
                max_segment_length = max_segment_length.max(seg_length);
            }
            
            // No segment should be longer than a reasonable threshold
            // Transform boundaries can subdivide to 150km segments, which is ~0.023 radians
            // Allow some margin for noise displacement
            assert!(
                max_segment_length < 0.015,
                "Boundary between plates {} and {} has a segment of length {} radians, which is too large. \
                 This indicates a gap or discontinuity.",
                boundary.plate_a, boundary.plate_b, max_segment_length
            );
        }
    }
}

#[test]
fn triple_junctions_remain_connected() {
    let mut world = TectonicWorld::new(6, 42);
    
    // Evolve the simulation
    for step in 0..20 {
        world.tick(0.1, true);
        
        // Find all boundary endpoints and check for triple junctions
        let mut endpoints: Vec<(usize, bool, SphericalPoint)> = Vec::new(); // (boundary_index, is_start, position)
        
        for (i, boundary) in world.boundaries().iter().enumerate() {
            if boundary.geometry.is_empty() {
                continue;
            }
            endpoints.push((i, true, boundary.geometry[0]));
            endpoints.push((i, false, *boundary.geometry.last().unwrap()));
        }
        
        // Check that endpoints that should be connected (triple junctions) are actually close
        for (idx1, is_start1, pos1) in &endpoints {
            let mut group = vec![(*idx1, *is_start1)];
            
            for (idx2, is_start2, pos2) in &endpoints {
                if idx1 == idx2 && is_start1 == is_start2 {
                    continue; // Same endpoint
                }
                
                let dist = distance(*pos1, *pos2);
                if dist < 0.001 { // Within ~6km on Earth
                    group.push((*idx2, *is_start2));
                }
            }
            
            if group.len() >= 3 {
                // This is a triple junction - verify all points are close
                for (idx, is_start) in &group {
                    let boundary = &world.boundaries()[*idx];
                    let endpoint_pos = if *is_start {
                        boundary.geometry[0]
                    } else {
                        *boundary.geometry.last().unwrap()
                    };
                    
                    let dist = distance(*pos1, endpoint_pos);
                    assert!(
                        dist < 0.01,
                        "Step {}: Triple junction not maintained! Endpoints of boundaries should be co-located but are {} radians apart",
                        step, dist
                    );
                }
            }
        }
    }
}

#[test]
fn boundary_count_remains_stable() {
    let mut world = TectonicWorld::new(8, 42);
    
    let initial_boundary_count = world.boundaries().len();
    
    // Evolve the simulation
    for _ in 0..30 {
        world.tick(0.1, true);
        
        let current_boundary_count = world.boundaries().len();
        assert_eq!(
            initial_boundary_count, current_boundary_count,
            "Boundary count should remain constant during evolution. \
             Started with {} boundaries, now have {}",
            initial_boundary_count, current_boundary_count
        );
    }
}

#[test]
fn no_boundary_segments_exceed_max_length() {
    let mut world = TectonicWorld::new(6, 42);
    
    // Evolve the simulation for many steps
    for step in 0..100 {
        world.tick(0.1, true);
        
        // Check segment lengths in all boundaries
        for boundary in world.boundaries() {
            if boundary.geometry.len() < 2 {
                continue;
            }
            
            for i in 0..boundary.geometry.len() - 1 {
                let seg_length = distance(boundary.geometry[i], boundary.geometry[i + 1]);
                
                // Even with evolution, segments should be subdivided if they get too long
                // Maximum should be around the max_segment_length threshold
                assert!(
                    seg_length < 0.02,
                    "Step {}: Segment {} in boundary {}-{} has length {} radians, which exceeds maximum. \
                     Evolution should subdivide long segments.",
                    step, i, boundary.plate_a, boundary.plate_b, seg_length
                );
            }
        }
    }
}

#[test]
fn boundary_endpoints_match_voronoi_vertices() {
    let mut world = TectonicWorld::new(6, 42);
    
    // Evolve a bit
    for _ in 0..10 {
        world.tick(0.1, true);
    }
    
    // Get the Voronoi vertices from the current configuration
    let seeds: Vec<SphericalPoint> = world.plates().iter().map(|p| p.seed).collect();
    let voronoi_diagram = SphericalVoronoiDiagram::new(seeds);
    
    // Check that boundary endpoints are near Voronoi vertices
    for boundary in world.boundaries() {
        if boundary.geometry.is_empty() {
            continue;
        }
        
        let start = boundary.geometry[0];
        let end = *boundary.geometry.last().unwrap();
        
        // Find the closest Voronoi vertex to each endpoint
        let mut start_has_close_vertex = false;
        let mut end_has_close_vertex = false;
        
        for vertex in &voronoi_diagram.vertices {
            if distance(start, vertex.position) < 0.01 {
                start_has_close_vertex = true;
            }
            if distance(end, vertex.position) < 0.01 {
                end_has_close_vertex = true;
            }
        }
        
        assert!(
            start_has_close_vertex,
            "Boundary start point should be near a Voronoi vertex. \
             This indicates the boundary endpoints are not properly constrained."
        );
        assert!(
            end_has_close_vertex,
            "Boundary end point should be near a Voronoi vertex. \
             This indicates the boundary endpoints are not properly constrained."
        );
    }
}

#[test]
fn evolution_preserves_topology() {
    let mut world = TectonicWorld::new(6, 42);
    
    // Record initial adjacency relationships
    let mut initial_adjacencies: std::collections::HashSet<(PlateId, PlateId)> = 
        std::collections::HashSet::new();
    
    for boundary in world.boundaries() {
        let pair = if boundary.plate_a < boundary.plate_b {
            (boundary.plate_a, boundary.plate_b)
        } else {
            (boundary.plate_b, boundary.plate_a)
        };
        initial_adjacencies.insert(pair);
    }
    
    // Evolve the simulation
    for _ in 0..30 {
        world.tick(0.1, true);
    }
    
    // Check that the same plate pairs are adjacent
    let mut current_adjacencies: std::collections::HashSet<(PlateId, PlateId)> = 
        std::collections::HashSet::new();
    
    for boundary in world.boundaries() {
        let pair = if boundary.plate_a < boundary.plate_b {
            (boundary.plate_a, boundary.plate_b)
        } else {
            (boundary.plate_b, boundary.plate_a)
        };
        current_adjacencies.insert(pair);
    }
    
    assert_eq!(
        initial_adjacencies, current_adjacencies,
        "The topology (which plates are adjacent) should remain constant during evolution"
    );
}

