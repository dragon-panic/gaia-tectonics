//! Comprehensive tests for spherical Voronoi diagram implementation

use gaia_tectonics::*;
use std::f64::consts::PI;

const EPSILON: f64 = 1e-6;

#[test]
fn test_tetrahedron_voronoi() {
    // 4 points forming a regular tetrahedron
    // This is the simplest non-degenerate case
    let seeds = vec![
        SphericalPoint::new(0.0, 0.0),                    // North pole region
        SphericalPoint::new(0.0, 2.0 * PI / 3.0),         // 120° around equator
        SphericalPoint::new(0.0, 4.0 * PI / 3.0),         // 240° around equator
        SphericalPoint::new(PI / 2.0, 0.0),               // Equator
    ];
    
    let diagram = SphericalVoronoiDiagram::new(seeds);
    
    println!("Tetrahedron: {} vertices, {} edges", diagram.vertices.len(), diagram.edges.len());
    
    // Euler's formula for sphere: V - E + F = 2
    // For tetrahedron: 4 faces, expected 6 edges, 4 vertices
    assert!(diagram.vertices.len() >= 4, "Should have at least 4 vertices");
    assert!(diagram.edges.len() >= 6, "Should have at least 6 edges (one per pair of seeds)");
    
    // Every seed should have edges
    for i in 0..4 {
        let has_edge = diagram.edges.iter().any(|e| e.plate_a == i || e.plate_b == i);
        assert!(has_edge, "Seed {} should have at least one edge", i);
    }
    
    // Check edge connectivity
    verify_edge_connectivity(&diagram);
}

#[test]
fn test_five_seeds_regular_distribution() {
    // 5 seeds in a somewhat regular pattern
    let seeds = vec![
        SphericalPoint::new(PI / 4.0, 0.0),
        SphericalPoint::new(PI / 4.0, PI / 2.0),
        SphericalPoint::new(PI / 4.0, PI),
        SphericalPoint::new(PI / 4.0, 3.0 * PI / 2.0),
        SphericalPoint::new(-PI / 4.0, PI / 4.0),
    ];
    
    let diagram = SphericalVoronoiDiagram::new(seeds);
    
    println!("5 seeds: {} vertices, {} edges", diagram.vertices.len(), diagram.edges.len());
    
    // Should have edges
    assert!(diagram.edges.len() >= 5, "Should have at least 5 edges");
    
    // Every seed should participate
    for i in 0..5 {
        let edge_count = diagram.edges.iter()
            .filter(|e| e.plate_a == i || e.plate_b == i)
            .count();
        assert!(edge_count >= 2, "Seed {} should have at least 2 edges, has {}", i, edge_count);
    }
    
    verify_edge_connectivity(&diagram);
}

#[test]
fn test_clustered_seeds() {
    // Seeds that are close together (stress test)
    let seeds = vec![
        SphericalPoint::new(0.0, 0.0),
        SphericalPoint::new(0.1, 0.0),
        SphericalPoint::new(0.0, 0.1),
        SphericalPoint::new(0.1, 0.1),
        SphericalPoint::new(1.0, 1.0), // One far away
    ];
    
    let diagram = SphericalVoronoiDiagram::new(seeds);
    
    println!("Clustered: {} vertices, {} edges", diagram.vertices.len(), diagram.edges.len());
    
    // Should still create valid tessellation
    assert!(!diagram.edges.is_empty(), "Should have edges even with clustered seeds");
    
    // The far seed should definitely have edges
    let far_seed_edges = diagram.edges.iter()
        .filter(|e| e.plate_a == 4 || e.plate_b == 4)
        .count();
    assert!(far_seed_edges >= 1, "Far seed should have edges");
}

#[test]
fn test_antipodal_seeds() {
    // Seeds that are opposite on the sphere
    let seeds = vec![
        SphericalPoint::new(PI / 4.0, 0.0),
        SphericalPoint::new(-PI / 4.0, PI),      // Opposite side
        SphericalPoint::new(0.0, PI / 2.0),
        SphericalPoint::new(0.0, 3.0 * PI / 2.0),
    ];
    
    let diagram = SphericalVoronoiDiagram::new(seeds);
    
    println!("Antipodal: {} vertices, {} edges", diagram.vertices.len(), diagram.edges.len());
    
    verify_edge_connectivity(&diagram);
}

#[test]
fn test_all_edges_have_valid_geometry() {
    let seeds = vec![
        SphericalPoint::new(0.0, 0.0),
        SphericalPoint::new(0.5, 0.0),
        SphericalPoint::new(0.0, 0.5),
        SphericalPoint::new(0.5, 0.5),
        SphericalPoint::new(0.25, 0.25),
        SphericalPoint::new(-0.25, 0.25),
    ];
    
    let diagram = SphericalVoronoiDiagram::new(seeds);
    
    for (i, edge) in diagram.edges.iter().enumerate() {
        // Edge should have geometry
        assert!(!edge.geometry.is_empty(), "Edge {} has no geometry", i);
        assert!(edge.geometry.len() >= 2, "Edge {} has only {} points", i, edge.geometry.len());
        
        // Start and end should be different
        let start = edge.geometry[0];
        let end = *edge.geometry.last().unwrap();
        let edge_length = distance(start, end);
        
        assert!(edge_length > EPSILON, 
                "Edge {} (plates {}-{}) is degenerate: length = {}", 
                i, edge.plate_a, edge.plate_b, edge_length);
        
        // Geometry should be continuous
        for j in 0..edge.geometry.len() - 1 {
            let seg_len = distance(edge.geometry[j], edge.geometry[j + 1]);
            assert!(seg_len > 0.0 && seg_len < 0.1, 
                    "Edge {} segment {} has invalid length: {}", i, j, seg_len);
        }
    }
}

#[test]
fn test_vertices_are_equidistant() {
    let seeds = vec![
        SphericalPoint::new(0.0, 0.0),
        SphericalPoint::new(0.5, 0.0),
        SphericalPoint::new(0.0, 0.5),
        SphericalPoint::new(0.5, 0.5),
    ];
    
    let diagram = SphericalVoronoiDiagram::new(seeds.clone());
    
    // Each Voronoi vertex should be equidistant to its incident seeds
    for (i, vertex) in diagram.vertices.iter().enumerate() {
        if vertex.incident_seeds.len() < 3 {
            println!("Warning: Vertex {} only has {} incident seeds", 
                     i, vertex.incident_seeds.len());
            continue;
        }
        
        let distances: Vec<f64> = vertex.incident_seeds.iter()
            .map(|&seed_idx| distance(vertex.position, seeds[seed_idx]))
            .collect();
        
        let first_dist = distances[0];
        for (j, &dist) in distances.iter().enumerate() {
            assert!((dist - first_dist).abs() < 0.01,
                    "Vertex {} not equidistant: seed {} dist={:.6}, seed 0 dist={:.6}",
                    i, j, dist, first_dist);
        }
    }
}

#[test]
fn test_every_pair_has_at_most_one_edge() {
    let seeds = vec![
        SphericalPoint::new(0.0, 0.0),
        SphericalPoint::new(0.5, 0.0),
        SphericalPoint::new(0.0, 0.5),
        SphericalPoint::new(0.5, 0.5),
        SphericalPoint::new(0.25, 0.25),
    ];
    
    let diagram = SphericalVoronoiDiagram::new(seeds.clone());
    
    // Check that each pair of seeds has at most one edge
    for i in 0..seeds.len() {
        for j in (i + 1)..seeds.len() {
            let matching_edges: Vec<_> = diagram.edges.iter()
                .filter(|e| {
                    (e.plate_a == i && e.plate_b == j) ||
                    (e.plate_a == j && e.plate_b == i)
                })
                .collect();
            
            assert!(matching_edges.len() <= 1,
                    "Seeds {} and {} have {} edges (should be 0 or 1)",
                    i, j, matching_edges.len());
        }
    }
}

/// Verify that edges form proper connectivity (triple junctions)
fn verify_edge_connectivity(diagram: &SphericalVoronoiDiagram) {
    println!("\nVerifying edge connectivity:");
    println!("  Total edges: {}", diagram.edges.len());
    
    let mut isolated_endpoints = Vec::new();
    
    for (i, edge) in diagram.edges.iter().enumerate() {
        if edge.geometry.is_empty() {
            println!("  WARNING: Edge {} has empty geometry", i);
            continue;
        }
        
        let start = edge.geometry[0];
        let end = *edge.geometry.last().unwrap();
        
        // Count connections at start
        let start_connections = count_endpoint_connections(&diagram.edges, start, i);
        
        // Count connections at end
        let end_connections = count_endpoint_connections(&diagram.edges, end, i);
        
        if start_connections < 2 {
            println!("  WARNING: Edge {} start only has {} connections", i, start_connections);
            isolated_endpoints.push((i, true, start, start_connections));
        }
        
        if end_connections < 2 {
            println!("  WARNING: Edge {} end only has {} connections", i, end_connections);
            isolated_endpoints.push((i, false, end, end_connections));
        }
    }
    
    if !isolated_endpoints.is_empty() {
        println!("\nIsolated endpoints found:");
        for (edge_idx, is_start, pos, connections) in &isolated_endpoints {
            println!("  Edge {} {}: ({:.6}, {:.6}) - {} connections",
                     edge_idx,
                     if *is_start { "start" } else { "end" },
                     pos.lat, pos.lon,
                     connections);
        }
        
        // This is actually a bug - fail the test
        panic!("{} isolated endpoints found (see output above)", isolated_endpoints.len());
    }
    
    println!("  ✓ All endpoints properly connected");
}

fn count_endpoint_connections(edges: &[VoronoiEdge], point: SphericalPoint, skip_edge: usize) -> usize {
    let mut count = 0;
    
    for (i, edge) in edges.iter().enumerate() {
        if i == skip_edge || edge.geometry.is_empty() {
            continue;
        }
        
        let other_start = edge.geometry[0];
        let other_end = *edge.geometry.last().unwrap();
        
        if distance(point, other_start) < 0.001 || distance(point, other_end) < 0.001 {
            count += 1;
        }
    }
    
    count
}

