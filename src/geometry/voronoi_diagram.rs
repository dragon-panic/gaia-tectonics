//! Complete spherical Voronoi diagram implementation
//! 
//! This module computes true Voronoi diagrams on a unit sphere, including:
//! - Voronoi vertices (points equidistant to 3+ seeds)
//! - Voronoi edges (great circle arcs between vertices)
//! - Proper adjacency detection

use super::point::{SphericalPoint, PlateId, distance};
use crate::constants::{EPSILON, DEDUPLICATION_EPSILON, TARGET_SEGMENT_LENGTH_RADIANS};

/// A Voronoi vertex - a point equidistant to 3 or more seeds
#[derive(Debug, Clone)]
pub struct VoronoiVertex {
    pub position: SphericalPoint,
    /// The 3+ seeds this vertex is equidistant to
    pub incident_seeds: Vec<PlateId>,
}

/// A Voronoi edge - an arc on a great circle between two vertices
#[derive(Debug, Clone)]
pub struct VoronoiEdge {
    pub plate_a: PlateId,
    pub plate_b: PlateId,
    pub vertices: (usize, usize), // Indices into the vertex list
    pub geometry: Vec<SphericalPoint>,
}

/// Complete spherical Voronoi diagram
#[derive(Debug, Clone)]
pub struct SphericalVoronoiDiagram {
    pub seeds: Vec<SphericalPoint>,
    pub vertices: Vec<VoronoiVertex>,
    pub edges: Vec<VoronoiEdge>,
}

impl SphericalVoronoiDiagram {
    /// Compute the complete Voronoi diagram for the given seeds
    pub fn new(seeds: Vec<SphericalPoint>) -> Self {
        if seeds.len() < 4 {
            // Need at least 4 points for a non-degenerate spherical Voronoi diagram
            return Self {
                seeds,
                vertices: Vec::new(),
                edges: Vec::new(),
            };
        }

        let vertices = Self::compute_vertices(&seeds);
        let edges = Self::compute_edges(&seeds, &vertices);

        Self {
            seeds,
            vertices,
            edges,
        }
    }

    /// Find all Voronoi vertices (points equidistant to 3+ seeds)
    fn compute_vertices(seeds: &[SphericalPoint]) -> Vec<VoronoiVertex> {
        let mut vertices = Vec::new();
        let n = seeds.len();

        // Check all triples of seeds
        for i in 0..n {
            for j in (i + 1)..n {
                for k in (j + 1)..n {
                    // Find the two points equidistant to seeds i, j, k
                    if let Some((p1, p2)) = Self::circumcenter_triple(seeds[i], seeds[j], seeds[k]) {
                        // Check if p1 is a valid Voronoi vertex
                        if Self::is_voronoi_vertex(p1, &[i, j, k], seeds) {
                            let mut incident = vec![i, j, k];
                            // Check for additional equidistant seeds
                            let dist = distance(p1, seeds[i]);
                            for (m, seed) in seeds.iter().enumerate().take(n) {
                                if m != i && m != j && m != k
                                    && (distance(p1, *seed) - dist).abs() < EPSILON {
                                    incident.push(m);
                                }
                            }
                            vertices.push(VoronoiVertex {
                                position: p1,
                                incident_seeds: incident,
                            });
                        }

                        // Check if p2 is a valid Voronoi vertex (p2 is antipodal to p1)
                        if Self::is_voronoi_vertex(p2, &[i, j, k], seeds) {
                            let mut incident = vec![i, j, k];
                            let dist = distance(p2, seeds[i]);
                            for (m, seed) in seeds.iter().enumerate().take(n) {
                                if m != i && m != j && m != k
                                    && (distance(p2, *seed) - dist).abs() < EPSILON {
                                    incident.push(m);
                                }
                            }
                            vertices.push(VoronoiVertex {
                                position: p2,
                                incident_seeds: incident,
                            });
                        }
                    }
                }
            }
        }

        // Remove duplicate vertices
        Self::deduplicate_vertices(vertices)
    }

    /// Find the two points equidistant to three seeds on a sphere
    /// These lie at the intersection of the perpendicular bisector planes
    fn circumcenter_triple(a: SphericalPoint, b: SphericalPoint, c: SphericalPoint) -> Option<(SphericalPoint, SphericalPoint)> {
        // Convert to Cartesian
        let cart_a = a.to_cartesian();
        let cart_b = b.to_cartesian();
        let cart_c = c.to_cartesian();

        // The perpendicular bisector plane of a and b has normal (b - a)
        let normal_ab = [
            cart_b[0] - cart_a[0],
            cart_b[1] - cart_a[1],
            cart_b[2] - cart_a[2],
        ];

        // The perpendicular bisector plane of a and c has normal (c - a)
        let normal_ac = [
            cart_c[0] - cart_a[0],
            cart_c[1] - cart_a[1],
            cart_c[2] - cart_a[2],
        ];

        // The intersection of these two planes (through origin) is a line
        // The direction of this line is perpendicular to both normals
        let direction = Self::cross_product(normal_ab, normal_ac);
        
        // Normalize the direction
        let length = (direction[0] * direction[0] + 
                     direction[1] * direction[1] + 
                     direction[2] * direction[2]).sqrt();
        
        if length < EPSILON {
            return None; // Degenerate case: all three points collinear
        }

        let normalized = [
            direction[0] / length,
            direction[1] / length,
            direction[2] / length,
        ];

        // The two circumcenters are at +/- normalized direction
        let p1 = SphericalPoint::from_cartesian(normalized);
        let p2 = SphericalPoint::from_cartesian([
            -normalized[0],
            -normalized[1],
            -normalized[2],
        ]);

        Some((p1, p2))
    }

    /// Check if a point is a valid Voronoi vertex
    /// (closer to incident seeds than to any other seed)
    fn is_voronoi_vertex(point: SphericalPoint, incident_indices: &[usize], all_seeds: &[SphericalPoint]) -> bool {
        if incident_indices.is_empty() {
            return false;
        }

        let reference_dist = distance(point, all_seeds[incident_indices[0]]);

        for (i, seed) in all_seeds.iter().enumerate() {
            if incident_indices.contains(&i) {
                // Should be equidistant
                if (distance(point, *seed) - reference_dist).abs() > EPSILON {
                    return false;
                }
            } else {
                // Should be farther away
                if distance(point, *seed) < reference_dist - EPSILON {
                    return false;
                }
            }
        }

        true
    }

    /// Remove duplicate vertices (within epsilon distance)
    /// CRITICAL: Merges incident_seeds when vertices are duplicates
    fn deduplicate_vertices(vertices: Vec<VoronoiVertex>) -> Vec<VoronoiVertex> {
        let mut unique: Vec<VoronoiVertex> = Vec::new();

        for vertex in vertices {
            let mut found_duplicate = false;
            
            // Check if this vertex is close to any existing unique vertex
            for existing in &mut unique {
                if distance(vertex.position, existing.position) < DEDUPLICATION_EPSILON {
                    // Merge the incident seeds
                    for &seed_id in &vertex.incident_seeds {
                        if !existing.incident_seeds.contains(&seed_id) {
                            existing.incident_seeds.push(seed_id);
                        }
                    }
                    found_duplicate = true;
                    break;
                }
            }
            
            if !found_duplicate {
                unique.push(vertex);
            }
        }

        unique
    }

    /// Compute all Voronoi edges from vertices
    fn compute_edges(seeds: &[SphericalPoint], vertices: &[VoronoiVertex]) -> Vec<VoronoiEdge> {
        let mut edges = Vec::new();

        // For each pair of seeds, check if they share a Voronoi edge
        for i in 0..seeds.len() {
            for j in (i + 1)..seeds.len() {
                // Find vertices that are incident to both i and j
                let shared_vertices: Vec<usize> = vertices
                    .iter()
                    .enumerate()
                    .filter_map(|(idx, v)| {
                        if v.incident_seeds.contains(&i) && v.incident_seeds.contains(&j) {
                            Some(idx)
                        } else {
                            None
                        }
                    })
                    .collect();

                // If exactly 2 vertices, we have an edge
                if shared_vertices.len() == 2 {
                    let vertex1 = vertices[shared_vertices[0]].position;
                    let vertex2 = vertices[shared_vertices[1]].position;
                    
                    // Skip degenerate edges (where vertices are too close)
                    let edge_length = distance(vertex1, vertex2);
                    if edge_length < EPSILON {
                        continue; // Skip this degenerate edge
                    }
                    
                    let geometry = Self::sample_edge_geometry(
                        vertex1,
                        vertex2,
                        seeds[i],
                        seeds[j],
                    );

                    // Verify the geometry is valid (not all points the same)
                    if geometry.len() >= 2 {
                        let geom_length = distance(geometry[0], *geometry.last().unwrap());
                        if geom_length > EPSILON {
                            edges.push(VoronoiEdge {
                                plate_a: i,
                                plate_b: j,
                                vertices: (shared_vertices[0], shared_vertices[1]),
                                geometry,
                            });
                        }
                    }
                }
            }
        }

        edges
    }

    /// Sample points along a Voronoi edge (great circle arc between two vertices)
    fn sample_edge_geometry(
        vertex1: SphericalPoint,
        vertex2: SphericalPoint,
        _seed_a: SphericalPoint,
        _seed_b: SphericalPoint,
    ) -> Vec<SphericalPoint> {
        // Calculate the arc length to determine how many samples we need
        let arc_length = distance(vertex1, vertex2);
        
        // Target: ~75km per segment on Earth
        // This ensures no segment exceeds our evolution thresholds
        let num_samples = ((arc_length / TARGET_SEGMENT_LENGTH_RADIANS).ceil() as usize).max(20);
        
        let mut geometry = Vec::with_capacity(num_samples);

        // Convert to Cartesian
        let cart1 = vertex1.to_cartesian();
        let cart2 = vertex2.to_cartesian();

        // Calculate the arc length
        let dot = cart1[0] * cart2[0] + cart1[1] * cart2[1] + cart1[2] * cart2[2];
        let arc_angle = dot.clamp(-1.0, 1.0).acos();

        // Sample along the great circle arc using slerp (spherical linear interpolation)
        for i in 0..num_samples {
            let t = i as f64 / (num_samples - 1) as f64;
            
            // Slerp formula
            let sin_arc = arc_angle.sin();
            if sin_arc < EPSILON {
                // Points are too close, just use linear interpolation
                geometry.push(vertex1);
                continue;
            }

            let a = ((1.0 - t) * arc_angle).sin() / sin_arc;
            let b = (t * arc_angle).sin() / sin_arc;

            let point_cart = [
                a * cart1[0] + b * cart2[0],
                a * cart1[1] + b * cart2[1],
                a * cart1[2] + b * cart2[2],
            ];

            geometry.push(SphericalPoint::from_cartesian(point_cart));
        }

        geometry
    }

    /// Cross product of two 3D vectors
    fn cross_product(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    }

    /// Get all edges for a specific plate
    pub fn edges_for_plate(&self, plate_id: PlateId) -> Vec<&VoronoiEdge> {
        self.edges
            .iter()
            .filter(|e| e.plate_a == plate_id || e.plate_b == plate_id)
            .collect()
    }

    /// Check if two plates are actually adjacent (share an edge)
    pub fn are_adjacent(&self, plate_a: PlateId, plate_b: PlateId) -> bool {
        self.edges.iter().any(|e| {
            (e.plate_a == plate_a && e.plate_b == plate_b) ||
            (e.plate_a == plate_b && e.plate_b == plate_a)
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tetrahedron_has_correct_topology() {
        // Place 4 seeds at vertices of a regular tetrahedron
        // Each should have 3 edges, 6 total edges, 4 vertices
        let seeds = vec![
            SphericalPoint::from_cartesian([1.0, 1.0, 1.0]),
            SphericalPoint::from_cartesian([1.0, -1.0, -1.0]),
            SphericalPoint::from_cartesian([-1.0, 1.0, -1.0]),
            SphericalPoint::from_cartesian([-1.0, -1.0, 1.0]),
        ];

        let diagram = SphericalVoronoiDiagram::new(seeds);

        // For a tetrahedron: 4 faces, 6 edges, 4 vertices
        // In Voronoi: 4 vertices (dual of faces), 6 edges (dual of edges)
        assert_eq!(diagram.vertices.len(), 4, "Tetrahedron should have 4 Voronoi vertices");
        assert_eq!(diagram.edges.len(), 6, "Tetrahedron should have 6 Voronoi edges");

        // Each seed should be incident to exactly 3 edges
        for i in 0..4 {
            let edge_count = diagram.edges_for_plate(i).len();
            assert_eq!(edge_count, 3, "Each tetrahedron vertex should have 3 edges");
        }
    }

    #[test]
    fn edge_points_are_equidistant() {
        let seeds = vec![
            SphericalPoint::new(0.0, 0.0),
            SphericalPoint::new(0.5, 0.0),
            SphericalPoint::new(0.0, 0.5),
            SphericalPoint::new(0.5, 0.5),
        ];

        let diagram = SphericalVoronoiDiagram::new(seeds);

        for edge in &diagram.edges {
            for point in &edge.geometry {
                let dist_a = distance(*point, diagram.seeds[edge.plate_a]);
                let dist_b = distance(*point, diagram.seeds[edge.plate_b]);

                // Points on the edge should be equidistant to both seeds
                assert!(
                    (dist_a - dist_b).abs() < 0.01,
                    "Edge point not equidistant: dist_a={:.6}, dist_b={:.6}, diff={:.6}",
                    dist_a, dist_b, (dist_a - dist_b).abs()
                );

                // And farther from all other seeds
                for (i, seed) in diagram.seeds.iter().enumerate() {
                    if i != edge.plate_a && i != edge.plate_b {
                        let dist_other = distance(*point, *seed);
                        assert!(
                            dist_other >= dist_a - 0.01,
                            "Edge point closer to non-incident seed"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn vertices_are_equidistant_to_incident_seeds() {
        let seeds = vec![
            SphericalPoint::new(0.0, 0.0),
            SphericalPoint::new(1.0, 0.0),
            SphericalPoint::new(0.0, 1.0),
            SphericalPoint::new(1.0, 1.0),
            SphericalPoint::new(0.5, 0.5),
        ];

        let diagram = SphericalVoronoiDiagram::new(seeds);

        for vertex in &diagram.vertices {
            assert!(
                vertex.incident_seeds.len() >= 3,
                "Voronoi vertex should be incident to at least 3 seeds"
            );

            let first_seed = diagram.seeds[vertex.incident_seeds[0]];
            let reference_dist = distance(vertex.position, first_seed);

            for &seed_idx in &vertex.incident_seeds {
                let dist = distance(vertex.position, diagram.seeds[seed_idx]);
                assert!(
                    (dist - reference_dist).abs() < 0.001,
                    "Vertex not equidistant to incident seeds: {:.6} vs {:.6}",
                    dist, reference_dist
                );
            }
        }
    }

    #[test]
    fn every_edge_has_two_endpoints() {
        let seeds = vec![
            SphericalPoint::new(0.0, 0.0),
            SphericalPoint::new(0.8, 0.0),
            SphericalPoint::new(0.0, 0.8),
            SphericalPoint::new(0.8, 0.8),
            SphericalPoint::new(0.4, 0.4),
            SphericalPoint::new(-0.4, 0.4),
        ];

        let diagram = SphericalVoronoiDiagram::new(seeds);

        for edge in &diagram.edges {
            let v1_idx = edge.vertices.0;
            let v2_idx = edge.vertices.1;

            assert!(v1_idx < diagram.vertices.len());
            assert!(v2_idx < diagram.vertices.len());
            assert_ne!(v1_idx, v2_idx);

            // Both vertices should be incident to both plates
            let v1 = &diagram.vertices[v1_idx];
            let v2 = &diagram.vertices[v2_idx];

            assert!(v1.incident_seeds.contains(&edge.plate_a));
            assert!(v1.incident_seeds.contains(&edge.plate_b));
            assert!(v2.incident_seeds.contains(&edge.plate_a));
            assert!(v2.incident_seeds.contains(&edge.plate_b));
        }
    }

    #[test]
    fn adjacency_detection_is_symmetric() {
        let seeds = vec![
            SphericalPoint::new(0.0, 0.0),
            SphericalPoint::new(0.5, 0.0),
            SphericalPoint::new(0.0, 0.5),
            SphericalPoint::new(0.5, 0.5),
        ];

        let diagram = SphericalVoronoiDiagram::new(seeds);

        for i in 0..diagram.seeds.len() {
            for j in 0..diagram.seeds.len() {
                assert_eq!(
                    diagram.are_adjacent(i, j),
                    diagram.are_adjacent(j, i),
                    "Adjacency should be symmetric"
                );
            }
        }
    }

    #[test]
    fn edge_count_follows_eulers_formula() {
        // For a spherical polyhedron: V - E + F = 2 (Euler's formula)
        // In Voronoi: vertices = dual faces, edges = dual edges, faces = seeds
        // So: voronoi_vertices - voronoi_edges + num_seeds = 2
        let seeds = vec![
            SphericalPoint::new(0.0, 0.0),
            SphericalPoint::new(1.0, 0.0),
            SphericalPoint::new(0.5, 1.0),
            SphericalPoint::new(-0.5, 1.0),
            SphericalPoint::new(0.0, -1.0),
            SphericalPoint::new(1.0, 1.0),
        ];

        let diagram = SphericalVoronoiDiagram::new(seeds);

        let v = diagram.vertices.len() as i32;
        let e = diagram.edges.len() as i32;
        let f = diagram.seeds.len() as i32;

        let euler = v - e + f;
        assert_eq!(
            euler, 2,
            "Euler's formula violated: V={}, E={}, F={}, V-E+F={}",
            v, e, f, euler
        );
    }

    #[test]
    fn all_seeds_have_at_least_one_edge() {
        let seeds = vec![
            SphericalPoint::new(0.0, 0.0),
            SphericalPoint::new(0.5, 0.0),
            SphericalPoint::new(0.0, 0.5),
            SphericalPoint::new(0.5, 0.5),
            SphericalPoint::new(0.25, 0.25),
        ];

        let diagram = SphericalVoronoiDiagram::new(seeds);

        for i in 0..diagram.seeds.len() {
            let edges = diagram.edges_for_plate(i);
            assert!(
                !edges.is_empty(),
                "Seed {} has no edges",
                i
            );
        }
    }
}

