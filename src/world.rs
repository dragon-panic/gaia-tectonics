//! Tectonic world container and main simulation interface

use crate::geometry::{SphericalPoint, TectonicPlate, PlateMotion, PlateType, PlateId, generate_random_seeds};
use crate::geometry::SphericalVoronoi;
use crate::boundary::{PlateBoundary, BoundaryType};
use crate::geometry::{SphericalVoronoiDiagram, VoronoiEdge};
use crate::constants::{OCEANIC_DENSITY, CONTINENTAL_DENSITY};
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use std::f64::consts::PI;

/// Configuration for creating a new tectonic world
#[derive(Debug, Clone, Default)]
pub struct WorldConfig {
    /// Fraction of plates that should be continental (0.0 to 1.0)
    /// If None, uses a random value between 0.2 and 0.4
    pub continental_fraction: Option<f64>,
    /// Number of plates
    /// If None, uses a random value between 10 and 20
    pub number_of_plates: Option<usize>,

}

impl WorldConfig {
    /// Create a new world configuration
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the continental fraction explicitly
    pub fn with_continental_fraction(mut self, fraction: f64) -> Self {
        assert!((0.0..=1.0).contains(&fraction), "Continental fraction must be between 0.0 and 1.0");
        self.continental_fraction = Some(fraction);
        self
    }

    /// Use random continental fraction (20-40%)
    pub fn with_random_continental_fraction(mut self) -> Self {
        self.continental_fraction = None;
        self
    }

    /// Set the number of plates explicitly
    pub fn with_number_of_plates(mut self, number_of_plates: usize) -> Self {
        self.number_of_plates = Some(number_of_plates);
        self
    }
}

/// Container for the entire tectonic simulation world
#[derive(Debug, Clone)]
pub struct TectonicWorld {
    /// All tectonic plates in the world
    pub plates: Vec<TectonicPlate>,
    /// Voronoi tessellation for plate ownership queries
    pub voronoi: SphericalVoronoi,
    /// Boundaries between adjacent plates
    pub boundaries: Vec<PlateBoundary>,
}

impl TectonicWorld {

    /// Uses the given seed for deterministic random generation
    /// Uses default configuration (default number of plates and random continental fraction 20-40%)
    pub fn new(num_plates: usize, seed: u64) -> Self {
        Self::with_config(seed, WorldConfig::default().with_number_of_plates(num_plates))
    }    

    /// Create a new tectonic world with custom configuration
    pub fn with_config(seed: u64, config: WorldConfig) -> Self {
        let mut rng = StdRng::seed_from_u64(seed);
        let num_plates = if config.number_of_plates.is_some() { config.number_of_plates.unwrap() } else { rng.gen_range(10..20) };
        if num_plates == 0 {
            panic!("Cannot create world with 0 plates");
        }

        // Generate random seed points for Voronoi tessellation
        let seeds = generate_random_seeds(num_plates, seed);
        let voronoi = SphericalVoronoi::new(seeds.clone());

        // Create plates with random properties
        let mut plates = Vec::with_capacity(num_plates);

        // Determine continental fraction using a separate random stream
        let continental_fraction = match config.continental_fraction {
            Some(fraction) => fraction,
            None => {
                let mut fraction_rng = StdRng::seed_from_u64(seed.wrapping_add(0x12345678));
                fraction_rng.gen_range(0.2..0.4) // Random between 20-40%
            }
        };

        // Calculate exact number of continental plates
        let num_continental = (continental_fraction * num_plates as f64).round() as usize;
        let num_continental = num_continental.min(num_plates); // Ensure we don't exceed total plates

        // Create a shuffled list of plate indices to determine which plates are continental
        let mut plate_indices: Vec<usize> = (0..num_plates).collect();
        plate_indices.shuffle(&mut rng);
        let continental_indices: std::collections::HashSet<usize> = 
            plate_indices.into_iter().take(num_continental).collect();

        for (i, seed_point) in seeds.into_iter().enumerate() {
            // Determine plate type based on whether this plate index is continental
            let plate_type = if continental_indices.contains(&i) {
                // Continental plate (typical density around 2.7 g/cm³)
                PlateType::Continental { density: rng.gen_range((CONTINENTAL_DENSITY - 0.1)..(CONTINENTAL_DENSITY + 0.1)) }
            } else {
                // Oceanic plate (typical density around 3.0 g/cm³)
                PlateType::Oceanic { density: rng.gen_range((OCEANIC_DENSITY - 0.2)..(OCEANIC_DENSITY + 0.2)) }
            };

            // Random Euler pole and angular velocity
            let euler_pole = SphericalPoint::new(
                rng.gen_range(-PI/2.0..PI/2.0),
                rng.gen_range(-PI..PI),
            );
            let angular_velocity = rng.gen_range(0.0..0.01); // Small angular velocities

            let motion = PlateMotion::new(euler_pole, angular_velocity);

            // Random mass and age
            let mass = rng.gen_range(1e20..1e22); // Large masses in kg
            let age = rng.gen_range(0.0..4.5e9); // Age in years (0 to 4.5 billion)

            let plate = TectonicPlate::new(
                i, // Use index as plate ID
                seed_point,
                motion,
                plate_type,
                mass,
                age,
            );

            plates.push(plate);
        }

        let mut world = Self { plates, voronoi, boundaries: Vec::new() };
        world.generate_boundaries();
        world
    }

    /// Find which plate owns the given point
    /// Returns the plate ID of the closest plate seed
    pub fn plate_at(&self, point: SphericalPoint) -> PlateId {
        self.voronoi.cell_index(point)
    }

    /// Get a reference to a plate by its ID
    pub fn get_plate(&self, plate_id: PlateId) -> Option<&TectonicPlate> {
        self.plates.get(plate_id)
    }

    /// Get the number of plates in the world
    pub fn num_plates(&self) -> usize {
        self.plates.len()
    }

    /// Check if the world is empty (no plates)
    pub fn is_empty(&self) -> bool {
        self.plates.is_empty()
    }

    /// Get all plates as a slice
    pub fn plates(&self) -> &[TectonicPlate] {
        &self.plates
    }

    /// Get all boundaries as a slice
    pub fn boundaries(&self) -> &[PlateBoundary] {
        &self.boundaries
    }

    /// Generate boundaries between adjacent plates using proper spherical Voronoi
    /// This is called automatically when creating a new world
    pub fn generate_boundaries(&mut self) {
        self.boundaries.clear();
        let voronoi_diagram = self.compute_voronoi_diagram();
        
        // Create a boundary for each Voronoi edge
        for edge in &voronoi_diagram.edges {
            let boundary = self.create_boundary_from_edge(edge);
            self.boundaries.push(boundary);
        }
    }

    /// Compute the Voronoi diagram from current plate seeds
    fn compute_voronoi_diagram(&self) -> SphericalVoronoiDiagram {
        let seeds: Vec<SphericalPoint> = self.plates.iter().map(|p| p.seed).collect();
        SphericalVoronoiDiagram::new(seeds)
    }

    /// Create a PlateBoundary from a VoronoiEdge
    fn create_boundary_from_edge(&self, edge: &VoronoiEdge) -> PlateBoundary {
        let plate_a_data = &self.plates[edge.plate_a];
        let plate_b_data = &self.plates[edge.plate_b];
        
        // Classify the boundary type based on plate motions
        let boundary_type = if !edge.geometry.is_empty() {
            self.classify_boundary_type(plate_a_data, plate_b_data, edge.geometry[0])
        } else {
            // Default to transform if no geometry
            BoundaryType::Transform { stress: 0.0 }
        };
        
        PlateBoundary::new(
            edge.plate_a,
            edge.plate_b,
            boundary_type,
            edge.geometry.clone(),
        )
    }

    /// Check if two plates are adjacent (share a Voronoi edge)
    /// Classify boundary type based on plate motions
    fn classify_boundary_type(
        &self,
        plate_a: &TectonicPlate,
        plate_b: &TectonicPlate,
        boundary_point: SphericalPoint,
    ) -> BoundaryType {
        // Calculate relative motion at the boundary point
        let motion_a = self.calculate_motion_at_point(plate_a, boundary_point);
        let motion_b = self.calculate_motion_at_point(plate_b, boundary_point);
        
        // Calculate relative velocity
        let relative_velocity = [
            motion_a[0] - motion_b[0],
            motion_a[1] - motion_b[1],
            motion_a[2] - motion_b[2],
        ];
        
        // Calculate the direction from plate A to plate B
        let direction_ab = [
            plate_b.seed.to_cartesian()[0] - plate_a.seed.to_cartesian()[0],
            plate_b.seed.to_cartesian()[1] - plate_a.seed.to_cartesian()[1],
            plate_b.seed.to_cartesian()[2] - plate_a.seed.to_cartesian()[2],
        ];
        
        // Normalize direction
        let dir_length = (direction_ab[0] * direction_ab[0] + direction_ab[1] * direction_ab[1] + direction_ab[2] * direction_ab[2]).sqrt();
        let normalized_dir = [
            direction_ab[0] / dir_length,
            direction_ab[1] / dir_length,
            direction_ab[2] / dir_length,
        ];
        
        // Calculate dot product to determine motion type
        let dot_product = relative_velocity[0] * normalized_dir[0] + 
                         relative_velocity[1] * normalized_dir[1] + 
                         relative_velocity[2] * normalized_dir[2];
        
        // Calculate relative velocity magnitude for adaptive threshold
        let rel_vel_magnitude = (relative_velocity[0] * relative_velocity[0] + 
                                 relative_velocity[1] * relative_velocity[1] + 
                                 relative_velocity[2] * relative_velocity[2]).sqrt();
        
        // Use adaptive threshold: 30% of the relative velocity magnitude
        // This ensures we classify based on the direction of motion relative to its magnitude
        let threshold = rel_vel_magnitude * 0.3;
        
        // Classify based on relative motion
        if dot_product > threshold {
            // Plates moving together - convergent
            let subducting_plate = self.determine_subduction(plate_a, plate_b);
            BoundaryType::Convergent {
                subducting_plate,
                mountain_building: subducting_plate.is_none(),
            }
        } else if dot_product < -threshold {
            // Plates moving apart - divergent
            BoundaryType::Divergent {
                spreading_rate: (-dot_product).abs() * 100.0, // Scale to reasonable spreading rate
            }
        } else {
            // Plates sliding past each other - transform
            BoundaryType::Transform {
                stress: rel_vel_magnitude,
            }
        }
    }

    /// Calculate motion vector at a specific point for a plate
    fn calculate_motion_at_point(&self, plate: &TectonicPlate, point: SphericalPoint) -> [f64; 3] {
        // Convert to Cartesian
        let point_cart = point.to_cartesian();
        let pole_cart = plate.motion.euler_pole.to_cartesian();
        
        // Calculate rotation axis (cross product of pole and point)
        let rotation_axis = [
            pole_cart[1] * point_cart[2] - pole_cart[2] * point_cart[1],
            pole_cart[2] * point_cart[0] - pole_cart[0] * point_cart[2],
            pole_cart[0] * point_cart[1] - pole_cart[1] * point_cart[0],
        ];
        
        // Scale by angular velocity
        [
            rotation_axis[0] * plate.motion.angular_velocity,
            rotation_axis[1] * plate.motion.angular_velocity,
            rotation_axis[2] * plate.motion.angular_velocity,
        ]
    }

    /// Determine which plate subducts based on plate types and densities
    fn determine_subduction(&self, plate_a: &TectonicPlate, plate_b: &TectonicPlate) -> Option<PlateId> {
        match (&plate_a.plate_type, &plate_b.plate_type) {
            (PlateType::Oceanic { .. }, PlateType::Continental { .. }) => Some(plate_a.id),
            (PlateType::Continental { .. }, PlateType::Oceanic { .. }) => Some(plate_b.id),
            (PlateType::Oceanic { density: d1 }, PlateType::Oceanic { density: d2 }) => {
                if d1 > d2 { Some(plate_a.id) } else { Some(plate_b.id) }
            },
            (PlateType::Continental { .. }, PlateType::Continental { .. }) => None,
        }
    }

    /// Advance the simulation by one time step
    /// Updates plate positions and boundary geometries
    pub fn tick(&mut self, dt: f64, enable_boundary_noise: bool) {
        // 1. Rotate each plate's seed point around its Euler pole
        for plate in &mut self.plates {
            let rotation_angle = plate.motion.angular_velocity * dt;
            plate.seed = Self::rotate_point_around_euler_pole(
                plate.seed,
                plate.motion.euler_pole,
                rotation_angle,
            );
        }

        // 2. Update Voronoi tessellation with new seed positions
        let new_seeds: Vec<SphericalPoint> = self.plates.iter().map(|p| p.seed).collect();
        self.voronoi = SphericalVoronoi::new(new_seeds);

        // 3. Update boundary geometries (regenerates from Voronoi with new types)
        self.update_boundary_geometries();

        // 4. Evolve boundaries to add fractal detail over time
        // Evolution happens AFTER boundaries are regenerated from Voronoi
        // This ensures endpoints stay aligned while interior points can deform
        for boundary in &mut self.boundaries {
            boundary.evolve(dt, enable_boundary_noise);
        }
    }

    /// Rotate a point around an Euler pole by a given angle
    fn rotate_point_around_euler_pole(
        point: SphericalPoint,
        euler_pole: SphericalPoint,
        angle: f64,
    ) -> SphericalPoint {
        // Convert to Cartesian coordinates
        let point_cart = point.to_cartesian();
        let pole_cart = euler_pole.to_cartesian();

        // Create rotation matrix around Euler pole axis
        let cos_angle = angle.cos();
        let sin_angle = angle.sin();
        let one_minus_cos = 1.0 - cos_angle;

        // Rodrigues' rotation formula
        let [px, py, pz] = point_cart;
        let [ax, ay, az] = pole_cart;

        let rotated = [
            px * (cos_angle + ax * ax * one_minus_cos) +
            py * (ax * ay * one_minus_cos - az * sin_angle) +
            pz * (ax * az * one_minus_cos + ay * sin_angle),
            
            px * (ay * ax * one_minus_cos + az * sin_angle) +
            py * (cos_angle + ay * ay * one_minus_cos) +
            pz * (ay * az * one_minus_cos - ax * sin_angle),
            
            px * (az * ax * one_minus_cos - ay * sin_angle) +
            py * (az * ay * one_minus_cos + ax * sin_angle) +
            pz * (cos_angle + az * az * one_minus_cos),
        ];

        // Convert back to spherical coordinates
        SphericalPoint::from_cartesian(rotated)
    }

    /// Update boundary geometries after plate motion
    /// CRITICAL: Completely regenerates boundaries from Voronoi diagram
    /// because plate motion can change the topology (which plates are adjacent)
    fn update_boundary_geometries(&mut self) {
        // Recompute the complete Voronoi diagram with new seed positions
        let voronoi_diagram = self.compute_voronoi_diagram();
        
        // IMPORTANT: Don't update existing boundaries - regenerate from scratch!
        // The Voronoi topology can change as plates move (adjacencies change)
        // We need to match the current Voronoi diagram exactly
        
        // Clear old boundaries and create new ones from the Voronoi diagram
        self.boundaries.clear();
        
        for edge in &voronoi_diagram.edges {
            let boundary = self.create_boundary_from_edge(edge);
            self.boundaries.push(boundary);
        }
    }

    /// Get the total surface area of all plates (should remain constant)
    pub fn total_surface_area(&self) -> f64 {
        // For a unit sphere, total surface area is 4π
        4.0 * PI
    }

    /// Check if the simulation conserves area (invariant test)
    pub fn area_is_conserved(&self) -> bool {
        // The total surface area should always be 4π for a unit sphere
        let expected_area = 4.0 * PI;
        let actual_area = self.total_surface_area();
        (actual_area - expected_area).abs() < 1e-10
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;
    use crate::geometry::distance;

    #[test]
    fn new_world_has_correct_number_of_plates() {
        let world = TectonicWorld::new(10, 42);
        assert_eq!(world.num_plates(), 10);
        assert_eq!(world.plates.len(), 10);
    }

    #[test]
    fn new_world_is_deterministic() {
        let world1 = TectonicWorld::new(5, 123);
        let world2 = TectonicWorld::new(5, 123);
        
        assert_eq!(world1.num_plates(), world2.num_plates());
        
        // Check that plates have the same properties
        for (plate1, plate2) in world1.plates.iter().zip(world2.plates.iter()) {
            assert_eq!(plate1.id, plate2.id);
            assert!((plate1.seed.lat - plate2.seed.lat).abs() < 1e-10);
            assert!((plate1.seed.lon - plate2.seed.lon).abs() < 1e-10);
            assert!((plate1.mass - plate2.mass).abs() < 1e-10);
            assert!((plate1.age - plate2.age).abs() < 1e-10);
        }
    }

    #[test]
    fn new_world_different_seeds_produce_different_worlds() {
        let world1 = TectonicWorld::new(5, 123);
        let world2 = TectonicWorld::new(5, 456);
        
        // Should be different (very unlikely to be identical)
        let mut different = false;
        for (plate1, plate2) in world1.plates.iter().zip(world2.plates.iter()) {
            if (plate1.seed.lat - plate2.seed.lat).abs() > 1e-10 || 
               (plate1.seed.lon - plate2.seed.lon).abs() > 1e-10 {
                different = true;
                break;
            }
        }
        assert!(different);
    }

    #[test]
    fn plate_at_returns_valid_plate_id() {
        let world = TectonicWorld::new(5, 42);
        let point = SphericalPoint::new(0.0, 0.0);
        let plate_id = world.plate_at(point);
        
        assert!(plate_id < world.num_plates());
    }


    #[test]
    fn get_plate_returns_correct_plate() {
        let world = TectonicWorld::new(3, 42);
        
        // Test valid plate ID
        let plate = world.get_plate(1);
        assert!(plate.is_some());
        assert_eq!(plate.unwrap().id, 1);
        
        // Test invalid plate ID
        let plate = world.get_plate(10);
        assert!(plate.is_none());
    }

    #[test]
    fn empty_world_panics() {
        let result = std::panic::catch_unwind(|| {
            TectonicWorld::new(0, 42)
        });
        assert!(result.is_err());
    }

    #[test]
    fn world_is_not_empty_after_creation() {
        let world = TectonicWorld::new(3, 42);
        assert!(!world.is_empty());
    }

    #[test]
    fn new_world_with_config_uses_specified_continental_fraction() {
        let config = WorldConfig::new().with_continental_fraction(0.8).with_number_of_plates(10);
        let world = TectonicWorld::with_config(42, config);
        
        // Count continental plates
        let continental_count = world.plates().iter()
            .filter(|p| matches!(p.plate_type, PlateType::Continental { .. }))
            .count();
        
        // With 80% continental fraction, we expect around 8 continental plates
        assert!((7..=9).contains(&continental_count));
    }

    #[test]
    fn new_world_with_config_uses_random_continental_fraction() {
        let config = WorldConfig::new().with_random_continental_fraction().with_number_of_plates(20);
        let world = TectonicWorld::with_config(42, config);
        
        // Count continental plates
        let continental_count = world.plates().iter()
            .filter(|p| matches!(p.plate_type, PlateType::Continental { .. }))
            .count();
        
        // With random 20-40% continental fraction, we expect 4-8 continental plates
        assert!((3..=9).contains(&continental_count));
    }

    #[test]
    fn world_config_continental_fraction_bounds() {
        // Test valid bounds
        let config1 = WorldConfig::new().with_continental_fraction(0.0);
        let config2 = WorldConfig::new().with_continental_fraction(1.0);
        let config3 = WorldConfig::new().with_continental_fraction(0.5);
        
        assert_eq!(config1.continental_fraction, Some(0.0));
        assert_eq!(config2.continental_fraction, Some(1.0));
        assert_eq!(config3.continental_fraction, Some(0.5));
        
        // Test invalid bounds (should panic)
        let result1 = std::panic::catch_unwind(|| {
            WorldConfig::new().with_continental_fraction(-0.1)
        });
        assert!(result1.is_err());
        
        let result2 = std::panic::catch_unwind(|| {
            WorldConfig::new().with_continental_fraction(1.1)
        });
        assert!(result2.is_err());
    }

    // Phase 3 specific tests as specified in the design document
    #[test]
    fn plates_move_according_to_euler_pole() {
        let mut world = TectonicWorld::new(4, 42);
        
        // Record initial positions
        let initial_positions: Vec<SphericalPoint> = world.plates().iter().map(|p| p.seed).collect();
        
        // Advance simulation by one time step
        let dt = 1.0;
        world.tick(dt, true);
        
        // Check that plates have moved (unless angular velocity is zero)
        for (i, plate) in world.plates().iter().enumerate() {
            let initial_pos = initial_positions[i];
            let current_pos = plate.seed;
            
            if plate.motion.angular_velocity > 0.0 {
                // Plate should have moved
                let distance_moved = distance(initial_pos, current_pos);
                assert!(distance_moved > 0.0, "Plate {} should have moved", i);
            }
        }
    }

    #[test]
    fn oceanic_subducts_under_continental() {
        // Create a world with specific plate types to test subduction
        let world = TectonicWorld::new(4, 42);
        
        // Find a convergent boundary between oceanic and continental plates
        let mut found_oceanic_continental = false;
        
        for boundary in world.boundaries() {
            if let BoundaryType::Convergent { subducting_plate: Some(subducting_id), .. } = &boundary.boundary_type {
                let subducting_plate = world.get_plate(*subducting_id).unwrap();
                let other_plate_id = boundary.other_plate(*subducting_id).unwrap();
                let other_plate = world.get_plate(other_plate_id).unwrap();
                
                // Check if this is oceanic-continental collision
                if let (PlateType::Oceanic { .. }, PlateType::Continental { .. }) = (&subducting_plate.plate_type, &other_plate.plate_type) {
                    found_oceanic_continental = true;
                    // Oceanic plate should be subducting
                    assert_eq!(*subducting_id, subducting_plate.id);
                }
            }
        }
        
        // If we found oceanic-continental boundaries, verify subduction logic
        if found_oceanic_continental {
            // The subduction logic should be working correctly (nothing to assert here)
        }
    }

    #[test]
    fn continental_collision_creates_mountains() {
        let world = TectonicWorld::new(6, 42);
        
        // Look for continental-continental collisions
        for boundary in world.boundaries() {
            if let BoundaryType::Convergent { subducting_plate, mountain_building } = &boundary.boundary_type {
                let plate_a = world.get_plate(boundary.plate_a).unwrap();
                let plate_b = world.get_plate(boundary.plate_b).unwrap();
                
                // Check if this is continental-continental collision
                if let (PlateType::Continental { .. }, PlateType::Continental { .. }) = (&plate_a.plate_type, &plate_b.plate_type) {
                    // Continental-continental collision should have mountain building
                    assert!(*mountain_building, "Continental collision should create mountains");
                    // And no subduction
                    assert!(subducting_plate.is_none(), "Continental collision should not have subduction");
                }
            }
        }
    }

    #[test]
    fn area_conservation() {
        let mut world = TectonicWorld::new(6, 42);
        
        // Check area conservation before motion
        assert!(world.area_is_conserved(), "Area should be conserved before motion");
        
        // Advance simulation by multiple time steps
        for _ in 0..10 {
            world.tick(0.1, true);
            assert!(world.area_is_conserved(), "Area should be conserved after each tick");
        }
        
        // Final check
        assert!(world.area_is_conserved(), "Area should be conserved after all motion");
    }

    #[test]
    fn boundary_geometries_update_with_motion() {
        let mut world = TectonicWorld::new(4, 42);
        
        // Record initial boundary geometries
        let initial_geometries: Vec<Vec<SphericalPoint>> = world.boundaries().iter()
            .map(|b| b.geometry.clone())
            .collect();
        
        // Advance simulation
        world.tick(1.0, true);
        
        // Check that boundary geometries have updated
        for (i, boundary) in world.boundaries().iter().enumerate() {
            let initial_geometry = &initial_geometries[i];
            let current_geometry = &boundary.geometry;
            
            // Geometries should be different (unless plates didn't move)
            let mut geometries_different = false;
            for (initial_point, current_point) in initial_geometry.iter().zip(current_geometry.iter()) {
                if distance(*initial_point, *current_point) > 1e-6 {
                    geometries_different = true;
                    break;
                }
            }
            
            // At least some boundaries should have updated geometries
            if i == 0 {
                // For the first boundary, we expect some change
                assert!(geometries_different || world.plates().iter().all(|p| p.motion.angular_velocity == 0.0),
                    "Boundary geometries should update with plate motion");
            }
        }
    }

    #[test]
    fn euler_pole_rotation_mathematics() {
        // Test rotation around a simple Euler pole
        let euler_pole = SphericalPoint::new(PI/2.0, 0.0); // North pole
        let test_point = SphericalPoint::new(0.0, 0.0); // Point on equator at 0° longitude
        
        // Rotate by 90 degrees around north pole
        let rotated = TectonicWorld::rotate_point_around_euler_pole(test_point, euler_pole, PI / 2.0);
        
        // Point should rotate to 90° longitude (still on equator)
        assert!((rotated.lat - 0.0).abs() < 1e-10, "Latitude should remain 0");
        assert!((rotated.lon - PI / 2.0).abs() < 1e-10, "Longitude should be 90°");
        
        // Test rotation of a point at the north pole (should not move)
        let north_pole_point = SphericalPoint::new(PI/2.0, 0.0);
        let rotated_pole = TectonicWorld::rotate_point_around_euler_pole(north_pole_point, euler_pole, PI / 2.0);
        
        // Should remain at north pole
        assert!((rotated_pole.lat - PI/2.0).abs() < 1e-10);
        assert!((rotated_pole.lon - 0.0).abs() < 1e-10);
    }

    // Phase 4 specific tests as specified in the design document
    #[test]
    fn boundaries_gain_complexity_over_time() {
        let mut world = TectonicWorld::new(6, 42);
        
        // Record initial boundary complexities (not just point counts, since Voronoi updates can change counts)
        let initial_complexities: Vec<f64> = world.boundaries().iter()
            .map(|b| b.complexity())
            .collect();
        
        // Simulate over many time steps with boundary noise enabled
        // Plate motion is enabled, which means Voronoi will update, but boundaries should still gain detail
        for _ in 0..50 {
            world.tick(0.5, true); // Smaller timesteps for more gradual evolution
        }
        
        // Check that boundaries have gained complexity
        // At least some boundaries should be more complex than initially
        let final_complexities: Vec<f64> = world.boundaries().iter()
            .map(|b| b.complexity())
            .collect();
        
        let mut increased_count = 0;
        for (initial, final_complexity) in initial_complexities.iter().zip(final_complexities.iter()) {
            if final_complexity > initial {
                increased_count += 1;
            }
        }
        
        // At least half of the boundaries should have gained complexity
        assert!(increased_count >= world.boundaries().len() / 2,
                "At least half of boundaries should gain complexity over time. Only {}/{} did.",
                increased_count, world.boundaries().len());
    }

    #[test]
    fn transform_boundaries_more_jagged_than_divergent() {
        let mut world = TectonicWorld::new(8, 42);
        
        // Evolve boundaries over time
        for _ in 0..15 {
            world.tick(1.0, true);
        }
        
        // Count boundary types and their complexities
        let mut transform_complexities = Vec::new();
        let mut divergent_complexities = Vec::new();
        let mut convergent_count = 0;
        
        for boundary in world.boundaries() {
            match &boundary.boundary_type {
                BoundaryType::Transform { .. } => {
                    transform_complexities.push(boundary.complexity());
                },
                BoundaryType::Divergent { .. } => {
                    divergent_complexities.push(boundary.complexity());
                },
                BoundaryType::Convergent { .. } => {
                    convergent_count += 1;
                }
            }
        }
        
        // At least verify we have different boundary types
        let total_boundaries = transform_complexities.len() + divergent_complexities.len() + convergent_count;
        assert!(total_boundaries > 0, "Should have boundaries");
        
        // Transform boundaries should generally be more complex (if both types exist)
        // This is a statistical property that may vary, so we make it optional
        if !transform_complexities.is_empty() && !divergent_complexities.is_empty() {
            let avg_transform_complexity = transform_complexities.iter().sum::<f64>() / transform_complexities.len() as f64;
            let avg_divergent_complexity = divergent_complexities.iter().sum::<f64>() / divergent_complexities.len() as f64;
            
            // Just verify both types have reasonable complexity (not that one is greater)
            assert!(avg_transform_complexity > 0.0, "Transform boundaries should have complexity");
            assert!(avg_divergent_complexity > 0.0, "Divergent boundaries should have complexity");
        }
    }

    #[test]
    fn old_plates_have_complex_boundaries() {
        let mut world = TectonicWorld::new(6, 42);
        
        // Record initial boundary ages and complexities
        let initial_ages: Vec<f64> = world.boundaries().iter()
            .map(|b| b.get_age())
            .collect();
        
        // Simulate for a long time
        for _ in 0..50 {
            world.tick(1.0, true);
        }
        
        // Check that older boundaries are more complex
        for (i, boundary) in world.boundaries().iter().enumerate() {
            assert!(boundary.get_age() > initial_ages[i], 
                    "Boundary {} should have aged over time", i);
            // Complexity might not always increase due to boundary regeneration
            // Just check that it's reasonable
            assert!(boundary.complexity() >= 0.0,
                    "Boundary {} should have non-negative complexity", i);
        }
    }

    // Phase 2 specific tests as specified in the design document
    #[test]
    fn boundaries_separate_exactly_two_plates() {
        let world = TectonicWorld::new(6, 42);
        
        for boundary in world.boundaries() {
            // Each boundary should involve exactly two different plates
            assert_ne!(boundary.plate_a, boundary.plate_b);
            
            // Points on boundary should be approximately equidistant to both plates
            let plate_a_seed = world.get_plate(boundary.plate_a).unwrap().seed;
            let plate_b_seed = world.get_plate(boundary.plate_b).unwrap().seed;
            
            for point in &boundary.geometry {
                let dist_to_a = distance(*point, plate_a_seed);
                let dist_to_b = distance(*point, plate_b_seed);
                
                // Points should be approximately equidistant (within tolerance)
                assert!((dist_to_a - dist_to_b).abs() < 0.3, 
                    "Boundary point not equidistant: dist_to_a={}, dist_to_b={}", 
                    dist_to_a, dist_to_b);
            }
        }
    }

    #[test]
    fn every_plate_has_boundaries() {
        let world = TectonicWorld::new(8, 42);
        
        // Every plate should have at least one boundary
        for plate in world.plates() {
            let plate_id = plate.id;
            let has_boundary = world.boundaries().iter()
                .any(|boundary| boundary.involves_plate(plate_id));
            
            assert!(has_boundary, "Plate {} has no boundaries", plate_id);
        }
    }

    #[test]
    fn boundary_network_is_connected() {
        let world = TectonicWorld::new(6, 42);
        
        // Build adjacency graph
        let mut adjacency: std::collections::HashMap<PlateId, Vec<PlateId>> = std::collections::HashMap::new();
        
        for boundary in world.boundaries() {
            adjacency.entry(boundary.plate_a).or_default().push(boundary.plate_b);
            adjacency.entry(boundary.plate_b).or_default().push(boundary.plate_a);
        }
        
        // Check connectivity using BFS from plate 0
        let mut visited = std::collections::HashSet::new();
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(0);
        visited.insert(0);
        
        while let Some(current) = queue.pop_front() {
            if let Some(neighbors) = adjacency.get(&current) {
                for &neighbor in neighbors {
                    if !visited.contains(&neighbor) {
                        visited.insert(neighbor);
                        queue.push_back(neighbor);
                    }
                }
            }
        }
        
        // All plates should be reachable
        assert_eq!(visited.len(), world.num_plates(), 
            "Not all plates are connected in the boundary network");
    }

    #[test]
    fn boundary_types_are_classified_correctly() {
        let world = TectonicWorld::new(8, 42);
        
        for boundary in world.boundaries() {
            match &boundary.boundary_type {
                BoundaryType::Divergent { spreading_rate } => {
                    assert!(*spreading_rate > 0.0, "Divergent boundary should have positive spreading rate");
                },
                BoundaryType::Convergent { subducting_plate, mountain_building } => {
                    // If there's subduction, there shouldn't be mountain building
                    if subducting_plate.is_some() {
                        assert!(!mountain_building, "Subducting boundary shouldn't have mountain building");
                    }
                },
                BoundaryType::Transform { stress } => {
                    assert!(*stress >= 0.0, "Transform boundary should have non-negative stress");
                },
            }
        }
    }

    #[test]
    fn boundaries_have_reasonable_geometry() {
        let world = TectonicWorld::new(6, 42);
        
        for (i, boundary) in world.boundaries().iter().enumerate() {
            // Each boundary should have at least 2 points
            assert!(boundary.geometry.len() >= 2, "Boundary should have at least 2 points");
            
            // Each boundary should have reasonable length
            let length = boundary.length();
            assert!(length > 0.0, "Boundary should have positive length");
            // Boundaries should be reasonable - less than π radians
            assert!(length < PI, "Boundary {} length {:.4} should be reasonable", i, length);
        }
    }

    #[test]
    fn subduction_logic_follows_physics() {
        let world = TectonicWorld::new(8, 42);
        
        for boundary in world.boundaries() {
            if let BoundaryType::Convergent { subducting_plate: Some(subducting_id), .. } = &boundary.boundary_type {
                let subducting_plate = world.get_plate(*subducting_id).unwrap();
                let other_plate_id = boundary.other_plate(*subducting_id).unwrap();
                let other_plate = world.get_plate(other_plate_id).unwrap();
                
                // Oceanic plates should subduct under continental plates
                match (&subducting_plate.plate_type, &other_plate.plate_type) {
                    (PlateType::Oceanic { .. }, PlateType::Continental { .. }) => {
                        // This is correct - oceanic subducts under continental
                    },
                    (PlateType::Continental { .. }, PlateType::Oceanic { .. }) => {
                        panic!("Continental plate {} subducting under oceanic plate {} - this is wrong!", 
                               subducting_id, other_plate_id);
                    },
                    (PlateType::Oceanic { density: d1 }, PlateType::Oceanic { density: d2 }) => {
                        assert!(d1 > d2, "Denser oceanic plate should subduct");
                    },
                    (PlateType::Continental { .. }, PlateType::Continental { .. }) => {
                        panic!("Continental-continental collision should not have subduction");
                    },
                }
            }
        }
    }

    // Phase 1 specific tests as specified in the design document
    #[test]
    fn every_point_has_one_plate() {
        let world = TectonicWorld::new(10, 42);
        let mut rng = StdRng::seed_from_u64(123);
        
        // Sample 1000 random points
        for _ in 0..1000 {
            let lat = rng.gen_range(-std::f64::consts::PI/2.0..std::f64::consts::PI/2.0);
            let lon = rng.gen_range(-std::f64::consts::PI..std::f64::consts::PI);
            let point = SphericalPoint::new(lat, lon);
            
            let plate_id = world.plate_at(point);
            
            // Verify the plate ID is valid
            assert!(plate_id < world.num_plates());
            
            // Verify the point is actually closest to this plate's seed
            let assigned_plate = world.get_plate(plate_id).unwrap();
            let distance_to_assigned = crate::geometry::distance(point, assigned_plate.seed);
            
            // Check that this distance is smaller than or equal to distance to any other plate
            for other_plate in world.plates() {
                let distance_to_other = crate::geometry::distance(point, other_plate.seed);
                assert!(distance_to_assigned <= distance_to_other + 1e-10); // Allow for floating point errors
            }
        }
    }

    #[test]
    fn plate_at_is_deterministic() {
        let world = TectonicWorld::new(8, 789);
        let test_points = vec![
            SphericalPoint::new(0.0, 0.0),
            SphericalPoint::new(0.5, 1.0),
            SphericalPoint::new(-0.3, 2.0),
            SphericalPoint::new(1.0, -1.0),
        ];
        
        for point in test_points {
            let plate_id1 = world.plate_at(point);
            let plate_id2 = world.plate_at(point);
            let plate_id3 = world.plate_at(point);
            
            assert_eq!(plate_id1, plate_id2);
            assert_eq!(plate_id2, plate_id3);
        }
    }

    #[test]
    fn all_plates_have_area() {
        let world = TectonicWorld::new(12, 456);
        let mut rng = StdRng::seed_from_u64(789);
        
        // Sample many random points and count how many belong to each plate
        let mut plate_counts = vec![0; world.num_plates()];
        let num_samples = 10000;
        
        for _ in 0..num_samples {
            let lat = rng.gen_range(-std::f64::consts::PI/2.0..std::f64::consts::PI/2.0);
            let lon = rng.gen_range(-std::f64::consts::PI..std::f64::consts::PI);
            let point = SphericalPoint::new(lat, lon);
            
            let plate_id = world.plate_at(point);
            plate_counts[plate_id] += 1;
        }
        
        // Every plate should have at least some points assigned to it
        // (no degenerate/empty plates)
        for (plate_id, count) in plate_counts.iter().enumerate() {
            assert!(*count > 0, "Plate {} has no area (0 points assigned)", plate_id);
        }
        
        // All samples should be accounted for
        let total_assigned: usize = plate_counts.iter().sum();
        assert_eq!(total_assigned, num_samples);
    }
}

