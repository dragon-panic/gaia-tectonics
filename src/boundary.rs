//! Boundary types and structures for tectonic plate interactions

use crate::geometry::{PlateId, SphericalPoint, distance};
use noise::{NoiseFn, Perlin};

/// Type of tectonic boundary between two plates
#[derive(Debug, Clone, PartialEq)]
pub enum BoundaryType {
    /// Divergent boundary where plates move apart
    /// Creates spreading centers, rifts, and new oceanic crust
    Divergent { 
        /// Spreading rate in meters per year
        spreading_rate: f64 
    },
    /// Convergent boundary where plates move together
    /// Can result in subduction or mountain building
    Convergent { 
        /// Which plate subducts (if any)
        subducting_plate: Option<PlateId>,
        /// Whether mountain building occurs
        mountain_building: bool,
    },
    /// Transform boundary where plates slide past each other
    /// Creates strike-slip faults
    Transform { 
        /// Shear stress magnitude
        stress: f64 
    },
}

impl BoundaryType {
    /// Get a human-readable description of the boundary type
    pub fn description(&self) -> &'static str {
        match self {
            BoundaryType::Divergent { .. } => "Divergent (Spreading)",
            BoundaryType::Convergent { .. } => "Convergent (Collision)",
            BoundaryType::Transform { .. } => "Transform (Sliding)",
        }
    }

    /// Get the spreading rate for divergent boundaries
    pub fn spreading_rate(&self) -> Option<f64> {
        match self {
            BoundaryType::Divergent { spreading_rate } => Some(*spreading_rate),
            _ => None,
        }
    }

    /// Get the stress for transform boundaries
    pub fn stress(&self) -> Option<f64> {
        match self {
            BoundaryType::Transform { stress } => Some(*stress),
            _ => None,
        }
    }

    /// Check if this is a convergent boundary with subduction
    pub fn has_subduction(&self) -> bool {
        match self {
            BoundaryType::Convergent { subducting_plate, .. } => subducting_plate.is_some(),
            _ => false,
        }
    }

    /// Check if this is a convergent boundary with mountain building
    pub fn has_mountain_building(&self) -> bool {
        match self {
            BoundaryType::Convergent { mountain_building, .. } => *mountain_building,
            _ => false,
        }
    }
}

/// A boundary between two tectonic plates
#[derive(Debug, Clone)]
pub struct PlateBoundary {
    /// First plate involved in the boundary
    pub plate_a: PlateId,
    /// Second plate involved in the boundary
    pub plate_b: PlateId,
    /// Type of boundary interaction
    pub boundary_type: BoundaryType,
    /// Geometry of the boundary as a series of points
    /// Points are ordered along the boundary
    pub geometry: Vec<SphericalPoint>,
    /// Age of the boundary in time units
    pub age: f64,
    /// Noise generator for boundary evolution
    noise: Perlin,
}

impl PlateBoundary {
    /// Create a new plate boundary
    pub fn new(
        plate_a: PlateId,
        plate_b: PlateId,
        boundary_type: BoundaryType,
        geometry: Vec<SphericalPoint>,
    ) -> Self {
        Self {
            plate_a,
            plate_b,
            boundary_type,
            geometry,
            age: 0.0,
            noise: Perlin::new(42), // Use deterministic seed
        }
    }

    /// Get the other plate ID given one plate ID
    pub fn other_plate(&self, plate_id: PlateId) -> Option<PlateId> {
        if plate_id == self.plate_a {
            Some(self.plate_b)
        } else if plate_id == self.plate_b {
            Some(self.plate_a)
        } else {
            None
        }
    }

    /// Check if this boundary involves the given plate
    pub fn involves_plate(&self, plate_id: PlateId) -> bool {
        self.plate_a == plate_id || self.plate_b == plate_id
    }

    /// Get the length of the boundary in radians (approximate)
    pub fn length(&self) -> f64 {
        if self.geometry.len() < 2 {
            return 0.0;
        }

        let mut total_length = 0.0;
        for i in 0..self.geometry.len() - 1 {
            total_length += crate::geometry::distance(self.geometry[i], self.geometry[i + 1]);
        }
        total_length
    }

    /// Get the number of points in the boundary geometry
    pub fn point_count(&self) -> usize {
        self.geometry.len()
    }

    /// Evolve the boundary by adding fractal detail over time
    pub fn evolve(&mut self, dt: f64, enable_noise: bool) {
        self.age += dt;
        
        // Only evolve if we have enough points
        if self.geometry.len() < 2 {
            return;
        }
        
        // PHASE 1: Subdivide long segments first (before applying noise)
        // This ensures we don't have large gaps even from initial Voronoi geometry
        let max_length = self.max_segment_length();
        self.subdivide_long_segments(max_length);
        
        if !enable_noise {
            return;
        }
        
        // PHASE 2: Apply noise displacement to interior points
        // This creates continuous deformation along the entire boundary
        // CRITICAL: Never modify endpoints - they must stay at triple junctions
        for i in 1..self.geometry.len() - 1 {
            let displaced = self.apply_noise(self.geometry[i], i);
            self.geometry[i] = displaced;
        }
        
        // PHASE 3: Subdivide again after noise application
        // Noise can stretch segments, so we need to check again
        self.subdivide_long_segments(max_length);
    }
    
    /// Subdivide segments that exceed the maximum length
    fn subdivide_long_segments(&mut self, max_length: f64) {
        let mut i = 0;
        let mut total_subdivisions = 0;
        let max_subdivisions_per_step = 100; // Limit to prevent infinite loops
        
        while i < self.geometry.len() - 1 && total_subdivisions < max_subdivisions_per_step {
            let segment_length = distance(self.geometry[i], self.geometry[i + 1]);
            
            // Always subdivide segments that exceed the maximum length
            if segment_length > max_length {
                // Insert midpoint using spherical interpolation
                let midpoint = self.interpolate_spherical(self.geometry[i], self.geometry[i + 1], 0.5);
                self.geometry.insert(i + 1, midpoint);
                total_subdivisions += 1;
                // Don't increment i - check both new segments
            } else {
                i += 1;
            }
        }
    }

    /// Get the maximum allowed segment length based on boundary type
    fn max_segment_length(&self) -> f64 {
        // Convert from meters to radians (approximate)
        // Earth radius ≈ 6,400,000 meters
        let earth_radius_m = 6_400_000.0;
        
        // Keep segments short to prevent gaps from developing
        // These values ensure no segment can exceed ~100km even with noise
        let base_threshold = match &self.boundary_type {
            BoundaryType::Divergent { .. } => 80_000.0 / earth_radius_m,   // ~80km - Rifts stay smooth
            BoundaryType::Transform { .. } => 60_000.0 / earth_radius_m,   // ~60km - Transform faults get jagged  
            BoundaryType::Convergent { .. } => 70_000.0 / earth_radius_m,  // ~70km - Moderate complexity
        };
        
        // Gradually decrease threshold with age to add more detail over time
        // Start at base_threshold, asymptote to 50% of base_threshold
        let age_factor = 0.5 + 0.5 * (-self.age * 0.05).exp();
        base_threshold * age_factor
    }

    /// Apply noise to a point to create fractal detail
    /// Uses multi-scale noise for realistic tectonic-style deformation
    fn apply_noise(&self, point: SphericalPoint, index: usize) -> SphericalPoint {
        let cart = point.to_cartesian();
        
        // Multi-scale noise: combine large-scale and small-scale features
        // This creates more realistic, self-similar fractal boundaries
        let scale_large = 5.0;
        let scale_medium = 15.0;
        let scale_small = 40.0;
        
        let noise_large = self.noise.get([cart[0] * scale_large, cart[1] * scale_large, cart[2] * scale_large]);
        let noise_medium = self.noise.get([cart[0] * scale_medium, cart[1] * scale_medium, cart[2] * scale_medium]);
        let noise_small = self.noise.get([cart[0] * scale_small, cart[1] * scale_small, cart[2] * scale_small]);
        
        // Combine with different weights (fractal behavior: 1/f weighting)
        let combined_noise = noise_large * 0.5 + noise_medium * 0.3 + noise_small * 0.2;
        
        // Base noise scale depends on boundary type
        // These values are kept small to prevent creating gaps between adjacent points
        let base_scale = match &self.boundary_type {
            BoundaryType::Divergent { .. } => 0.0008,    // Very smooth rifts
            BoundaryType::Transform { stress } => 0.002 * (1.0 + stress * 0.5), // Jagged, stress-dependent
            BoundaryType::Convergent { .. } => 0.0012,   // Moderate complexity
        };
        
        // Age factor: boundaries become more complex over time
        // Start slow, then accelerate, then plateau
        let age_factor = ((self.age * 0.05).tanh() * 2.0).min(1.5);
        
        // Apply displacement with temporal variation (simulate tectonic activity pulses)
        let time_variation = (self.age * 0.3 + index as f64 * 0.1).sin() * 0.3 + 0.7;
        let noise_magnitude = combined_noise * base_scale * age_factor * time_variation;
        
        // Calculate perpendicular displacement direction
        // Find local tangent by looking at neighboring points
        let tangent = if index > 0 && index < self.geometry.len() - 1 {
            let prev = self.geometry[index - 1].to_cartesian();
            let next = self.geometry[index + 1].to_cartesian();
            [next[0] - prev[0], next[1] - prev[1], next[2] - prev[2]]
        } else {
            // Fallback for endpoints
            [1.0, 0.0, 0.0]
        };
        
        // Calculate perpendicular direction (cross product with position vector)
        let perp = [
            tangent[1] * cart[2] - tangent[2] * cart[1],
            tangent[2] * cart[0] - tangent[0] * cart[2],
            tangent[0] * cart[1] - tangent[1] * cart[0],
        ];
        
        // Normalize perpendicular vector
        let perp_mag = (perp[0] * perp[0] + perp[1] * perp[1] + perp[2] * perp[2]).sqrt();
        let perp_norm = if perp_mag > 0.0001 {
            [perp[0] / perp_mag, perp[1] / perp_mag, perp[2] / perp_mag]
        } else {
            [0.0, 1.0, 0.0] // Fallback
        };
        
        // Apply displacement in perpendicular direction
        let displaced_cart = [
            cart[0] + perp_norm[0] * noise_magnitude,
            cart[1] + perp_norm[1] * noise_magnitude,
            cart[2] + perp_norm[2] * noise_magnitude,
        ];
        
        // Re-normalize to stay on unit sphere
        let mag = (displaced_cart[0] * displaced_cart[0] + 
                   displaced_cart[1] * displaced_cart[1] + 
                   displaced_cart[2] * displaced_cart[2]).sqrt();
        let normalized_cart = [
            displaced_cart[0] / mag,
            displaced_cart[1] / mag,
            displaced_cart[2] / mag,
        ];
        
        SphericalPoint::from_cartesian(normalized_cart)
    }


    /// Interpolate between two spherical points
    fn interpolate_spherical(&self, a: SphericalPoint, b: SphericalPoint, t: f64) -> SphericalPoint {
        // Convert to Cartesian for interpolation
        let cart_a = a.to_cartesian();
        let cart_b = b.to_cartesian();
        
        // Interpolate in Cartesian space
        let interpolated = [
            cart_a[0] + t * (cart_b[0] - cart_a[0]),
            cart_a[1] + t * (cart_b[1] - cart_a[1]),
            cart_a[2] + t * (cart_b[2] - cart_a[2]),
        ];
        
        // Convert back to spherical
        SphericalPoint::from_cartesian(interpolated)
    }

    /// Get the age of the boundary
    pub fn get_age(&self) -> f64 {
        self.age
    }

    /// Get the complexity measure (number of points per unit length)
    pub fn complexity(&self) -> f64 {
        let length = self.length();
        if length > 0.0 {
            self.point_count() as f64 / length
        } else {
            0.0
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn boundary_type_descriptions() {
        let divergent = BoundaryType::Divergent { spreading_rate: 2.0 };
        let convergent = BoundaryType::Convergent { 
            subducting_plate: Some(0), 
            mountain_building: true 
        };
        let transform = BoundaryType::Transform { stress: 1.5 };

        assert_eq!(divergent.description(), "Divergent (Spreading)");
        assert_eq!(convergent.description(), "Convergent (Collision)");
        assert_eq!(transform.description(), "Transform (Sliding)");
    }

    #[test]
    fn boundary_type_properties() {
        let divergent = BoundaryType::Divergent { spreading_rate: 2.0 };
        let convergent = BoundaryType::Convergent { 
            subducting_plate: Some(1), 
            mountain_building: true 
        };
        let transform = BoundaryType::Transform { stress: 1.5 };

        assert_eq!(divergent.spreading_rate(), Some(2.0));
        assert_eq!(transform.stress(), Some(1.5));
        assert!(convergent.has_subduction());
        assert!(convergent.has_mountain_building());
        assert!(!divergent.has_subduction());
        assert!(!transform.has_mountain_building());
    }

    #[test]
    fn plate_boundary_operations() {
        let geometry = vec![
            SphericalPoint::new(0.0, 0.0),
            SphericalPoint::new(0.1, 0.0),
            SphericalPoint::new(0.2, 0.0),
        ];
        let boundary = PlateBoundary::new(
            0, 1,
            BoundaryType::Divergent { spreading_rate: 2.0 },
            geometry,
        );

        assert_eq!(boundary.other_plate(0), Some(1));
        assert_eq!(boundary.other_plate(1), Some(0));
        assert_eq!(boundary.other_plate(2), None);
        assert!(boundary.involves_plate(0));
        assert!(boundary.involves_plate(1));
        assert!(!boundary.involves_plate(2));
        assert_eq!(boundary.point_count(), 3);
    }

    #[test]
    fn plate_boundary_length() {
        let geometry = vec![
            SphericalPoint::new(0.0, 0.0),
            SphericalPoint::new(0.1, 0.0),
            SphericalPoint::new(0.2, 0.0),
        ];
        let boundary = PlateBoundary::new(
            0, 1,
            BoundaryType::Divergent { spreading_rate: 2.0 },
            geometry,
        );

        let length = boundary.length();
        assert!(length > 0.0);
        assert!(length < 1.0); // Should be reasonable for small distances
    }

    #[test]
    fn empty_boundary_length() {
        let boundary = PlateBoundary::new(
            0, 1,
            BoundaryType::Divergent { spreading_rate: 2.0 },
            vec![],
        );

        assert_eq!(boundary.length(), 0.0);
        assert_eq!(boundary.point_count(), 0);
    }

    #[test]
    fn boundary_evolution_adds_complexity() {
        let mut boundary = PlateBoundary::new(
            0, 1,
            BoundaryType::Transform { stress: 1.0 },
            vec![
                SphericalPoint::new(0.0, 0.0),
                SphericalPoint::new(0.1, 0.0), // Large segment
            ],
        );

        let initial_points = boundary.point_count();
        
        // Evolve over time
        boundary.evolve(1.0, true);
        boundary.evolve(1.0, true);
        boundary.evolve(1.0, true);
        
        // Should have more points due to subdivision
        assert!(boundary.point_count() > initial_points);
        assert!(boundary.get_age() > 0.0);
    }

    #[test]
    fn transform_boundaries_more_complex_than_divergent() {
        let mut transform_boundary = PlateBoundary::new(
            0, 1,
            BoundaryType::Transform { stress: 1.0 },
            vec![
                SphericalPoint::new(0.0, 0.0),
                SphericalPoint::new(0.1, 0.0),
            ],
        );

        let mut divergent_boundary = PlateBoundary::new(
            0, 1,
            BoundaryType::Divergent { spreading_rate: 2.0 },
            vec![
                SphericalPoint::new(0.0, 0.0),
                SphericalPoint::new(0.1, 0.0),
            ],
        );

        // Evolve both boundaries
        for _ in 0..5 {
            transform_boundary.evolve(1.0, true);
            divergent_boundary.evolve(1.0, true);
        }

        // Transform boundaries should have more points (more complex)
        assert!(transform_boundary.point_count() >= divergent_boundary.point_count());
    }

    #[test]
    fn boundary_age_increases_over_time() {
        let mut boundary = PlateBoundary::new(
            0, 1,
            BoundaryType::Convergent { 
                subducting_plate: Some(0), 
                mountain_building: true 
            },
            vec![
                SphericalPoint::new(0.0, 0.0),
                SphericalPoint::new(0.1, 0.0),
            ],
        );

        assert_eq!(boundary.get_age(), 0.0);
        
        boundary.evolve(2.5, true);
        assert_eq!(boundary.get_age(), 2.5);
        
        boundary.evolve(1.0, true);
        assert_eq!(boundary.get_age(), 3.5);
    }

    #[test]
    fn boundary_complexity_measurement() {
        let mut boundary = PlateBoundary::new(
            0, 1,
            BoundaryType::Transform { stress: 1.0 },
            vec![
                SphericalPoint::new(0.0, 0.0),
                SphericalPoint::new(0.1, 0.0),
            ],
        );

        let initial_complexity = boundary.complexity();
        
        // Evolve to add more points
        for _ in 0..3 {
            boundary.evolve(1.0, true);
        }
        
        // Complexity should increase (more points per unit length)
        assert!(boundary.complexity() >= initial_complexity);
    }
}
