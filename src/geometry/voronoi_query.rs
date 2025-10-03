//! Spherical Voronoi tessellation implementation for fast point queries

use super::point::{SphericalPoint, distance};

/// Spherical Voronoi tessellation for plate ownership queries
#[derive(Debug, Clone)]
pub struct SphericalVoronoi {
    /// Seed points for the Voronoi cells
    pub seeds: Vec<SphericalPoint>,
}

impl SphericalVoronoi {
    /// Create a new spherical Voronoi with the given seed points
    pub fn new(seeds: Vec<SphericalPoint>) -> Self {
        Self { seeds }
    }

    /// Find the index of the Voronoi cell that contains the given point
    /// Returns the index of the closest seed point
    pub fn cell_index(&self, point: SphericalPoint) -> usize {
        if self.seeds.is_empty() {
            panic!("Cannot query empty Voronoi diagram");
        }

        let mut closest_index = 0;
        let mut closest_distance = distance(point, self.seeds[0]);

        for (i, seed) in self.seeds.iter().enumerate().skip(1) {
            let dist = distance(point, *seed);
            if dist < closest_distance {
                closest_distance = dist;
                closest_index = i;
            }
        }

        closest_index
    }


    /// Get the number of seed points
    pub fn len(&self) -> usize {
        self.seeds.len()
    }

    /// Check if the Voronoi diagram is empty
    pub fn is_empty(&self) -> bool {
        self.seeds.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::generate_random_seeds;

    #[test]
    fn cell_index_returns_valid_index() {
        let seeds = vec![
            SphericalPoint::new(0.0, 0.0),
            SphericalPoint::new(0.5, 0.0),
            SphericalPoint::new(-0.5, 0.0),
        ];
        let voronoi = SphericalVoronoi::new(seeds);
        
        let point = SphericalPoint::new(0.1, 0.0);
        let index = voronoi.cell_index(point);
        
        assert!(index < voronoi.len());
    }

    #[test]
    fn cell_index_finds_closest_seed() {
        let seeds = vec![
            SphericalPoint::new(0.0, 0.0),    // index 0
            SphericalPoint::new(1.0, 0.0),    // index 1
            SphericalPoint::new(-1.0, 0.0),   // index 2
        ];
        let voronoi = SphericalVoronoi::new(seeds);
        
        // Point closer to seed 0
        let point = SphericalPoint::new(0.1, 0.0);
        let index = voronoi.cell_index(point);
        assert_eq!(index, 0);
        
        // Point closer to seed 1
        let point = SphericalPoint::new(0.6, 0.0);
        let index = voronoi.cell_index(point);
        assert_eq!(index, 1);
    }

    #[test]
    fn generate_random_seeds_produces_correct_count() {
        let seeds = generate_random_seeds(10, 42);
        assert_eq!(seeds.len(), 10);
    }

    #[test]
    fn generate_random_seeds_is_deterministic() {
        let seeds1 = generate_random_seeds(5, 123);
        let seeds2 = generate_random_seeds(5, 123);
        
        assert_eq!(seeds1.len(), seeds2.len());
        for (s1, s2) in seeds1.iter().zip(seeds2.iter()) {
            assert!((s1.lat - s2.lat).abs() < 1e-10);
            assert!((s1.lon - s2.lon).abs() < 1e-10);
        }
    }

    #[test]
    fn generate_random_seeds_different_seeds_produce_different_results() {
        let seeds1 = generate_random_seeds(5, 123);
        let seeds2 = generate_random_seeds(5, 456);
        
        // Should be different (very unlikely to be identical)
        let mut different = false;
        for (s1, s2) in seeds1.iter().zip(seeds2.iter()) {
            if (s1.lat - s2.lat).abs() > 1e-10 || (s1.lon - s2.lon).abs() > 1e-10 {
                different = true;
                break;
            }
        }
        assert!(different);
    }

    #[test]
    fn empty_voronoi_panics_on_query() {
        let voronoi = SphericalVoronoi::new(vec![]);
        let point = SphericalPoint::new(0.0, 0.0);
        
        // This should panic
        let result = std::panic::catch_unwind(|| {
            voronoi.cell_index(point)
        });
        assert!(result.is_err());
    }

    #[test]
    fn voronoi_length_and_is_empty() {
        let empty_voronoi = SphericalVoronoi::new(vec![]);
        assert!(empty_voronoi.is_empty());
        assert_eq!(empty_voronoi.len(), 0);
        
        let seeds = vec![SphericalPoint::new(0.0, 0.0)];
        let non_empty_voronoi = SphericalVoronoi::new(seeds);
        assert!(!non_empty_voronoi.is_empty());
        assert_eq!(non_empty_voronoi.len(), 1);
    }
}

