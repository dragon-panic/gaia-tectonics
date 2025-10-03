//! Core data structures for the tectonic plate system

use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;

/// Unique identifier for a tectonic plate
pub type PlateId = usize;

/// A point on the unit sphere in spherical coordinates
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SphericalPoint {
    /// Latitude in radians (-π/2 to π/2)
    pub lat: f64,
    /// Longitude in radians (-π to π)
    pub lon: f64,
}

impl SphericalPoint {
    /// Create a new spherical point from latitude and longitude in radians
    pub fn new(lat: f64, lon: f64) -> Self {
        Self { lat, lon }
    }

    /// Convert spherical coordinates to 3D Cartesian coordinates on unit sphere
    pub fn to_cartesian(&self) -> [f64; 3] {
        let x = self.lat.cos() * self.lon.cos();
        let y = self.lat.cos() * self.lon.sin();
        let z = self.lat.sin();
        [x, y, z]
    }

    /// Convert 3D Cartesian coordinates to spherical coordinates
    pub fn from_cartesian(cartesian: [f64; 3]) -> Self {
        let [x, y, z] = cartesian;
        let r = (x * x + y * y + z * z).sqrt();
        
        // Normalize to unit sphere
        let x = x / r;
        let y = y / r;
        let z = z / r;
        
        let lat = z.asin();
        let lon = y.atan2(x);
        
        Self { lat, lon }
    }
}

/// Motion parameters for a tectonic plate
#[derive(Debug, Clone, Copy)]
pub struct PlateMotion {
    /// Euler pole (rotation axis) in spherical coordinates
    pub euler_pole: SphericalPoint,
    /// Angular velocity in radians per unit time
    pub angular_velocity: f64,
}

impl PlateMotion {
    /// Create a new plate motion with given Euler pole and angular velocity
    pub fn new(euler_pole: SphericalPoint, angular_velocity: f64) -> Self {
        Self {
            euler_pole,
            angular_velocity,
        }
    }
}

/// Type of tectonic plate (oceanic or continental)
#[derive(Debug, Clone, Copy)]
pub enum PlateType {
    /// Oceanic plate with given density
    Oceanic { density: f64 },
    /// Continental plate with given density
    Continental { density: f64 },
}

impl PlateType {
    /// Get the density of the plate
    pub fn density(&self) -> f64 {
        match self {
            PlateType::Oceanic { density } => *density,
            PlateType::Continental { density } => *density,
        }
    }
}

/// A tectonic plate with its properties and motion
#[derive(Debug, Clone)]
pub struct TectonicPlate {
    /// Unique identifier for this plate
    pub id: PlateId,
    /// Seed point for Voronoi tessellation
    pub seed: SphericalPoint,
    /// Motion parameters
    pub motion: PlateMotion,
    /// Type of plate (oceanic/continental)
    pub plate_type: PlateType,
    /// Mass of the plate (for physics calculations)
    pub mass: f64,
    /// Age of the plate in time units
    pub age: f64,
}

impl TectonicPlate {
    /// Create a new tectonic plate
    pub fn new(
        id: PlateId,
        seed: SphericalPoint,
        motion: PlateMotion,
        plate_type: PlateType,
        mass: f64,
        age: f64,
    ) -> Self {
        Self {
            id,
            seed,
            motion,
            plate_type,
            mass,
            age,
        }
    }
}

/// Calculate the great circle distance between two points on the unit sphere
pub fn distance(a: SphericalPoint, b: SphericalPoint) -> f64 {
    let [x1, y1, z1] = a.to_cartesian();
    let [x2, y2, z2] = b.to_cartesian();
    
    // Dot product of unit vectors gives cosine of angle between them
    let dot_product = x1 * x2 + y1 * y2 + z1 * z2;
    
    // Clamp to avoid numerical errors
    let dot_product = dot_product.clamp(-1.0, 1.0);
    
    // Arc cosine gives the angle in radians
    dot_product.acos()
}

/// Generate random seed points on the unit sphere
/// Uses a uniform distribution over the sphere surface
pub fn generate_random_seeds(num_seeds: usize, seed: u64) -> Vec<SphericalPoint> {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut seeds = Vec::with_capacity(num_seeds);

    for _ in 0..num_seeds {
        // Generate uniform random points on unit sphere
        // Using the method: generate random normal vector, then normalize
        let x: f64 = rng.gen_range(-1.0..1.0);
        let y: f64 = rng.gen_range(-1.0..1.0);
        let z: f64 = rng.gen_range(-1.0..1.0);
        
        let r = (x * x + y * y + z * z).sqrt();
        if r > 0.0 {
            let cartesian = [x / r, y / r, z / r];
            seeds.push(SphericalPoint::from_cartesian(cartesian));
        } else {
            // Fallback to a known point if we get the origin
            seeds.push(SphericalPoint::new(0.0, 0.0));
        }
    }

    seeds
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn spherical_to_cartesian_round_trip() {
        let original = SphericalPoint::new(0.5, 1.0);
        let cartesian = original.to_cartesian();
        let converted = SphericalPoint::from_cartesian(cartesian);
        
        // Should be approximately equal (within floating point precision)
        assert!((original.lat - converted.lat).abs() < 1e-10);
        assert!((original.lon - converted.lon).abs() < 1e-10);
    }

    #[test]
    fn cartesian_to_spherical_round_trip() {
        // Use a unit vector for the test
        let original = [0.5_f64, 0.3, 0.8];
        let r = (original[0] * original[0] + original[1] * original[1] + original[2] * original[2]).sqrt();
        let unit_vector = [original[0] / r, original[1] / r, original[2] / r];
        
        let spherical = SphericalPoint::from_cartesian(unit_vector);
        let converted = spherical.to_cartesian();
        
        // Should be approximately equal (within floating point precision)
        for i in 0..3 {
            assert!((unit_vector[i] - converted[i]).abs() < 1e-10);
        }
    }

    #[test]
    fn distance_is_symmetric() {
        let a = SphericalPoint::new(0.5, 1.0);
        let b = SphericalPoint::new(-0.3, 2.0);
        
        assert_eq!(distance(a, b), distance(b, a));
    }

    #[test]
    fn distance_to_self_is_zero() {
        let point = SphericalPoint::new(0.5, 1.0);
        assert_eq!(distance(point, point), 0.0);
    }

    #[test]
    fn distance_between_antipodal_points_is_pi() {
        let a = SphericalPoint::new(0.0, 0.0);
        let b = SphericalPoint::new(0.0, PI);
        
        assert!((distance(a, b) - PI).abs() < 1e-10);
    }

    #[test]
    fn plate_type_density() {
        let oceanic = PlateType::Oceanic { density: 3.0 };
        let continental = PlateType::Continental { density: 2.7 };
        
        assert_eq!(oceanic.density(), 3.0);
        assert_eq!(continental.density(), 2.7);
    }
}
