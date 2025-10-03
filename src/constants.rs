//! Constants used throughout the tectonic simulation

/// Epsilon value for floating point comparisons
pub const EPSILON: f64 = 1e-9;

/// Larger epsilon for deduplicating vertices
/// (vertices within this distance are considered the same)
pub const DEDUPLICATION_EPSILON: f64 = EPSILON * 10000.0;

/// Target segment length in radians (~75km on Earth)
/// This ensures no boundary segment exceeds evolution thresholds
pub const TARGET_SEGMENT_LENGTH_RADIANS: f64 = 0.012;

/// Maximum segment length before subdivision (in radians)
/// Used in boundary evolution to determine when to add detail
pub const MAX_SEGMENT_LENGTH: f64 = 0.1;

/// Typical oceanic plate density (g/cm³)
pub const OCEANIC_DENSITY: f64 = 3.0;

/// Typical continental plate density (g/cm³)
pub const CONTINENTAL_DENSITY: f64 = 2.7;

