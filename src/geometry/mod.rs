//! Geometric data structures and algorithms for spherical geometry
//! 
//! This module contains all the geometric primitives and algorithms used
//! in the tectonic simulation, including:
//! - Point representations on a sphere
//! - Distance calculations
//! - Voronoi tessellation for spatial queries
//! - Complete Voronoi diagram construction

pub mod point;
pub mod voronoi_query;
pub mod voronoi_diagram;

// Re-export core types from point module
pub use point::{
    PlateId,
    SphericalPoint,
    PlateMotion,
    PlateType,
    TectonicPlate,
    distance,
    generate_random_seeds,
};

// Re-export voronoi query types
pub use voronoi_query::SphericalVoronoi;

// Re-export voronoi diagram types
pub use voronoi_diagram::{
    VoronoiVertex,
    VoronoiEdge,
    SphericalVoronoiDiagram,
};

