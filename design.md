# Tectonic Plate System Implementation Plan

## Executive Summary

Build a Voronoi-based tectonic plate simulation library for a fractal-scale Earth simulation. The system uses spherical Voronoi tessellation with evolving boundary complexity, integrates with an existing LOD chunk system (LOD 0 = 12-voxel diameter Earth, LOD 1 = 194-voxel diameter), and supports narrative generation.

## Core Design Principles

1. **Dual Representation**: Voronoi seeds (topology) + boundary polylines (geometry)
2. **Spherical Coordinates**: Work on unit sphere, adapter converts to/from Cartesian WorldPosition
3. **Evolutionary Complexity**: Boundaries start simple, gain points over geological time
4. **LOD-Aware Detail**: LOD 0 = simplified schematic (6-12 points), LOD 1 = full detail (100-500 points)
5. **Separation of Concerns**: Tectonics library separate from climate/erosion

---

## Phase 1: Core Data Structures & Voronoi

**Goal**: Get basic plate ownership queries working

### Types to Implement
```rust
// Core coordinates
struct SphericalPoint { lat: f64, lon: f64 }

// Plate definition
struct TectonicPlate {
    id: PlateId,
    seed: SphericalPoint,
    motion: PlateMotion,
    plate_type: PlateType,
    mass: f64,
    age: f64,
}

struct PlateMotion {
    euler_pole: SphericalPoint,
    angular_velocity: f64,
}

enum PlateType {
    Oceanic { density: f64 },
    Continental { density: f64 },
}

// World container
struct TectonicWorld {
    plates: Vec<TectonicPlate>,
    voronoi: SphericalVoronoi,
}
```

### Voronoi Implementation
```rust
struct SphericalVoronoi {
    seeds: Vec<SphericalPoint>,
}

impl SphericalVoronoi {
    // Convert spherical to 3D Cartesian, find closest seed
    fn cell_index(&self, point: SphericalPoint) -> usize;
}
```

### Key Functions
- `SphericalPoint::to_cartesian() -> [f64; 3]`
- `SphericalPoint::from_cartesian([f64; 3]) -> Self`
- `distance(SphericalPoint, SphericalPoint) -> f64` (great circle distance)
- `TectonicWorld::new(num_plates: usize, seed: u64) -> Self`
- `TectonicWorld::plate_at(point: SphericalPoint) -> PlateId`

### Tests (Phase 1)
```rust
#[test]
fn every_point_has_one_plate() {
    // Sample 1000 random points, verify closest seed assignment
}

#[test]
fn plate_at_is_deterministic() {
    // Same point always returns same plate
}

#[test]
fn all_plates_have_area() {
    // No degenerate/empty plates
}
```

**Deliverable**: Can create world with N plates, query which plate any point belongs to.

---

## Phase 2: Boundary Generation

**Goal**: Generate boundary polylines between adjacent plates

### Types to Add
```rust
enum BoundaryType {
    Divergent { spreading_rate: f64 },
    Convergent { 
        subduction: Option<SubductionZone>,
        orogeny: Option<MountainBuilding>,
    },
    Transform { stress: f64 },
}

struct PlateBoundary {
    plate_a: PlateId,
    plate_b: PlateId,
    boundary_type: BoundaryType,
    geometry: Vec<SphericalPoint>,  // Start simple: 20-50 points
}
```

### Boundary Generation Algorithm
```rust
impl TectonicWorld {
    fn generate_boundaries(&mut self) {
        // For each pair of adjacent plates:
        // 1. Find Voronoi edge between them
        // 2. Sample points along edge (20-50 points initially)
        // 3. Determine boundary type from plate motions
        // 4. Store as PlateBoundary
    }
    
    fn classify_boundary_type(
        plate_a: &TectonicPlate,
        plate_b: &TectonicPlate,
        boundary_point: SphericalPoint
    ) -> BoundaryType {
        // Calculate relative motion at boundary
        // Divergent: moving apart (dot product negative)
        // Convergent: moving together (dot product positive)
        // Transform: sliding (perpendicular motion)
    }
}
```

### Tests (Phase 2)
```rust
#[test]
fn boundaries_separate_exactly_two_plates() {
    // Points on boundary equidistant to both plates
}

#[test]
fn every_plate_has_boundaries() {
    // No isolated plates
}

#[test]
fn boundary_network_is_connected() {
    // Can reach any plate from any other via boundaries
}
```

**Deliverable**: Can generate and query boundaries between all adjacent plates.

---

## Phase 3: Plate Motion & Physics

**Goal**: Plates move according to Euler pole rotation

### Motion Implementation
```rust
impl TectonicWorld {
    fn tick(&mut self, dt: f64) {
        // For each plate:
        // 1. Rotate seed point around Euler pole
        // 2. Update boundary geometries
        // 3. Detect collisions/interactions
        // 4. Update boundary types based on new relative motions
    }
}

fn rotate_point_around_euler_pole(
    point: SphericalPoint,
    euler_pole: SphericalPoint,
    angle: f64
) -> SphericalPoint {
    // Convert to Cartesian
    // Create rotation matrix around Euler pole axis
    // Rotate and convert back
}
```

### Physics Rules
```rust
fn determine_subduction(plate_a: &TectonicPlate, plate_b: &TectonicPlate) 
    -> Option<PlateId> 
{
    match (plate_a.plate_type, plate_b.plate_type) {
        (Oceanic{..}, Continental{..}) => Some(plate_a.id), // Oceanic subducts
        (Continental{..}, Oceanic{..}) => Some(plate_b.id),
        (Oceanic{d1}, Oceanic{d2}) => {
            if d1 > d2 { Some(plate_a.id) } else { Some(plate_b.id) } // Denser subducts
        },
        (Continental{..}, Continental{..}) => None, // Mountain building instead
    }
}
```

### Tests (Phase 3)
```rust
#[test]
fn plates_move_according_to_euler_pole() {
    // Verify rotation math
}

#[test]
fn oceanic_subducts_under_continental() {
    // Set up collision, verify subduction assignment
}

#[test]
fn continental_collision_creates_mountains() {
    // No subduction, mountain building instead
}

#[test]
fn area_conservation() {
    // Total surface area remains 4π after motion
}
```

**Deliverable**: Plates move realistically, boundaries update, physics rules enforced.

---

## Phase 4: Boundary Evolution & Complexity

**Goal**: Boundaries gain fractal detail over time

### Evolution Algorithm
```rust
impl PlateBoundary {
    fn evolve(&mut self, dt: f64) {
        for i in 0..self.geometry.len()-1 {
            let segment_length = distance(geometry[i], geometry[i+1]);
            let max_length = self.max_segment_length();
            
            if segment_length > max_length {
                // Insert midpoint with noise
                let midpoint = interpolate(geometry[i], geometry[i+1], 0.5);
                let displaced = self.apply_noise(midpoint);
                self.geometry.insert(i+1, displaced);
            }
        }
    }
    
    fn max_segment_length(&self) -> f64 {
        match self.boundary_type {
            Divergent => 50_000.0,   // Rifts stay smooth
            Transform => 10_000.0,    // Transform faults get jagged
            Convergent => 20_000.0,   // Moderate complexity
        }
    }
    
    fn apply_noise(&self, point: SphericalPoint) -> SphericalPoint {
        // Perlin/Simplex noise perpendicular to boundary
        // Scale by boundary age and type
    }
}
```

### Tests (Phase 4)
```rust
#[test]
fn boundaries_gain_complexity_over_time() {
    // Point count increases over many ticks
}

#[test]
fn transform_boundaries_more_jagged_than_divergent() {
    // Measure fractal dimension
}

#[test]
fn old_plates_have_complex_boundaries() {
    // Age correlates with complexity
}
```

**Deliverable**: Boundaries evolve realistic fractal detail based on geological activity.

---

## Phase 5: LOD System & Simplification

**Goal**: Support different detail levels for different scales

### LOD Implementation
```rust
impl PlateBoundary {
    // Store multiple LOD versions
    simplified_lod0: Vec<SphericalPoint>,  // 6-12 points
    geometry: Vec<SphericalPoint>,          // Full detail (100-500+)
    
    fn geometry_at_lod(&self, lod: u8) -> &[SphericalPoint] {
        match lod {
            0 => &self.simplified_lod0,
            1..=7 => &self.geometry,
            _ => &self.geometry,
        }
    }
    
    fn update_lod_representations(&mut self) {
        // Douglas-Peucker simplification for LOD 0
        self.simplified_lod0 = douglas_peucker(
            &self.geometry,
            epsilon: 500_000.0  // 500km tolerance
        );
    }
}

// Douglas-Peucker line simplification
fn douglas_peucker(
    points: &[SphericalPoint],
    epsilon: f64
) -> Vec<SphericalPoint> {
    // Recursive simplification algorithm
    // Keep points that deviate > epsilon from straight line
}
```

### Tests (Phase 5)
```rust
#[test]
fn lod0_boundaries_are_simplified() {
    // LOD 0 has 6-12 points per boundary
}

#[test]
fn lod1_boundaries_have_full_detail() {
    // LOD 1 preserves all evolved complexity
}

#[test]
fn simplification_preserves_topology() {
    // Simplified lines don't cross or change connectivity
}
```

**Deliverable**: Can query boundaries at appropriate detail for rendering scale.

---

## Phase 6: Integration Adapter

**Goal**: Connect to existing chunk system

### Adapter Implementation
```rust
struct TectonicChunkAdapter {
    world: TectonicWorld,
    earth_radius: f64,  // 6_400_000.0 meters
}

impl TectonicChunkAdapter {
    fn world_pos_to_spherical(&self, pos: [f64; 3]) -> SphericalPoint {
        // Cartesian [x,y,z] from Earth center -> lat/lon
        let r = (pos[0].powi(2) + pos[1].powi(2) + pos[2].powi(2)).sqrt();
        let lat = (pos[2] / r).asin();
        let lon = pos[1].atan2(pos[0]);
        SphericalPoint { lat, lon }
    }
    
    fn spherical_to_world_pos(&self, point: SphericalPoint) -> [f64; 3] {
        // Lat/lon -> Cartesian on Earth surface
        let x = self.earth_radius * point.lat.cos() * point.lon.cos();
        let y = self.earth_radius * point.lat.cos() * point.lon.sin();
        let z = self.earth_radius * point.lat.sin();
        [x, y, z]
    }
    
    fn elevation_for_voxel(&self, world_pos: [f64; 3]) -> f64 {
        let spherical = self.world_pos_to_spherical(world_pos);
        
        // Base elevation from distance to nearest boundary
        let dist_to_boundary = self.world.distance_to_boundary(spherical);
        
        // Modify based on boundary type
        // Convergent: mountains (high elevation)
        // Divergent: rifts (low elevation)
        // Transform: neutral
        
        // Add noise for terrain variation
    }
    
    fn narrative_for_chunk(&self, chunk_coord: ChunkCoord) -> NarrativeSummary {
        // Sample chunk center, find plate and nearby boundaries
        // Generate description based on tectonic activity
    }
}
```

### Tests (Phase 6)
```rust
#[test]
fn lod0_chunk_contains_earth() {
    // All Earth surface points within chunk (0,0,0) bounds
}

#[test]
fn coordinate_conversion_round_trips() {
    // world_pos -> spherical -> world_pos preserves position
}

#[test]
fn elevation_continuous_across_chunks() {
    // No discontinuities at chunk boundaries
}

#[test]
fn narrative_reflects_plate_boundaries() {
    // Convergent boundaries mention mountains/volcanism
}
```

**Deliverable**: Tectonic system fully integrated with existing chunk architecture.

---

## Phase 7: Comprehensive Testing

**Goal**: Complete test coverage per TDD design

### Test Categories

**Invariant Tests** (must always hold):
- Every point belongs to exactly one plate
- Boundaries separate exactly two plates
- Total surface area conserved (4π)
- No degenerate plates
- Boundary network connected

**Physics Tests** (realism):
- Oceanic subducts under continental
- Motion follows Euler pole
- Divergent boundaries create spreading
- Mountain building at continental collisions

**Evolution Tests** (emergence):
- Complexity increases over time
- Transform boundaries more jagged than divergent
- Age correlates with complexity

**Integration Tests** (system compatibility):
- LOD 0 chunk contains Earth
- Elevation continuous across chunks
- Narrative reflects tectonic activity
- Deterministic for same seed

**Snapshot Tests** (regression):
- Deterministic evolution
- Known good state preservation

### Test Helpers to Implement
```rust
fn random_spherical_point() -> SphericalPoint;
fn calculate_fractal_dimension(points: &[SphericalPoint]) -> f64;
fn is_graph_connected(world: &TectonicWorld) -> bool;
fn assert_worlds_equal(a: &TectonicWorld, b: &TectonicWorld);
```

**Deliverable**: Full test suite, all tests passing.

---

## Implementation Notes

### Dependencies
- `rand` - Random number generation
- `noise` - Perlin/Simplex noise for boundary evolution
- `serde` - Serialization for snapshots (optional)

### Performance Considerations
- Voronoi query: O(n) naive, optimize with spatial indexing later
- Boundary evolution: Only evolve boundaries near areas of interest
- LOD simplification: Cache simplified versions, update periodically

### Future Extensions (Not Phase 1-7)
- Spatial indexing (k-d tree) for faster Voronoi queries
- Microplates (hierarchical Voronoi)
- Heat flux calculation for climate integration
- Volcanic activity zones
- Earthquake simulation at transform faults

---

## Success Criteria

After Phase 7 complete:
1. ✅ Can create 10-20 plate Earth simulation
2. ✅ Plates move realistically via Euler pole rotation
3. ✅ Boundaries show correct interaction types (divergent/convergent/transform)
4. ✅ Boundaries evolve fractal complexity over geological time
5. ✅ LOD 0 shows simplified schematic (6-12 points per boundary)
6. ✅ LOD 1 shows full detail (100-500+ points per boundary)
7. ✅ Integrates with existing chunk system via adapter
8. ✅ All invariant, physics, evolution, and integration tests pass
9. ✅ Deterministic simulation for same seed
10. ✅ Can generate narrative summaries for chunks based on tectonic activity

---

## Phase Priorities

**Must Complete**: Phases 1-6 (core functionality + integration)
**Should Complete**: Phase 7 (comprehensive testing)
**Nice to Have**: Future extensions

Each phase should be completed with its tests passing before moving to next phase. Use TDD approach: write tests first, then implement to make tests pass.