# Tectonics

A Voronoi-based tectonic plate simulation library for spherical surfaces.

[![Crates.io](https://img.shields.io/crates/v/tectonics.svg)](https://crates.io/crates/tectonics)
[![Documentation](https://docs.rs/tectonics/badge.svg)](https://docs.rs/tectonics)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/License-MIT%20OR%20Apache--2.0-blue.svg)](LICENSE)

## Overview

This library implements a spherical Voronoi tessellation system for simulating tectonic plate dynamics. Plates move according to Euler pole rotation mechanics, automatically generating and classifying boundaries (convergent, divergent, transform) that evolve realistic fractal detail over geological time.

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
tectonics = "0.1"
```

## Try the GUI Demo

Want to see the simulation in action? Download pre-built binaries from the [latest release](https://github.com/dragon-panic/tectonics/releases/latest):

### Downloads

| Platform | Architecture | Download |
|----------|-------------|----------|
| 🐧 Linux | x86_64 | `simple_gui-linux-x86_64` |
| 🪟 Windows | x86_64 | `simple_gui-windows-x86_64.exe` |
| 🍎 macOS | Intel (x86_64) | `simple_gui-macos-x86_64` |
| 🍎 macOS | Apple Silicon (ARM64) | `simple_gui-macos-arm64` |

### Running the GUI Demo

**Linux/macOS:**
```bash
# Make the binary executable
chmod +x simple_gui-*

# Run it
./simple_gui-linux-x86_64
# or
./simple_gui-macos-x86_64
```

**Windows:**
```cmd
# Just double-click or run from command line
simple_gui-windows-x86_64.exe
```

The GUI provides:
- Real-time plate motion visualization
- Interactive simulation controls (play/pause/step)
- Multiple projection modes (orthographic globe view, equirectangular flat map)
- Live plate statistics and boundary analysis
- Toggle features on/off to explore the simulation

## Quick Start

```rust
use tectonics::{TectonicWorld, SphericalPoint};

// Create a world with 12 tectonic plates
let mut world = TectonicWorld::new(12, 42);

// Query which plate owns a point
let point = SphericalPoint::new(0.0, 0.0);
let plate_id = world.plate_at(point);

// Simulate 1 million years
world.tick(1e6, true);

// Examine the boundaries
for boundary in world.boundaries() {
    println!("Boundary type: {:?}", boundary.boundary_type);
}
```

## Features

- 🌍 **Spherical Voronoi tessellation** - Accurate plate boundaries on a sphere
- 🔄 **Euler pole rotation** - Realistic plate motion mechanics
- 🏔️ **Boundary classification** - Automatic detection of convergent/divergent/transform boundaries
- 📈 **Fractal evolution** - Boundaries gain geological detail over time
- 🎯 **Efficient queries** - Fast spatial lookups and plate ownership
- 🔬 **Physics validation** - Area conservation and motion accuracy

## Core Concepts

### Tectonic Plates

Each plate has:
- **Type**: Oceanic or Continental (with density)
- **Motion**: Euler pole rotation (axis and angular velocity)
- **Properties**: Mass, age, seed position

```rust
let config = WorldConfig::new()
    .with_number_of_plates(15)
    .with_continental_fraction(0.3);  // 30% continental plates

let world = TectonicWorld::with_config(42, config);
```

### Plate Boundaries

Three types of boundaries, automatically classified from plate motions:

- **Divergent**: Plates moving apart → Mid-ocean ridges, rifts
- **Convergent**: Plates colliding → Subduction zones, mountain building
- **Transform**: Plates sliding past → Strike-slip faults

```rust
for boundary in world.boundaries() {
    match boundary.boundary_type {
        BoundaryType::Convergent { subducting_plate, mountain_building } => {
            if mountain_building {
                println!("Mountains forming!");
            }
        }
        BoundaryType::Divergent { spreading_rate } => {
            println!("Spreading at {} m/year", spreading_rate);
        }
        BoundaryType::Transform { stress } => {
            println!("Transform fault with stress {}", stress);
        }
    }
}
```

### Boundary Evolution

Boundaries gain fractal complexity over time using Perlin noise:

```rust
// Advance simulation by 1 million years with boundary evolution
world.tick(1e6, true);

// Check boundary complexity
for boundary in world.boundaries() {
    println!("Boundary age: {} years", boundary.get_age());
    println!("Complexity: {:.2} points/radian", boundary.complexity());
}
```

## Key Features

### Spherical Voronoi Tessellation
- Accurate plate boundaries on unit sphere
- Efficient spatial queries for plate ownership
- Deterministic generation from seed values

### Euler Pole Rotation
- Physically accurate plate motion
- Configurable rotation axis and angular velocity
- Automatic Voronoi updates as plates move

### Boundary Classification
- Automatic detection based on relative plate motions
- Physics-based subduction logic (oceanic plates subduct under continental)
- Continental collision detection (mountain building)

### Fractal Detail Evolution
- Boundaries gain complexity over geological time
- Type-specific characteristics (transform boundaries are most jagged)
- Perlin noise-based deformation for realistic appearance

### Physics Validation
- Area conservation (total surface area remains 4π)
- Motion accuracy verification
- Comprehensive test suite (57 tests)

## Running Examples

```bash
# Run the basic usage example
cargo run --example basic_usage

# Run the GUI visualization (requires GUI feature)
cargo run --example simple_gui --features gui

# Run all tests
cargo test

# Run tests with verbose output
cargo test -- --nocapture
```

### GUI Visualization

The `simple_gui` example provides an interactive visualization:

- **Real-time Simulation**: Watch plates move and boundaries evolve
- **Interactive Controls**: Step through time, adjust time step, pause/resume
- **Plate Statistics**: View plate types, positions, motion, and age
- **Boundary Analysis**: See boundary types, complexity, and evolution
- **Visual Representation**: Colored circles for plates, colored lines for boundaries
- **Plate Selection**: Click "Select" to highlight specific plates

## Dependencies

- **Core**: `rand`, `noise`
- **GUI examples** (optional): `eframe`, `egui` (enabled with `--features gui`)

## Use Cases

This library is designed for:
- **Procedural world generation** in games and simulations
- **Geological modeling** and visualization
- **Educational tools** for plate tectonics
- **Research** in spherical geometry and spatial algorithms

For terrain generation (heightmaps, elevation), consider building on top of this library or using a companion crate like `tectonics-terrain`.

## Design Principles

- **Deterministic**: Same seed always produces same world
- **Spherical accuracy**: All calculations on unit sphere (not flat projections)
- **Physically motivated**: Euler pole rotation, subduction rules based on density
- **Separation of concerns**: Pure plate dynamics, no terrain/climate coupling
- **Evolving complexity**: Boundaries gain detail over geological time
