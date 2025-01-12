# HydroMech3D

**HydroMech3D** is a cutting-edge 3D Hydro-Mechanical Boundary Element Method (BEM) code designed to study fault stability and poromechanics during injection activities related to **geothermal energy** and **carbon (CO₂) storage**. Written in **C++** and powered by **Armadillo**, a high-performance linear algebra library, this code incorporates advanced modeling techniques such as **rate-and-state friction laws**, **fluid flow simulations**, and **GPU-accelerated computations** using **NVIDIA cuDSS Sparse Linear Systems Solver**.

---

## Table of Contents
1. [Introduction](#introduction)
2. [Features](#features)
3. [Applications](#applications)
4. [Requirements](#requirements)
5. [Installation](#installation)
6. [Usage](#usage)
7. [Code Structure](#code-structure)
8. [Acknowledgements](#acknowledgements)
9. [References](#references)

---

## Introduction

As the world strives to combat climate change, **Carbon Capture and Storage (CCS)** has emerged as a crucial technology for achieving net-zero emissions. CCS involves capturing CO₂ from industrial processes or directly from the atmosphere and securely storing it underground. This technology supports global decarbonization strategies by enabling the storage of CO₂ in deep geological formations, such as depleted oil fields, saline aquifers, and basalt reservoirs.

The safety and efficacy of CCS depend on accurate predictions of CO₂ migration and pressure buildup over hundreds of years. Traditional flow simulators face computational limitations, making the high-throughput screening of storage sites and injection schemes impractical. HydroMech3D addresses these challenges by providing a highly efficient, GPU-accelerated framework for simulating fault mechanics and poromechanics, ensuring the safe and effective storage of CO₂.

---

## Features

- **3D Hydro-Mechanical BEM**: Models coupled hydro-mechanical processes in subsurface reservoirs.
- **Rate-and-State Friction Laws**: Simulates fault stability under dynamic stress conditions.
- **Fluid Flow Modeling**: Accounts for fluid injection and its impact on reservoir pressure.
- **GPU Acceleration**: Leverages **NVIDIA cuDSS** for efficient sparse linear system solvers.
- **High Performance**: Built using **Armadillo**, ensuring fast and robust matrix computations.
- **Modular Design**: Flexible code structure, enabling easy modification and extension.

---

## Applications

HydroMech3D is designed for a variety of use cases, including:
1. **Carbon Capture and Storage (CCS)**:
   - Predicting CO₂ plume migration.
   - Managing reservoir pressure to prevent fault reactivation or seismic hazards.
2. **Geothermal Energy**:
   - Studying fluid injection impacts on fault stability.
3. **Seismic Hazard Assessment**:
   - Evaluating induced seismicity from fluid injection activities.
4. **Subsurface Resource Management**:
   - Optimizing well placement and injection rates.

---

## Requirements

To run HydroMech3D, the following software and hardware are required:

- **Operating System**: Linux, macOS, or Windows.
- **C++ Compiler**: GCC or Clang (supporting C++11 or later).
- **Armadillo Library**: For linear algebra and matrix operations.
- **NVIDIA GPU** (optional, for acceleration): Requires CUDA toolkit.
- **CMake**: For building the project.

---

## Installation

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/lucadalzilio/HydroMech3D.git
   cd HydroMech3D