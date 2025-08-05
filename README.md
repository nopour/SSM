
![A](https://github.com/nopour/SSM/blob/main/Figure1.png)

# SSM Reduction for a 2-DOF Nonlinear Shaw-Pierre System

[![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/navendu-pottekkat/awesome-readme?include_prereleases)](https://img.shields.io/github/v/release/navendu-pottekkat/awesome-readme?include_prereleases)
[![GitHub](https://img.shields.io/github/license/navendu-pottekkat/awesome-readme)](https://img.shields.io/github/license/navendu-pottekkat/awesome-readme)

This script applies the equation-driven Spectral Submanifold (SSM) model reduction methodology to a classic 2-DOF spring-mass-damper system with a cubic nonlinearity. It performs the reduction from first principles using MATLAB's Symbolic Math Toolbox, without relying on the external SSMTool package.

# WorkFlow


![GuideLine](https://github.com/nopour/SSM/blob/main/Figure2.png)

% 1. System Definition: Construct the 4D Full Order Model (FOM).

% 2. Modal Analysis: Find the slow and fast subspaces of the linear system.

% 3. Symbolic SSM Calculation: Derive the 3rd-order manifold z = h(y).

% 4. Symbolic ROM Formulation: Derive the 2D nonlinear reduced dynamics.

% 5. Numerical Simulation: Solve both FOM and ROM ODEs.

% 6. Visualization

# License

[MIT license](./LICENSE)
![GuideLine](https://github.com/nopour/SSM/blob/main/Figure3.png)

