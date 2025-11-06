# Antenna Array Design (MATLAB + CVX)

This repository contains MATLAB scripts for antenna array weight design using an MM (Majorization–Minimization) epigraph formulation. It supports both continuous-modulus and discrete-phase constraints and includes a sweep of the penalty parameter to compare performance.

## Overview

- Core script: `mm_epgigraph.m`
- Array: Uniform Linear Array (ULA)
- Method: MM epigraph optimization to maximize beamforming power over a target angular sector
- Constraints:
  - Continuous modulus (`|x_i| <= 1`)
  - Discrete phase (L-PSK via convex hull and projection)
- Outputs: Diagnostic figures and selected weights saved to a `.mat` file

## Prerequisites

- MATLAB (R2018a or later recommended)
- [CVX](http://cvxr.com/cvx/) convex optimization toolbox installed and added to MATLAB path
  - Follow CVX installation instructions and run `cvx_setup` once

## Quick Start

1. Ensure CVX is installed and configured:
   ```matlab
   % Only needed once after installing CVX
   cvx_setup
   ```
2. Open and run the main script:
   ```matlab
   % From MATLAB
   open mm_epgigraph.m
   % Click Run, or execute directly:
   run('mm_epgigraph.m');
   ```

The script will generate beam pattern figures and save a results file like `lambda=<index>.mat` containing the selected weights.

## Key Parameters (in `mm_epgigraph.m`)

- `N`: number of transmit antennas (default `64`)
- `width`: target sector width in degrees (default `60`)
- `target_DoA`: boresight angle in degrees (default `0`)
- `fc`: carrier frequency (default `3.2e9` Hz)
- `spacing`: element spacing (default `0.5 * lambda`)
- `L`: discrete phase levels for PSK (default `4`)
- `ll`: lambda sweep vector for penalty term, e.g. `0.1:0.02:3`

You can modify these at the top-level of the script to explore different array sizes, sectors, and discrete quantization levels.

## Important Functions

The main script references helper functions included in the repository:

- `getH.m`: builds structured vectors/matrices for target sector handling
- `getw.m`: transforms complex vector to real representation used in MM update
- `projection.m`: projects continuous solution onto discrete phase set `s`
- `projection_con.m`: projects onto constant modulus (continuous phase)

Ensure these files are present in the same folder when running the script.

## Outputs

- Beam pattern plots for:
  - Continuous modulus solution
  - Discrete phase (before/after projection)
  - Comparison plots across lambda sweep
- Results file:
  - `lambda=<index>.mat` with variables `w_cont`, `w_before`, `w_after`, `w_after_c`

## Notes

- CVX runs silently (`cvx_begin quiet`) in the script. If you need debug output, remove `quiet`.
- Large `.mat`, `.fig`, and `.eps` files are typical byproducts; consider excluding them from Git version control using a `.gitignore`.

## Acknowledgments

- CVX: Grant and Boyd (http://cvxr.com/cvx/)
