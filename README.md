# Spline tie-decay temporal networks

This repository provides the codes for constructing kernel functions of spline tie-decay networks used in\
[Chanon Thongprayoon, Naoki Masuda. Preprint arXiv:2408.11913](https://arxiv.org/abs/2408.11913).

## Python packages
- `import sympy as sp`

## Variable set up
- `ta`: event arrival time
- `a`: kernel's value at time `ta`
- `h`: time delay
- `e1`: kernel's derivative at time `ta`
- `k`: kernel's value at time `ta+h`
- `e2`: kernel's derivative at time `ta+h`
  
## Classes
- `CubicSpline`: creates cubic spline polynomials.
- `ExponentialDecay`: a subclass of `CubicSpline`. Creates exponential decays.
