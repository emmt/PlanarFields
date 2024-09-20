# PlanarFields [![Build Status](https://github.com/emmt/PlanarFields.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/emmt/PlanarFields.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Coverage](https://codecov.io/gh/emmt/PlanarFields.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/emmt/PlanarFields.jl)

`PlanarFields` is a Julia package to deal with fields sampled on 2-dimensional grids of
equally spaced nodes. This package extends the types defined in
[`TwoDimensional`](https://github.com/emmt/TwoDimensional.jl)

The following types of object are provided:

- `GridAxis{T}` for abstract ranges whose values of type `T` are multiple of a step.

- `Grid{T}` for abstract 2-dimensional arrays of equally spaced nodes with coordinates of
  type `T` and elements of type `Point{T}`.

- `PlanarField{T,S}` for abstract 2-dimensional arrays storing samples of a field in a
  grid of nodes. `T` and `S` are the respective types of the stored values and of the grid
  step.

A typical usage of `PlanarField` objects is to represent a wave in transverse planes
relative to the direction of propagation as well as transmission masks whose properties
depend on the position.
