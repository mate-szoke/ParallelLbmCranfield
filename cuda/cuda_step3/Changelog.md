Changelog
=========

Notable changes to the CUDA LBM Solver project

## [1.8.2] - 2015-08-09
### Added
  * Add combined, faster TRT kernel

### Changed
  * rename variable and function names

### Fixes
  * Fix TRT collision model
  * Fix wall condition

## [1.8.1] - 2015-07-25
### Added
  * Add support for building on Windows
  * Add resource freeing to failing tests (teardown)
  * Add more unittests
  * Add result check for CUDA API calls

### Changed
  * Split residual step into norm and drag/lift computation

### Fixed
  * Fix failing unittests
  * FIx unittests not checking anything

## [1.8.0] - 2015-07-20
### Added
  * Add doc and latexdoc target to Makefile
  * Add Doxygen documentation comments
  * Add new struct for input parameters
  * Add this changelog

### Changed
  * Refactor variable names
  * Separate headers
  * Separate input handling from main

## [1.7.1] - 2015-07-10
### Added
  * O3 to release target (to little effect)
  * Add comparing unittests for boundary conditions

### Fixed
  * Fix outlet handling in macro step

## [1.7.0] - 2015-07-06
### Added
  * Add command line interface

### Fixed
  * Fix problem with float/double switch

### Changed
  * Makefile options for float/double

## [1.6.1] - 2015-07-03
### Changed
  * Combined MRT steps
  * Rearrange MRT collision model to use less resurces

### Added
  * Add release target to Makefile
  * Add debug target to Makefile

## [1.6.0] - 2015-07-02
### Changed
  * Unfold kernels for sum, streaming and collision

### Added
  * Add utility functions for unittests

## [1.5.0] - 2015-06-29
### Changed
  * Rearranged boundary conditions to speed up process and save memory

### Added
  * Add unittests for init and boundary conditions
  * Add option to compare result files

## [1.4.0] - 2015-06-24
### Added
  * GPU residual computation using GPU reduce (sum)
  * Add unittests

### Changed
  * refactor init loop

## [1.3.0] - 2015-06-16
### Changed
  * Generalise float type (use macros)

### Removed
  * Remove gpu_update_f and integrate it into gpu_streaming

## [1.2.0] - 2015-06-16
### Changed
  * Move GPU constants to initialisation

### Added
  * Add final results write

## [1.1.0] - 2015-06-15
### Removed
  * Removed binaries and meshes from repository

### Added
  * Add header files

### Changed
  * Move GPU functions to separate source files
  * Speed up initialisation by rearranging the init loop

## [1.0.1] - 2015-06-11
### Fixed
  * Fix typo in Makefile

## [1.0] - 2015-06-11

### Added
  * Create repository
  * Fork lbm-solver
  * Add Makefile