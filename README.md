# Parallel LBM solver Project

This software is a two dimensional parallel Lattice Boltzmann Method solver implemented in CUDA and PGAS UPC

## Abstract

The Unified Parallel C (UPC) language from the Partitioned Global Address Space (PGAS) family unifies the advantages of shared and local memory spaces and offers a relatively straight forward code parallelisation with the Central Processing Unit (CPU). In contrast, the Computer Unified Device Architecture (CUDA) development kit gives a tool to make use of the Graphics Processing Unit (GPU). We provide a detailed comparison between these novel techniques through the parallelisation of a two-dimensional lattice Boltzmann method based fluid flow solver. Our comparison between the CUDA and UPC parallelisation takes into account the required conceptual effort, the performance gain, and the limitations of the approaches from the application oriented developers’ point of view. We demonstrated that UPC led to competitive efficiency with the local memory implementation. However, the performance of the shared memory code fell behind our expectations, and we concluded that the investigated UPC compilers could not treat efficiently the shared memory space. The CUDA implementation proved to be more complex compared to the UPC approach mainly because of the complicated memory structure of the graphics card which also makes GPUs suitable for the parallelisation of the lattice Boltzmann method.

## Project structure

- **cases**: Mesh data for a standart lid driven cavity with different fineness
- **cuda**: CUDA source code for the different phases of the project
- **serial**: Initial serial code from which the parallel code was derived.
- **upc**: PGAS UPC implementation for local and shared memory model

For building instructions look for the README in the individual directories.

## Authors
Máté Tibor Szőke

Tamás István Józsa

Ádám Koleszár

Irene Moulitsas

László Könözsy
