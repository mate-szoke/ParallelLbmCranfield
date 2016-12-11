# Parallel LBM solver Project

This software is a two dimensional parallel Lattice Boltzmann Method solver implemented in CUDA and PGAS UPC

## Abstract

The Unified Parallel C (UPC) language, as a member of the Partitioned Global Address Space (PGAS) family, unifies the advantages of shared and local memory spaces and offers a relatively straight forward code parallelisation with the Central Processing Unit (CPU). In contrast, the Computer Unified Device Architecture (CUDA) development kit gives a tool to make use of the Graphics Processing Unit (GPU) and manage its complex memory structure. We provide a detailed comparison between these novel techniques through the parallelisation of a two-dimensional lattice Boltzmann method based incompressible fluid flow solver. Our comparison between the CUDA and UPC parallelisation takes into account the required conceptual effort, the performance gain, and the limitations of the approaches from the application oriented developers’ point of view. UPC was proven to be applied for complex physical problems, as the local memory implementation led to competitive efficiency. However, the performance of the shared memory code fell behind our expectations, and we concluded that the investigated compilers could not treat efficiently the shared memory. The CUDA implementation proved to be more cumbersome compared to the UPC approach. However, the results confirmed the formerly documented experience of the community, based on which the structure of the graphics card makes it highly suitable for the parallelisation of the lattice Boltzmann method. Our work may serve as a guidance for programmers to avoid stumbling blocks during the parallelisation of the lattice Boltzmann method. Furthermore, the presented performance measurements might help to understand and identify the bottle neck of the method depending on the chosen parallelisation technique and architecture. We hope that the compiler developers will find our feedback useful, and they will be able to further develop the promising shared memory approach of the UPC language.

## Project structure
- **cases**: Mesh data for a standart lid driven cavity with different fineness
- **cuda**: CUDA source code for the different phases of the project
- **serial**: Initial serial code from which the parallel code was derived.
- **upc**: PGAS UPC implementation for local and shared memory model

For building instructions look for the README in the individual directories.

## Authors
Mate Tibor Szoke

Tamas Istvan Jozsa

Adam Koleszar
