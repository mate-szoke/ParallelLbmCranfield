# CUDA LBM SOLVER
## Cranfield University 2015
This software is a Lattice Boltzmann solver for simple geometries.

# Build
To build the program use the provided Makefile. By default it will create object files and build the code with single precision (float)

    $ make


To build with double precision

    $ make FLOAT_TYPE=USE_DOUBLE


To build a more compact (and slightly faster) version without object files

    $ make release


To create Doxygen documentation (HTML, LaTeX)

    $ make doc


To create Doxygen documentation and pdf

    $ made latexdoc


To run unittests

    $ make test


To create a debug build

    $ make debug


To clean the project

    $ make clean

## Build on Windows
The code is prepared to be built in Windows environment although success is not guaranteed. To build
it you will need GNU make. You can use the Windows compliant MinGW version or the POSIX compliant
Cygwin version, both should work without other changes the same way on Linux.
 * Visual C++ compiler suite (added to PATH) / also work with the free Express edition
 * CUDA Toolkit (added to PATH)
 * GNU make (added to PATH)

# Dependencies
## CuTest
CuTest is a single file ANSI C unittest framework, which compiles together with the code so it does
not need any special treatment.
## Argtable 2
Argtable is an ANSI C command line parser, whick use CMAKE for building configuration thus can be
built on any environment. On windows it does not support regex so some parameters are defined as
normal string parameters on windows. To build it, download and setup CMAKE than generate the
project according to your environment and build it.

The repository contains already built version of the library for 64bit linux and 64bit Windows.

# Running the program
## Input files
The program needs at least 3 files to run if no parameters passed to the executable
 * SetUpData.ini
 * Mesh/D2node.dat
 * Mesh/BCconnectors.dat

## Input parameters
To get the full list of parameters

    $ ./lbmsolver -h
    Usage: ./lbmsolver  [-ht] [-f <file>] [-n <file>] [-b <file>] [-o <file>] [-u <u>] [-v <v>] [-r <rho>] [-s <nu>] [--inlet=[yes|no|pulsatile]] [-c [BGKW|TRT|MRT]] [--curved] [-l [yes|second|first]] [-i <N>] [--every=<N>] [--after=<N>] [--format=[paraview|tecplot]] [-d <id>]
    Usage: ./lbmsolver  compare <file> [<file>]

| short | long                          | description                                                 |
|-------|-------------------------------|-------------------------------------------------------------|
| -h    | --help                        | Print help options                                          |
| -f    | --initfile=\<file\>           | Initialisation from file (default: SetUpData.ini)           |
| -t    | --test                        | Run unit tests                                              |
| -n    | --node=\<file\>               | Node file (default: Mesh/D2node.dat)                        |
| -b    | --bc=\<file\>                 | Boundary conditions file (default: Mesh/BCconnectors.dat)   |
| -o    | --output=\<file\>             | Output directory (default: Results)                         |
| -u    | --uavg=\<u\>                  | Mean U (x velocity)                                         |
| -v    | --vavg=\<v\>                  | Mean V (y velocity)                                         |
| -r    | --rho=\<rho\>                 | Density                                                     |
| -s    | --viscosity=\<nu\>            | Viscosity                                                   |
|       | --inlet=[yes\|no\|pulsatile]  | Inlet profile (default: no)                                 |
| -c    | --collision=[BGKW\|TRT\|MRT]  | Collision model (default: BGKW)                             |
|       | --curved                      | Curved boundaries                                           |
| -l    | --outlet=[yes\|second\|first] | Outlet profile (default: second)                            |
| -i    | --iter=\<N\>                  | Number of iterations (default: 1000)                        |
|       | --every=\<N\>                 | Autosave after every \<N\> iterations (default: 0)          |
|       | --after=\<N\>                 | Start autosaving after the \<N\>th iteration (default: 1000)|
|       | --format=[paraview\|tecplot]  | Output format (default: paraview)                           |
| -d    | --draglift=\<id\>             | Calculate drag/lift on \<id\> boundary (default: 0)         |

If an init file is passed with -f all other parameters are discarded.

## Compare files
To compare results in the default place with previous ones just run

    $ ./lbmsolver compare <path to result file>

To compare files at path different than the default

    $ ./lbmsolver compare <path to result file1> <path to result file2>

# Authors
Adam Koleszar

Tamas Istvan Jozsa

Mate Tibor Szoke