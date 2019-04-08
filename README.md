# FC Tables

----
FC Tables is a Fortran 90 program that computes FC-Gram tables. 

This program heavily utilizes the [MPFUN2015](https://www.davidhbailey.com/dhbsoftware/) library for arbitrary precision floating point arithmetic. This library comes in two variants and both are included. The first one, **MPFUN-Fort**, is a pure Fortran implementation, which is straightforward to compile. The second variant, **MPFUN-MPFR**, depends on the C libraries **GMP** and **MPFR** and, although trickier to compile, delivers about 5x faster performance. It is, however, very simply to choose which variant of MPFUN should FC Tables  try compile.

FC Tables can also take advantage of multicore processors by utilizing OpenMP to split loops across multiple threads.

More information about the FC-Gram algorithm to perform periodic continuations of data can be found at https://doi.org/10.1016/j.jcp.2009.11.020 and references therein.

----
## Compilation
1. Clone/download this repository.
2. `cd` into the project directory.
3. Edit the text file `Makefile.in` to set the options for the desired precision, MPFUN variant,  OpenMP usage, compilers and flags.
4. Execute the `make` command.

## Running FC Tables

1. Edit `parameter.inp` and select the desired table (i.e. a number of matching points —`d`— and continuation points —`c`—) and the output directory (`odir`).
2. Run `./fc_tables` in the project directory.

## Veryfing the tables
For convenience a simple program that compares numerical (using FC-Gram continuations) to analytical derivatives of the function $`J_0(x)e^x`$ is also supplied. 

To compile the test program simply run, after compiling `fc_tables`, the command `make test` and execute the newly generated binary with `./fc_test N`, where `N` the number of grid points to use (which must be a power of two due to the `fft` implementation)

The necessary table is automatically loaded based on the values present in the `parameter.inp` file.


### Author
Mauro Fontana - Department of Physics - University of Buenos Aires.
