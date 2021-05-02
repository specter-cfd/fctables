# FC Tables

----
FC Tables is a Fortran 2003 program that computes FC-Gram tables.

This program heavily utilizes the 
[MPFUN2020](https://www.davidhbailey.com/dhbsoftware/) library for arbitrary 
precision floating point arithmetic. This library comes  in two variants and 
both are included. The first one, **MPFUN-Fort**, is a pure Fortran 
implementation, which is straightforward to compile. The second variant, 
**MPFUN-MPFR**, depends on the C libraries **GMP** and **MPFR** and, although 
trickier to compile, delivers about 3x faster performance. It is, however, 
very simple to choose which variant of MPFUN should FC Tables try compile.

FC Tables can also take advantage of multicore processors by utilizing OpenMP 
to split loops across multiple threads.

More information about the FC-Gram algorithm to perform periodic continuations 
of data can be found at, https://doi.org/10.1016/j.cpc.2020.107482 and 
https://doi.org/10.1016/j.jcp.2009.11.020, and references therein.

----
## Compilation
1. Clone/download this repository.
2. `cd` into the project directory.
3. Edit the text file `Makefile.in` to set 
   1. The desired precision in decimal digits, `DIGITS`.
   2. Whether to build support for constructing blend-to-zero operators 
   (`BLEND = yes`) or not (`BLEND = no`).
   3. The `MPFUN` variant, `fort` or `mpfr`.
   4. Whether OpenMP support should be enabled (`OPENMP = yes`) or 
   not (`OPENMP = no`).
   5. Compiler family to use. Options are `GNU` for `gfortran` and `INTEL`
   for `ifort`. Note that building with `ifort` has not been thoroughly tested.
4. Execute the `make` command.
 
Note that if `MPFUN = mpfr` is selected, the first compilation can take up to 
20/30 minutes, as `mpfr` and `gmp` (its dependency) have very long compilation 
times. These libraries are only compiled one time unless `make distclean` 
is invoked.


## Usage
1. Edit `parameter.inp` and select the desired settings.
2. Run `./fc_tables` in the project directory.

To run in `N` cores, `OMP_NUM_THREADS=N ./fc_tables` should be used instead.


## Runtime options
The file `parameter.inp` is splitted in three sections, namely `required`, 
`optional` and `svd`. The `optional` and `svd` sections are only used if 
computing blend-to-zero operators `A`.
 
| Namelist   |  Parameter |   Default value    |         Description           |
|:----------:|:----------:|:------------------:|-------------------------------|
| `required` |   `odir`   |                    | Directory to write the `Q` and `A` arrays.|
| `required` |    `d`     |                    | Number of matching points to use in the periodic extension. Note that `d` points generate a `d-1` order extension |
| `required` |    `C`     |                    | Number of continuation points to use in the periodic extension.|
| `required` |  `tkind`   |                    | `=0` to construct Dirichlet projector `Q` and, if compiled with `BLEND=yes`, blend-to-zero operator `A`.`=1` to construct Neumann projector. `=2` to construct second normal derivative projector.|
|    `svd`   | `iters`    |         `200`      | Number of SVD iterations to perform. |
|    `svd`   | `resume`   |          `0`       | `0` to start a new SVD decomposition, `1` to resume from `svd_temp.dat`. |
|    `svd`   | `sstep`    |       `iters+1`    | Number of SVD iterations between printing partial errors and saving `svd_temp.dat` for resuming. |
|    `svd`   | `eigenmin` | `1.e-(0.9*DIGITS)` | Singular values `< eigenmin` are considered 0 for numerical stability. |
| `optional` |    `o`     |         `20`       | Oversampling factor for the dense grid. |
| `optional` |    `bw`    |          `4`       | Bandwidth reduction when adjusting a trigonometric polynomial. | 
| `optional` |    `z`     |         `C/2`      | Number of gridpoints for the smooth transition to zero. |
| `optional` |    `e`     |         `C`        | Number of gridpoints after the transition to zero region. |

## Output format
The resulting matrices are saved as double precision raw arrays in Fortran 
(i.e. column major) order and the endianness is fixed by the CPU architecture.

For the blend-to-zero operator the output filename is `AC-d`, where `C` is the 
number of continuation points and `d` the number of matching points. `Qd` is 
the output filename for the Dirichlet projector. 

In the case of projectors involving derivatives (`tkind = 1` or `tkind = 2`), 
the first element in the output file is the grid spacing and then the 
appropriate double precision projector array (`Qnd` or `Q2nd`, respectively). 
This grid spacing is useful for using the operators with different grid
spacings (via chain rule derivatives).

A sample Python script named `load_tables.py` showcases how to load each 
output in a Python environment. The included source code `fc_test.f08` 
shows how to load the operators in a Fortran environment.

## Verifying the tables
For convenience a simple Fortran 2008 program that compares numerical 
(using FC-Gram continuations) to analytical derivatives of the function 
`J_0(x)e^(3*x)` is also supplied. This program can verify Dirichlet operators
(that is, `Q` and `A`), and normal derivative operators (`Qn` and `Q2n`). Do
note, however, that checking for the derivative operators require the Dirichlet
operator for the same number of matching points, that is, to check `Qn7.dat`, 
`Q7.dat` is also required. If a Fortran 2008 compiler isn't available to 
you, it should be straight forward to edit `fc_test.f08`, edit the target 
function (and its derivative) and compile it as a Fortran 2003 program.

To compile the test program simply run, after compiling `fc_tables`, the 
command `make test` and execute the newly generated binary as `./fc_test N`,
where `N` the number of grid points to use (which must be a power of 2 due 
to the simplicity of the `fft` implementation).

The necessary table is automatically loaded and tested based on the values 
present in the `parameter.inp` file.

### Author
Mauro Fontana - Department of Physics - University of Buenos Aires.
