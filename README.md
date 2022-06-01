# Testing affine equivalence
## Build and run the project
### What you need
- [GCC](https://gcc.gnu.org/) to compile the program, version 10 or above.

### How to build and run the program
To build the program, make `compile.sh` executable and run `./compile.sh` at the command line in the root directory
of the project, which will generate three executables called  `ea_orthoderivative`, `affine` and `linear`.

The program uses flags, for help type `./ea_orthoderivative -h` in the command line.

Flags:
```text
> ./ea_orthoderivative -h
EA-equivalence
Check for EA-equivalence via their orthoderivatives.
Usage: ea [ea_options] [filenameF] [filenameG] 
Ea_options:
	-h 	- Print help
	-t 	- Print run time

	filenameF = the path to file of function F
	filenameG = the path to file of function G
```
```text
> ./affine -h
Affine test
Usage: affine [affine_options] [filenameF] [filenameG] 
Affine_options:
	-h 	- Print help
	-t 	- Print run time

	filenameF = the path to file of function F
	filenameG = the path to file of function G
```
```text
> ./linear -h
Linear equivalence test
Check for linear equivalences between two functions F and G.
Usage: linear [linear_options] [filenameF] [filenameG] 
Linear_options:
        -h      - Print help
        -t      - Print run time

        filenameF = the path to file of function F
        filenameG = the path to file of function G
```
The functions `F` and `G` is expected to be in the format:
```text
6
0 1 8 15 27 14 35 48 53 39 43 63 47 41 1 1 41 15 15 47 52 6 34 22 20 33 36 23 8 41 8 47 36 52 35 53 35 39 20 22 33 34 48 53 39 48 6 23 22 33 63 14 23 52 14 43 27 63 36 6 27 43 20 34 
```
This function is the Gold function of dimension 6, given in the representation as a truth table, where the first line is the dimension of the function, and the second line is all the elements of the function.

#### Examples
Example testing two functions `F` and `G` that are EA-equivalent:
```text
./ea_orthoderivative path/to/functionF path/to/functionG
L1:
0 1 12 13 36 37 40 41 25 24 21 20 61 60 49 48 28 29 16 17 56 57 52 53 5 4 9 8 33 32 45 44 50 51 62 63 22 23 26 27 43 42 39 38 15 14 3 2 46 47 34 35 10 11 6 7 55 54 59 58 19 18 31 30

L2:
0 49 36 21 17 32 53 4 47 30 11 58 62 15 26 43 54 7 18 35 39 22 3 50 25 40 61 12 8 57 44 29 28 45 56 9 13 60 41 24 51 2 23 38 34 19 6 55 42 27 14 63 59 10 31 46 5 52 33 16 20 37 48 1
```

## What the programs do
- `ea_orthoderivative`: Test for EA-equivalence between two function `F` and `G`;
- `affine`: Test for affine equivalence between two functions `F` and `G`;
- `linear`: Test for linear equivalence between two functions `F` and `G`.