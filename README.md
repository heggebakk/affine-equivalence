# Testing affine equivalence
## Build and run the project
### What you need
- [GCC](https://gcc.gnu.org/) to compile the program, version 10 or above.

### How to build and run the program
To build the program, make `compile.sh` executable and run `./compile.sh` at the command line in the root directory
of the project, which will generate an executable called `affine`.

The program uses flags, for help type `./affine -h` in the command line.

Flags:
```text
Affine
Usage: affine [affine_options] [filenameF] [filenameG] 
Affine_options:
	-h 	- Print help
	-t 	- Print run time

	filenameF = the path to file of function F
	filenameG = the path to file of function G

```
The `function F` is expected to be in the format:
```text
6
0 1 8 15 27 14 35 48 53 39 43 63 47 41 1 1 41 15 15 47 52 6 34 22 20 33 36 23 8 41 8 47 36 52 35 53 35 39 20 22 33 34 48 53 39 48 6 23 22 33 63 14 23 52 14 43 27 63 36 6 27 43 20 34 
```
This function is the Gold function of dimension 6, given in the representation as a truth table, where the first line is the dimension of the function, and the second line is all the elements of the function.

#### Examples
Example testing two functions `F` and `G`:
```shell
./affine -g path/to/functionG path/to/functionF
```

Example finding `A` for a function `F`:
```text
./affine -a path/to/functionF
```

### Results
The results are written to a file set by the flag `-w`(results.txt by default).

Suppose testing two functions `F` and `G`, and we get `A1 * ortho-derivative F * A2 = ortho-derivative G`. In that case, the program will give the two functions `A1` and `A2`.
If the file is empty, the program found no such relation.

If trying to find `A` where `A1 * F * A2 + A = G`,
the program will give the function `A`.
