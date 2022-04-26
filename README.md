# Testing affine equivalence
## Build and run the project
### What you need
- [GCC](https://gcc.gnu.org/) to compile the program, version 10 or above.

### How to build and run the program
To build the program, make `compile.sh` executable and run `./compile.sh` at the commandline in the root directory
of the project, which will generate an executable called `affine`.

The program uses flags, for help type `./affine -h` in the commandline.

Flags:
```text
Affine
Usage: ./affine [affine_options] [filename]
Affine_options:
	-a 	- Set this if you want to find the affine function A.
	-g 	- The root filename where the function G is found.	If not given, the program will compute a random G with respect to F.
	-h 	- Print help
	-w 	- The root filename where the results should be written to, default = "results.txt"

	filename = the root filename of function F
	-h override all other options
```
The `function F` is expected to be on the format:
```text
6
0 1 8 15 27 14 35 48 53 39 43 63 47 41 1 1 41 15 15 47 52 6 34 22 20 33 36 23 8 41 8 47 36 52 35 53 35 39 20 22 33 34 48 53 39 48 6 23 22 33 63 14 23 52 14 43 27 63 36 6 27 43 20 34 
```
This function is the `GF(6)` given in the representation as a truth table, where the first line is the dimension of the function, and the second line is all the elements of the function.

#### Examples
Example:
```shell
./affine -w test.txt path/to/functionF
```

Example testing two functions `F` and `G` 
```shell
./affine -g path/to/functionG path/to/functionF
```

Example finding `A` for a function `F`
```text
./affine -a path/to/functionF
```

### Results
The result is written to a file set by the flag `-w`(results.txt by default).

If testing two functions `F` and `G`, and we get `A1 * orthoderivative F * A2 = orthoderivative G` the program will give the two functions `A1` and `A2`.
If the file is empty, no such relation was found.

If trying to find `A` where `A1 * F * A2 + A = G`,
the program will give the function `A`
