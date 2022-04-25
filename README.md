# Testing affine equivalence

## Build and run the project
To build the program, run `compile.sh` at the commandline in the root directory
of the project, which will generate an executable called `affine.out`.

The program uses flags, for help type `affine -h` in the commanline.

Flags:
```text
Affine
Usage: affine [affine_options] [filename]
Affine_options:
        -h      - Print help
        -w      - The root filename where the results should be written to

        filename = the root filename of function F
        -h override all other options
```

Example:
```shell
./affine -w results.txt resources/dim6/q_6_1.tt
```

The `function F` is expected to be on the format:
```text
6
0 1 8 15 27 14 35 48 53 39 43 63 47 41 1 1 41 15 15 47 52 6 34 22 20 33 36 23 8 41 8 47 36 52 35 53 35 39 20 22 33 34 48 53 39 48 6 23 22 33 63 14 23 52 14 43 27 63 36 6 27 43 20 34 
```
Where the first line is the dimension of the function, and the second line is all the elements of the truth table.