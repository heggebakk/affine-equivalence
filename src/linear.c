#include <time.h>
#include "structures.h"
#include "equivalence.h"
#include "orthoderivative.h"

/**
 * Print out a list over all the flags that can be used in the program
 */
void printLinearHelp();

int main(int argc, char *argv[]) {
    size_t n; // Working n
    size_t *basis; // List of the standard basis, {b_1, ..., b_n}
    RunTimes *runTime;
    bool times = false;
    clock_t startTotalTime;
    TruthTable *functionF = NULL;
    TruthTable *functionG = NULL;

    // Check for flags
    if (argc < 2) {
        printLinearHelp();
        return 0;
    }
    startTotalTime = clock();
    runTime = initRunTimes();

    // Loop over the arguments given
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] == '-') {
            switch (argv[i][1]) {
                case 'h':
                    printLinearHelp();
                    return 0;
                case 't':
                    times = true;
                    i++;
                    continue;
            }
        } else {
            if (functionF == NULL) {
                functionF = parseFile(argv[i]);
                printf("%s\n", argv[i]);
            } else if (functionG == NULL) {
                functionG = parseFile(argv[i]);
            }
        }
    }
    if (functionF == NULL) {
        printf("Missing function F. \n");
        return 0;
    }
    n = functionF->n;

    if (functionG == NULL) {
        functionG = createLinearFunction(functionF); // Create a random function G with respect to F
        printf("G:\n");
        printTruthTable(functionG);
    }

    Partition *partitionF = partitionTt(functionF); // The partition of the orthoderivative of F
    basis = createStandardBasis(n); // Basis {b_1, ..., b_n}, here we use the standard basis.

    // Need to test for all possible constants, 0..2^n - 1.
    _Bool foundSolution = false; /* for breaking out of nested loops */
    Partition *partitionG = partitionTt(functionG);
    size_t *mapOfPreImages = mapPreImages(partitionF, partitionG); // Create a mapping between the pre-images of F and ODGc

    // Calculate outer permutation, A1
    foundSolution = outerPermutation(partitionF, partitionG, n, basis, mapOfPreImages, functionF, functionG, false);

    destroyPartition(partitionG);
    free(mapOfPreImages);

    destroyTruthTable(functionF);
    destroyTruthTable(functionG);
    destroyPartition(partitionF);
    free(basis);
    runTime->total = stopTime(runTime->total, startTotalTime);
    if (times) {
        printTimes(runTime);
    }
    destroyRunTimes(runTime);
    return 0;
}

void printLinearHelp() {
    printf("Linear equivalence test\n");
    printf("Check for linear equivalences between two functions F and G.\n");
    printf("Usage: linear [linear_options] [filenameF] [filenameG] \n");
    printf("Linear_options:\n");
    printf("\t-h \t- Print help\n");
    printf("\t-t \t- Print run time\n");
    printf("\n");
    printf("\tfilenameF = the path to file of function F\n");
    printf("\tfilenameG = the path to file of function G\n");
}
