#include <time.h>
#include "adjoint.h"
#include "structures.h"
#include "affine.h"
#include "orthoderivative.h"

/**
 * Print out a list over all the flags that can be used in the program
 */
void printHelp();

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
        printHelp();
        return 0;
    }
    startTotalTime = clock();
    runTime = initRunTimes();

    // Loop over the arguments given
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] == '-') {
            switch (argv[i][1]) {
                case 'h':
                    printHelp();
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
    // L1 and L2 are the linear functions from the creation of a random G with respect to F
    n = functionF->n;

    if (functionG == NULL) {
        functionG = createLinearTruthTable(functionF); // Create a random function G with respect to F
        printf("G:\n");
        printTruthTable(functionG);
    }
    TruthTable *orthoderivativeF = orthoderivative(functionF); // The orthoderivative of F
    TruthTable *orthoderivativeG = orthoderivative(functionG); // The orthoderivative of G

    Partition *partitionF = partitionTt(orthoderivativeF); // The partition of the orthoderivative of F
    basis = createStandardBasis(n); // Basis {b_1, ..., b_n}, here we use the standard basis.

    // Need to test for all possible constants, 0..2^n - 1.
    Partition *partitionG = partitionTt(orthoderivativeG);
    size_t *mapOfPreImages = mapPreImages(partitionF, partitionG); // Create a mapping between the pre-images of F and ODGc

    // Calculate outer permutation, A1
    outerPermutation(partitionF, partitionG, n, basis, mapOfPreImages, orthoderivativeF, orthoderivativeG,
                     false);

    destroyPartition(partitionG);
    free(mapOfPreImages);


    destroyTruthTable(functionF);
    destroyTruthTable(functionG);
    destroyTruthTable(orthoderivativeF);
    destroyTruthTable(orthoderivativeG);
    destroyPartition(partitionF);
    free(basis);
    runTime->total = stopTime(runTime->total, startTotalTime);
    if (times) {
        printTimes(runTime);
    }
    destroyRunTimes(runTime);
    return 0;
}

void printHelp() {
    printf("Linear test\n");
    printf("Usage: linear [linear_options] [filenameF] [filenameG] \n");
    printf("Linear_options:\n");
    printf("\t-h \t- Print help\n");
    printf("\t-t \t- Print run time\n");
    printf("\n");
    printf("\tfilenameF = the path to file of function F\n");
    printf("\tfilenameG = the path to file of function G\n");
}
