#include <time.h>
#include "adjoint.h"
#include "structures.h"
#include "affine.h"
#include "orthoderivative.h"

/**
 * Create a standard basis, {b_1, ..., b_n}
 * @param n The dimension
 * @return Standard basis
 */
size_t *createStandardBasis(size_t n);

/**
 * Add a constant c to a function F, s.t. F = F + c
 * @param F A function F to add a constant to
 * @param c The value of the constant
 */
void addConstant(TruthTable *F, size_t c);

/**
 * Check if a function F is affine
 * @param F The function F to check
 * @return True if the function F is affine, false otherwise
 */
bool isAffine(TruthTable *F);

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
        functionG = createTruthTable(functionF); // Create a random function G with respect to F
        printf("G:\n");
        printTruthTable(functionG);
    }
    TruthTable *orthoderivativeF = orthoderivative(functionF); // The orthoderivative of F
    TruthTable *orthoderivativeG = orthoderivative(functionG); // The orthoderivative of G

    Partition *partitionF = partitionTt(orthoderivativeF); // The partition of the orthoderivative of F
    basis = createStandardBasis(n); // Basis {b_1, ..., b_n}, here we use the standard basis.

    // Need to test for all possible constants, 0..2^n - 1.
    _Bool foundSolution = false; /* for breaking out of nested loops */
    for (size_t c1 = 0; c1 < 1L << n; ++c1) {
        TruthTable *ODGc = initTruthTable(n); // ODGc' = orthoderivativeG + c_1
        memcpy(ODGc->elements, orthoderivativeG->elements, sizeof(size_t) * 1L << n);
        addConstant(ODGc, c1); // Add the constant c1 to ODGc: ODGc' = ODGc + c_1
        Partition *partitionG = partitionTt(ODGc);
        size_t *mapOfPreImages = mapPreImages(partitionF, partitionG); // Create a mapping between the pre-images of F and ODGc

        // Calculate outer permutation, A1
        foundSolution = outerPermutation(partitionF, partitionG, n, basis, mapOfPreImages, orthoderivativeF, ODGc);

        destroyTruthTable(ODGc);
        destroyPartition(partitionG);
        free(mapOfPreImages);

        if (foundSolution) break;
    }

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
    printf("Affine\n");
    printf("Usage: affine [affine_options] [filenameF] [filenameG] \n");
    printf("Affine_options:\n");
    printf("\t-h \t- Print help\n");
    printf("\t-t \t- Print run time\n");
    printf("\n");
    printf("\tfilenameF = the path to file of function F\n");
    printf("\tfilenameG = the path to file of function G\n");
}

void addConstant(TruthTable *F, size_t c) {
    size_t dimension = F->n;
    for (size_t i = 0; i < 1L << dimension; ++i) {
        F->elements[i] ^= c;
    }
}

size_t *createStandardBasis(size_t n) {
    size_t *basis = malloc(sizeof(size_t) * n);
    for (size_t i = 0; i < n; ++i) {
        basis[i] = 1L << i;
    }
    return basis;
}

bool isAffine(TruthTable *F) {
    for (size_t a = 1; a < 1L << F->n; ++a) {
        for (size_t b = a + 1; b < 1L << F->n; ++b) {
            if (b > (a ^ b)) continue;
            size_t result = F->elements[0] ^ F->elements[a] ^ F->elements[b] ^ F->elements[a ^ b];
            if (result != 0) return false;
        }
    }
    return true;
}

