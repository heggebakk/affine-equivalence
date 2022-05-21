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
    char *pathF; // Filepath for the function F
    char *pathG = NULL; // Filepath for hte function G
    char *writePath = NULL; // Path to file for writing the results
    size_t n; // Working n
    size_t *basis; // List of the standard basis, {b_1, ..., b_n}
    bool computeAffineA = false; // Set to true if the user want to find A

    // Check for flags
    if (argc < 2) {
        printHelp();
        return 0;
    }
    // Loop over the arguments given
    pathF = argv[argc - 1];
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] == '-') {
            switch (argv[i][1]) {
                case 'a':
                    computeAffineA = true;
                    continue;
                case 'h':
                    printHelp();
                    return 0;
                case 'g':
                    i++;
                    pathG = argv[i];
                    continue;
                case 'w':
                    i++;
                    writePath = argv[i];
                    continue;
            }
        }
    }
    if (writePath == NULL) {
        // If the user have not sent in a pathF to write to, we create a default file
        writePath = "results.txt";
    }
    FILE *fp = fopen(writePath, "w+");
    fprintf(fp, "%s\n", pathF); // Write the pathF of the function F
    // L1 and L2 are the linear functions from the creation of a random G with respect to F
    TruthTable *functionF = parseFile(pathF); // Parsed truth table of function F
    TruthTable *functionG;
    n = functionF->n;
    if (pathG) {
        functionG = parseFile(pathG);
    } else {
        functionG = createTruthTable(functionF); // Create a random function G with respect to F
    }
    TruthTable *orthoderivativeF = orthoderivative(functionF); // The orthoderivative of F
    TruthTable *orthoderivativeG = orthoderivative(functionG); // The orthoderivative of G

    Partition *partitionF = partitionTt(orthoderivativeF); // The partition of the orthoderivative of F
    basis = createStandardBasis(n); // Basis {b_1, ..., b_n}, here we use the standard basis.

    // Need to test for all possible constants, 0..2^n - 1.
    _Bool foundSolution = false; /* for breaking out of nested loops */
    for (size_t c1 = 0; c1 < 1L << n; ++c1) {
        TruthTable *G = initTruthTable(n); // G' = orthoderivativeG + c_1
        memcpy(G->elements, orthoderivativeG->elements, sizeof(size_t) * 1L << n);
        addConstant(G, c1); // Add the constant c1 to G: G' = G + c_1
        Partition *partitionG = partitionTt(G);
        size_t *mapOfPreImages = mapPreImages(partitionF, partitionG); // Create a mapping between the pre-images of F and G

        // Calculate outer permutation, A1
        TtNode *A1 = outerPermutation(partitionF, partitionG, n, basis, mapOfPreImages);
        size_t numPermutations = countTtNodes(A1);

        // Go over all the possible permutations A1, and try for inner permutation L2
        for (size_t i = 0; i < numPermutations; ++i) {
            TruthTable *currentA1 = getTtNode(A1, i); // A1[i], the A1 we are testing
            TruthTable *A1Inverse = inverse(currentA1); // A1^{-1}
            TruthTable *GPrime = compose(A1Inverse, G); // A1^{-1} * G = G'
            TruthTable *A2 = initTruthTable(n);
            A2->elements[0] = 0; // We know that the function is linear => L[0] -> 0

            if (innerPermutation(orthoderivativeF, GPrime, basis, A2)) {
                /* At this point, we know (A1,A2) linear s.t. A1 * orthoderivativeF * A2 = orthoderivativeG */
                // We don't want to print out all the A1, A2 if the user looks for A
                if (!computeAffineA) {
                    foundSolution = true;
                    fprintf(fp, "A1:\n");
                    writeTruthTable(fp, currentA1);
                    fprintf(fp, "A2:\n");
                    writeTruthTable(fp, A2);
                }

                    /* Now, since the affine function A takes more time to compute, the user can use the flag -a to
                     * choose to compute A. */
                else {
                    /* If L1 * F * L2 + A = G for the actual functions F and G (as opposed to the ODs),
                    * then A1 = A1Inverse, and A2 = A2 */
                    TruthTable * A1Adjoint = adjoint(currentA1); // The adjoint of A1
                    TruthTable *A1AdjointInverse = inverse(A1Adjoint); // The inverse of the adjoint of A1 = L1
                    destroyTruthTable(A1Adjoint);

                    // Now compute L1 * F * L2 + G = A
                    TruthTable *fComposeL2 = compose(functionF, A2); // F * L2
                    TruthTable *A = compose(A1AdjointInverse, fComposeL2); // L*Inverse * F * L2
                    add(A, functionG); // L*Inverse * F * A2 + G = A
                    // Check if A is affine, if true, write result and quit
                    if (isAffine(A)) {
                        foundSolution = true;
                        fprintf(fp, "A:\n");
                        writeTruthTable(fp, A);
                    }

                    destroyTruthTable(fComposeL2);
                    destroyTruthTable(A);
                    destroyTruthTable(A1AdjointInverse);
                }
            }
            destroyTruthTable(A1Inverse);
            destroyTruthTable(GPrime);
            destroyTruthTable(A2);

            if (foundSolution) break;
        }
        destroyTtNode(A1);

        destroyTruthTable(G);
        destroyPartition(partitionG);
        free(mapOfPreImages);

        if (foundSolution) break;
    }

//    printf("Results found in \"%s\"\n", writePath);
    destroyTruthTable(functionF);
    destroyTruthTable(functionG);
    destroyTruthTable(orthoderivativeF);
    destroyTruthTable(orthoderivativeG);
    destroyPartition(partitionF);
    free(basis);
    fclose(fp);

    return 0;
}

void printHelp() {
    printf("Affine\n");
    printf("Usage: affine [affine_options] [filename]\n");
    printf("Affine_options:\n");
    printf("\t-a \t- Set this if you want to find the affine function A.\n");
    printf("\t-g \t- The root filename where the function G is found.\tIf not given, the program will compute a random G with respect to F.\n");
    printf("\t-h \t- Print help\n");
    printf("\t-w \t- The root filename where the results should be written to, default: \"results.txt\"\n");
    printf("\n");
    printf("\tfilename = the root filename of function F\n");
    printf("\t-h override all other options\n");
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

