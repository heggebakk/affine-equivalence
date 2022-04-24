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

int main(int argc, char *argv[]) {
    char *filename; // Filepath for the function F
    char *writePath; // Path to file for writing the results
    size_t n; // Working n
    size_t *basis; // List of the standard basis, {b_1, ..., b_n}

    // TODO: Create a "interface" for the users
    if (argc < 2) {
        filename = "resources/dim6/q_6_1.tt"; // Defaults to GF(6) if file is not given
        writePath = "results.txt";
    } else {
        printf("%s\n", argv[1]);
        filename = argv[1]; // First parameter is the path to the truth table F
        writePath = argv[2]; // Second parameter is the path to write the results to
    }
    FILE *fp = fopen(writePath, "w+");
    fprintf(fp, "%s\n", filename); // Write the filename of the function F

    TruthTable *functionF1 = parseFile(filename); // Parsed truth table of function F
    TruthTable *functionG1 = createTruthTable(functionF1); // Create a random function G with respect to F
    TruthTable *orthoderivativeF = orthoderivative(functionF1);
    TruthTable *orthoderivativeG = orthoderivative(functionG1);

    printf("Function F:\n");
    printTruthTable(functionF1);
    printf("Function G:\n");
    printTruthTable(functionG1);
    printf("Orthoderivatives F & G:\n");
    printTruthTable(orthoderivativeF);
    printTruthTable(orthoderivativeG);

    Partition *partitionF = partitionTt(orthoderivativeF);
    n = orthoderivativeF->n;
    basis = createStandardBasis(n); // Basis {b_1, ..., b_n}

    // Need to test for all possible constants, 0..2^n - 1.
    for (size_t c1 = 0; c1 < 1L << n; ++c1) {
        _Bool foundSolution = false; /* for breaking out of nested loops */
        TruthTable *G = initTruthTable(n); // G' = orthoderivativeG + c_1
        memcpy(G->elements, orthoderivativeG->elements, sizeof(size_t) * 1L << n);
        addConstant(G, c1); // Add the constant c1 to G: G' = G + c_1
        Partition *partitionG = partitionTt(G);
        BucketsMap *bucketsMap = mapBuckets(partitionF, partitionG); // Map F -> G

        for (size_t map = 0; map < bucketsMap->numOfMappings; ++map) {
            // Calculate outer permutation, A1
            TtNode *A1 = outerPermutation(partitionF, partitionG, n, basis, bucketsMap->mappings[map]);
            size_t numPermutations = countTtNodes(A1);

            // Go over all the possible permutations A1, and try for inner permutation L2
            for (size_t i = 0; i < numPermutations; ++i) {
                TruthTable *currentA1 = getTtNode(A1, i); // A1[i], the A1 we are testing
                TruthTable *A1Inverse = inverse(currentA1); // A1^{-1}
                TruthTable *GPrime = compose(A1Inverse, G); // A1^{-1} * G = G'
                TruthTable *A2 = initTruthTable(n);
                A2->elements[0] = 0; // We know that the function is linear => L[0] -> 0

                if (innerPermutation(orthoderivativeF, GPrime, basis, A2, fp)) {
                    foundSolution = true;

                    printf("Constant c1: %zu\n", c1);
                    fprintf(fp, "Constant c1: %zu\n", c1);
                    printf("A1: \n");
                    fprintf(fp, "A1:\n");
                    printTruthTable(currentA1);
                    writeTruthTable(currentA1, fp);
                    printf("A2:\n");
                    fprintf(fp, "A2:\n");
                    printTruthTable(A2);
                    writeTruthTable(A2, fp);

                    /* At this point, we know (A1,A2) linear s.t. A1 * orthoderivativeF * A2 = orthoderivativeG
		            *
		            * If A1 * F * A2 + A = G for the actual functions F and G (as opposed to the ODs),
		            * then A1 = A1Inverse, and A2 = A2
		            */
                    TruthTable *fComposeA2 = compose(functionF1, A2); // F * A2
                    TruthTable *A = compose(A1Inverse, fComposeA2); // A1Inverse * F * A2
                    add(A, functionG1); // A1Inverse * F * A2 + G = A
                    printf("A:\n");
                    printTruthTable(A);
                    printf("A is affine %s\n", isAffine(A) ? "True" : "False");

                    destroyTruthTable(fComposeA2);
                    destroyTruthTable(A);
                }
                destroyTruthTable(A1Inverse);
                destroyTruthTable(GPrime);
                destroyTruthTable(A2);

                if (foundSolution) break;
            }
            destroyTtNode(A1);

            if (foundSolution) break;
        }
        destroyTruthTable(G);
        destroyBucketsMap(bucketsMap);
        destroyPartition(partitionG);

        if (foundSolution) break;
    }

    destroyTruthTable(functionF1);
    destroyTruthTable(functionG1);
    destroyTruthTable(orthoderivativeF);
    destroyTruthTable(orthoderivativeG);
    destroyPartition(partitionF);
    free(basis);
    fclose(fp);

    return 0;
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

