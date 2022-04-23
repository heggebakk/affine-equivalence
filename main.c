#include "structures.h"
#include "affine.h"
#include "orthoderivative.h"

size_t *createBasis(size_t dimension);

void addConstant(TruthTable *tt, size_t c);

bool isAffine(TruthTable *f);

int main(int argc, char *argv[]) {
    char *filename;
    char *wp; // Path to file for writing the results
    size_t dimension; // Working dimension
    size_t *basis; // List of the standard basis
    if (argc < 2) {
        filename = "resources/dim6/q_6_1.tt"; // Default to GF(6)
        wp = "results.txt";
    } else {
        printf("%s\n", argv[1]);
        filename = argv[1]; // First parameter is the path to the truth table F
        wp = argv[2]; // Second parameter is the path to write the results to
    }
    FILE *fp = fopen(wp, "w+");
    fprintf(fp, "%s\n", filename);
    TruthTable *functionF1 = parseFile(filename); // Parsed truth table F
    TruthTable *functionG1 = createTruthTable(functionF1); // Affine function G
    printf("Function F:\n");
    printTruthTable(functionF1);
    printf("Function G:\n");
    printTruthTable(functionG1);
    TruthTable *orthoderivativeF = orthoderivative(functionF1);
    TruthTable *orthoderivativeG = orthoderivative(functionG1);
    printf("Orthoderivatives F & G:\n");
    printTruthTable(orthoderivativeF);
    printTruthTable(orthoderivativeG);
    TruthTable *functionF = orthoderivativeF;
    TruthTable *functionG = orthoderivativeG;
    Partition *partitionF = partitionTt(functionF);
    dimension = functionF->n;
    basis = createBasis(dimension);

    // Need to test for all possible constants, 0..2^n - 1.
    for (size_t c1 = 0; c1 < 1L << dimension; ++c1) {
        _Bool foundSolution = false; /* for breaking out of nested loops */
        TruthTable *gPrime = initTruthTable(dimension);
        memcpy(gPrime->elements, functionG->elements, sizeof(size_t) * 1L << dimension);
        addConstant(gPrime, c1); // Add the constant c1 to g: g' = g + c1
        Partition *partitionG = partitionTt(gPrime);
        BucketsMap *bucketsMap = mapBuckets(partitionF, partitionG);

        for (size_t map = 0; map < bucketsMap->numOfMappings; ++map) {
            // Calculate outer permutation
            TtNode *a1 = outerPermutation(partitionF, partitionG, dimension, basis, bucketsMap->mappings[map]);
            size_t numPermutations = countTtNodes(a1);

            for (size_t i = 0; i < numPermutations; ++i) {
                TruthTable *a1Prime = getTtNode(a1, i);
                TruthTable *a1Inverse = inverse(a1Prime);
                TruthTable *gDoublePrime = compose(a1Inverse, gPrime);
                TruthTable *a2 = initTruthTable(dimension);
                a2->elements[0] = 0;

                if (innerPermutation(functionF, gDoublePrime, basis, a2, fp)) {
                    printf("Constant c1: %zu\n", c1);
                    fprintf(fp, "Constant c1: %zu\n", c1);
                    printf("a1: \n");
                    fprintf(fp, "a1:\n");
                    printTruthTable(a1Prime);
                    writeTruthTable(a1Prime, fp);
                    printf("a2:\n");
                    fprintf(fp, "a2:\n");
                    printTruthTable(a2);
                    writeTruthTable(a2, fp);
                    foundSolution = true;

                    /* At this point, we know (A1,A2) linear s.t. A1 * functionF * A2 = gPrime, where
		            * functionF is the OD of F and gPrime is the OD of G
		            *
		            * If L1 * F * L2 + A = G for the actual functions F and G (as opposed to the ODs),
		            * then L1 = a1Inverse, and L2 = a2
		            */
                    TruthTable *fComposeA2 = compose(functionF1, a2); // F * A2
                    TruthTable *a = compose(a1Inverse, fComposeA2); // A1Inverse * F * A2
                    add(a, functionG1); // A1Inverse * F * A2 + G = A
                    destroyTruthTable(fComposeA2);
                    printf("A:\n");
                    printTruthTable(a);
                    printf("A is affine %s\n", isAffine(a) ? "True" : "False");
                    destroyTruthTable(a);
                }
                destroyTruthTable(a1Inverse);
                destroyTruthTable(gDoublePrime);
                destroyTruthTable(a2);

                if (foundSolution) break;
            }
            destroyTtNode(a1);

            if (foundSolution) break;
        }
        destroyTruthTable(gPrime);
        destroyBucketsMap(bucketsMap);
        destroyPartition(partitionG);

        if (foundSolution) break;
    }

    destroyTruthTable(functionF1);
    destroyTruthTable(functionG1);
    destroyTruthTable(functionF);
    destroyTruthTable(functionG);
    destroyPartition(partitionF);
    free(basis);
    fclose(fp);

    return 0;
}

void addConstant(TruthTable *tt, size_t c) {
    size_t dimension = tt->n;
    for (size_t i = 0; i < 1L << dimension; ++i) {
        tt->elements[i] ^= c;
    }
}

size_t *createBasis(size_t dimension) {
    size_t *basis = malloc(sizeof(size_t) * (dimension));
    for (size_t i = 0; i < dimension; ++i) {
        basis[i] = 1L << i;
    }
    return basis;
}

bool isAffine(TruthTable *f) {
    for (size_t a = 1; a < 1L << f->n; ++a) {
        for (size_t b = a + 1; b < 1L << f->n; ++b) {
            if (b > (a ^ b)) continue;
            size_t result = f->elements[0] ^ f->elements[a] ^ f->elements[b] ^ f->elements[a ^ b];
            if (result != 0) return false;
        }
    }
    return true;
}

