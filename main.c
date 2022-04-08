#include "structures.h"
#include "affine.h"
#include "orthoderivative.h"

size_t *createBasis(size_t dimension);

void addConstant(TruthTable *tt, size_t c);

int main() {
//    char *filename = "resources/dim6/gf/orthoderivative_GF.tt";
    char *filename = "resources/dim6/gf/q_6_1.tt";
    size_t dimension;
    size_t *basis;
//    TruthTable *functionF = parseFile(filename);
    TruthTable *functionF= orthoderivative(parseFile(filename));
//    TruthTable *functionG = createTruthTable(functionF);
//    TruthTable *functionG = parseFile("resources/dim6/gf/orthoderivative_g.tt");
//    TruthTable *functionG = parseFile("resources/dim6/gf/g.tt");
    TruthTable *functionG = orthoderivative(parseFile("resources/dim6/gf/g.tt"));
//    TruthTable *functionG = orthoderivative(createTruthTable(functionF));
    printTruthTable(functionF);
    printTruthTable(functionG);
    Partition *partitionF = partitionTt(functionF);
    dimension = functionF->dimension;
    basis = createBasis(dimension);

    // Need to test for all possible constants, 0..2^n - 1.
    for (size_t c1 = 0; c1 < 1L << dimension; ++c1) {
        _Bool foundSolution = false; /* for breaking out of nested loops */
        TruthTable *gPrime = initTruthTable(dimension);
        memcpy(gPrime->elements, functionG->elements, sizeof(size_t) * 1L << dimension);
        addConstant(gPrime, c1); // Add the constant c1 to g: g' = g + c1
        Partition *partitionG = partitionTt(gPrime);
        BucketsMap *bucketsMap = mapBuckets(partitionF, partitionG);
        printf("Partition F:\n");
        printPartition(partitionF);
        printf("Partition G:\n");
        printPartition(partitionG);

        for (size_t map = 0; map < bucketsMap->numOfMappings; ++map) {
            // Calculate outer permutation
            TtNode *a1 = outerPermutation(partitionF, partitionG, dimension, basis, bucketsMap->domains[map]);
            size_t numPermutations = countTtNodes(a1);
            printf("Number of permutations: %zu\n", numPermutations);

            for (size_t i = 0; i < numPermutations; ++i) {
                TruthTable *a1Prime = getTtNode(a1, i);
                TruthTable *a1Inverse = inverse(a1Prime);
                TruthTable *gDoublePrime = compose(a1Inverse, gPrime);
                TruthTable *a2 = initTruthTable(dimension);
                a2->elements[0] = 0;

                if (innerPermutation(functionF, gDoublePrime, basis, a2)) {
                    printf("Constant c1: %zu\n", c1);
                    printf("affine function:\n");
                    printTruthTable(a2);
                    foundSolution = true;
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

    destroyTruthTable(functionF);
    destroyTruthTable(functionG);
    destroyPartition(partitionF);
    free(basis);

    return 0;
}

void addConstant(TruthTable *tt, size_t c) {
    size_t dimension = tt->dimension;
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
