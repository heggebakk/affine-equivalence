#include "structures.h"
#include "affine.h"

size_t *createBasis(size_t dimension);

void addConstant(TruthTable *tt, size_t c);

int main(int argc, char *argv[]) {
    char *filename = "resources/q_6_1.tt";
    size_t dimension;
    size_t *basis;
    TruthTable *functionF = parseFile(filename);
//    TruthTable *functionG = createTruthTable(functionF);
    TruthTable *functionG = parseFile("resources/g.tt");
    printTruthTable(functionF);
    printTruthTable(functionG);
    Partition *partitionF = partitionTt(functionF);
    dimension = functionF->dimension;
    basis = createBasis(dimension);

    /* Cheating */
    // Need to test for all possible constants, 0..2^n - 1.
    for (int c = 0; c < 1L << dimension; ++c) {
        TruthTable *gPrime = initTruthTable(dimension);
        memcpy(gPrime->elements, functionG->elements, sizeof(size_t) * 1L << dimension);
        addConstant(gPrime, c);
        Partition *partitionG = partitionTt(gPrime);
        BucketsMap *bucketsMap = mapBuckets(partitionF, partitionG);
        destroyTruthTable(gPrime);

        for (size_t map = 0; map < bucketsMap->numOfMappings; ++map) {
            // Calculate outer permutation
            TtNode *a1 = outerPermutation(partitionF, partitionG, dimension, basis, bucketsMap->domains[map]);
            size_t numPermutations = countTtNodes(a1);
            bool foundSolution = false;

            for (size_t i = 0; i < numPermutations; ++i) {
                TruthTable *a1Prime = getTtNode(a1, i);
                TruthTable *a1Inverse = inverse(a1Prime);
                TruthTable *gDPrime = compose(a1Inverse, functionG);
                TruthTable *aPrime;
                TruthTable *a2 = initTruthTable(dimension);

                if (innerPermutation(functionF, gDPrime, basis, a2, aPrime)) {
                    printf("Hello!\n");
                    exit(0);
                }

                destroyTruthTable(a1Inverse);
                destroyTruthTable(gDPrime);
                destroyTruthTable(a2);
            }
            destroyTtNode(a1);
        }
        destroyBucketsMap(bucketsMap);
        destroyPartition(partitionG);
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
        printf("%zu ", basis[i]);
    }
    printf("\n");
    return basis;
}
