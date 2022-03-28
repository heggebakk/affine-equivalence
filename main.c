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
        printf("\n");
        BucketsMap *bucketsMap = mapBuckets(partitionF, partitionG);
        destroyTruthTable(gPrime);

        printf("Num of domains: %zu\n", bucketsMap->numOfMappings);

        for (size_t map = 0; map < bucketsMap->numOfMappings; ++map) {
            printPartition(partitionF);
            printPartition(partitionG);
            TtNode *l1 = initTtNode();
            bool foundSolution = false;

            // Calculate outer permutation
            outerPermutation(partitionF, partitionG, dimension, basis, l1, bucketsMap->domains[map]);
            size_t numPermutations = countTtNodes(l1);
            printf("Number of permutations = %zu\n", numPermutations);

            for (size_t i = 0; i < numPermutations; ++i) {
            }
            destroyTtNode(l1);
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
