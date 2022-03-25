#include "structures.h"
#include "affine.h"

size_t *createBasis(size_t dimension);

int main(int argc, char *argv[]) {
    char *filename = "resources/q_6_1.tt";
    size_t dimension;
    size_t *basis;
    TruthTable *functionF = parseFile(filename);
//    TruthTable *functionG = createTruthTable(functionF);
    TruthTable *functionG = parseFile("resources/g.tt");
    printTruthTable(functionF);
    printTruthTable(functionG);

    /* Cheating */
    for(size_t i = 0; i < (1L << functionG->dimension); ++i) {
      functionG->elements[i] ^= 10;
    }
    Partition *partitionF = partitionTt(functionF);
    Partition *partitionG = partitionTt(functionG);
    printf("\n");
    printPartition(partitionF);
    printPartition(partitionG);
    dimension = functionF->dimension;
    basis = createBasis(dimension);
    BucketsMap *bucketsMap = mapBuckets(partitionF, partitionG);

    for (size_t map = 0; map < bucketsMap->numOfMappings; ++map) {
        TtNode *l1 = initTtNode();
        bool foundSolution = false;

        // Calculate outer permutation
        outerPermutation(partitionF, partitionG, functionF->dimension, basis, l1, bucketsMap->domains[map]);


        destroyTtNode(l1);
    }

    destroyTruthTable(functionF);
    destroyTruthTable(functionG);
    destroyPartition(partitionF);
    destroyPartition(partitionG);
    free(basis);
    destroyBucketsMap(bucketsMap);
    return 0;
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
