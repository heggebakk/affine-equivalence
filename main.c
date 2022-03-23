#include "structures.h"
#include "affine.h"

TruthTable *createTruthTable(TruthTable *tt);

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
    Partition *partitionF = partitionTt(functionF);
    Partition *partitionG = partitionTt(functionG);
    printf("\n");
    printPartition(partitionF);
    printPartition(partitionG);
    dimension = functionF->dimension;
    basis = createBasis(dimension);
    BucketsMap *bucketsMap = mapBuckets(partitionF, partitionG, functionF->dimension);

    for (int map = 0; map < bucketsMap->numOfMappings; ++map) {
        TtNode *l1 = initTtNode();
        bool foundSolution = false;

        // Calculate outer permutation
//        outerPermutation(partitionF, partitionG, functionF->dimension, basis, l1, bucketsMap);
    }

    destroyTruthTable(functionF);
    destroyTruthTable(functionG);
    destroyPartition(partitionF);
    destroyPartition(partitionG);
    destroyBucketsMap(bucketsMap);
    return 0;
}

size_t *createBasis(size_t dimension) {
    size_t *basis = malloc(sizeof(size_t) * dimension + 1);
    basis[0] = 0;
    printf("Basis: %zu ", basis[0]);
    for (int i = 1; i < dimension + 1; ++i) {
        basis[i] = 1L << (i - 1);
        printf("%zu ", basis[i]);
    }
    printf("\n");
}

TruthTable *createTruthTable(TruthTable *tt) {
    size_t dimension = tt->dimension;
    TruthTable *a1 = randomLinearPermutation(dimension);
    TruthTable *a2 = randomLinearPermutation(dimension);
    size_t c = rand() % 1L << dimension;
    for (int i = 0; i < 1L << dimension; ++i) {
        a1->elements[i] = (a1->elements[i] + c) % (1L << dimension);
        a2->elements[i] = (a2->elements[i] + c) % (1L << dimension);
    }
    printTruthTable(a1);
    printTruthTable(a2);

    TruthTable *temp = compose(tt, a2);
    TruthTable *g = compose(a1, temp);

    destroyTruthTable(a1);
    destroyTruthTable(a2);
    destroyTruthTable(temp);
    return g;
}

