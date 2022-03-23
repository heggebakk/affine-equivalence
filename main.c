#include "structures.h"
#include "affine.h"

TruthTable *createTruthTable(TruthTable *tt);

int main(int argc, char *argv[]) {
    char *filename = "resources/q_6_1.tt";
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
    BucketsMap *bucketsMap = mapBuckets(partitionF, partitionG, functionF->dimension);

    destroyTruthTable(functionF);
    destroyPartition(partitionF);
    return 0;
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

