#include "structures.h"
#include "affine.h"

} Partition;

/**
 * Initialize a new TruthTable
 * @param dimension The dimension of the TruthTable
 * @return The new TruthTable
 */
TruthTable *initTruthTable(size_t dimension);

/**
 * Parse file containing a truth table
 * @param file The file path
 */
TruthTable * parseFile(char *file);

/**
 * Free memory for struct TruthTable
 * @param truthTable The truth table to be freed
 */
void destroyTruthTable(TruthTable *truthTable);

Partition *partitionTt(TruthTable *truthTable);

void calculateMultiplicities(size_t i, size_t *multiplicities, TruthTable *truthTable, size_t x, size_t value);

Partition *initPartition(size_t dimension);

int main(int argc, char *argv[]) {
    char *filename = "resources/q_6_3.tt";
    TruthTable *functionF = parseFile(filename);
    Partition *partition = partitionTt(functionF);
    BucketsMap *bucketsMap = mapBuckets(partition, partition, functionF->dimension);

    destroyTruthTable(functionF);
    destroyPartition(partition);
    return 0;
}

Partition *partitionTt(TruthTable *truthTable) {
    Partition *partition;
    size_t dimension = truthTable->dimension;
    size_t *multiplicities = malloc(sizeof(size_t) * 1L << dimension);
    memset(multiplicities, 0, sizeof (size_t) * 1L << dimension);
    calculateMultiplicities(0, multiplicities, truthTable, 0, 0);
    partition = initPartition(1L << dimension);
}

Partition *initPartition(size_t dimension) {
    Partition *partition = malloc(sizeof(Partition));
    partition->multiplicities = malloc(sizeof(size_t) * dimension);
    partition->bucketSizes = malloc(sizeof(size_t) * dimension);
    partition->buckets = malloc(sizeof(size_t **) * dimension);
}

void calculateMultiplicities(size_t i, size_t *multiplicities, TruthTable *truthTable, size_t x, size_t value) {
    size_t dimension = truthTable->dimension;
    if (i == k - 1) {
        size_t numOfMultiplicities = value ^ truthTable->elements[x];
        multiplicities[numOfMultiplicities] += 1;
        return;
    } else {
        for (int y = 0; y < 1L << dimension; ++y) {
            size_t newX = x ^ y;
            size_t newMultiplicity = value ^ truthTable->elements[y];
            calculateMultiplicities(i + 1, multiplicities, truthTable, newX, newMultiplicity);
        }
    }
}

TruthTable *parseFile(char *file) {
    size_t dimension;
    FILE *fp = fopen(file, "r");

    // Check if file is found
    if (fp == NULL) {
        printf("File, %s, not found", file);
        fclose(fp);
        exit(1);
    }

    fscanf(fp, "%zu", &dimension);
    TruthTable *tt = initTruthTable(dimension);
    for (int i = 0; i < 1L << dimension; ++i) {
        fscanf(fp, "%zu", &tt->elements[i]);
    }
    fclose(fp);
    return tt;

}

void destroyTruthTable(TruthTable *truthTable) {
    free(truthTable->elements);
    free(truthTable);
}

TruthTable *initTruthTable(size_t dimension) {
    TruthTable *newTt = malloc(sizeof(TruthTable));
    newTt->dimension = dimension;
    newTt->elements = malloc(sizeof(size_t) * 1L << dimension);
    return newTt;
}
