#ifndef AFFINE_AFFINE_H
#define AFFINE_AFFINE_H

/**
 * Parse file containing a truth table
 * @param file The file path
 */
struct TruthTable * parseFile(char *file);

struct BucketsMap *mapBuckets(struct Partition *f, struct Partition *g, size_t dimension);

size_t factorial(size_t value);

void calculateMultiplicities(TruthTable *truthTable, size_t *multiplicities);

#endif //AFFINE_AFFINE_H
