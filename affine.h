#ifndef AFFINE_AFFINE_H
#define AFFINE_AFFINE_H

/**
 * Parse file containing a truth table
 * @param file The file path
 */
struct TruthTable * parseFile(char *file);

struct BucketsMap *mapBuckets(struct Partition *f, struct Partition *g, size_t dimension);

void createBucketsMap(BucketsMap *bucketsMap, struct Node **domains, struct Partition *pG);

void selectRecursive(int i, bool *chosen, struct Node **domains, struct Partition *pG, BucketsMap *bucketsMap,
                     size_t *currentDomain);

void addDomain(BucketsMap *bucketsMap, size_t domainSize, size_t *domain);

size_t factorial(size_t value);

void calculateMultiplicities(TruthTable *truthTable, size_t *multiplicities);

TruthTable *randomLinearPermutation(size_t dimension);

TruthTable *randomLinearFunction(size_t dimension);

TruthTable * compose(TruthTable *dest, TruthTable *src);

void add(TruthTable *dest, TruthTable *src);

#endif //AFFINE_AFFINE_H
