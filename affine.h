#ifndef AFFINE_AFFINE_H
#define AFFINE_AFFINE_H

/**
 * Parse file containing a truth table
 * @param file The file path
 */
TruthTable * parseFile(char *file);

BucketsMap *mapBuckets(struct Partition *f, struct Partition *g);

void createBucketsMap(BucketsMap *bucketsMap, Node **domains, Partition *partition);

void selectRecursive(size_t i, bool *chosen, Node **domains, Partition *pG, BucketsMap *bucketsMap,
                     size_t *currentDomain);

void addDomain(BucketsMap *bucketsMap, size_t domainSize, size_t *domain);

size_t factorial(size_t value);

void calculateMultiplicities(TruthTable *truthTable, size_t *multiplicities);

TruthTable *randomLinearPermutation(size_t dimension);

TruthTable *randomLinearFunction(size_t dimension);

TruthTable * compose(TruthTable *dest, TruthTable *src);

void add(TruthTable *dest, TruthTable *src);

TruthTable *createTruthTable(TruthTable *tt);

void outerPermutation(Partition *f, Partition *g, size_t dimension, size_t *basis, TtNode *l1, size_t *domain);

void
guessValuesOfL(size_t k, size_t *basis, size_t *images, Partition *f, Partition *g, size_t dimension, size_t *generated,
               bool *generatedImages, TtNode *l1, size_t *fClass, size_t *gClass, size_t *domain);

size_t *createClassFromDomain(Partition *partition, size_t dimension);

#endif //AFFINE_AFFINE_H
