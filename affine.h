#ifndef AFFINE_AFFINE_H
#define AFFINE_AFFINE_H

/**
 * Parse file containing the elements of a truth table. The first line is the dimension of the truth table. The
 * second line contains all the elements of the truth table:
 * For the GF(6):
 * 6
 * 0 1 8 15 27 14 35 48 53 39 43 63 47 41 1 1 41 15 15 47 52 6 34 22 20 33 36 23 8 41 8 47 36 52 35 53 35 39 20 22 33 34 48 53 39 48 6 23 22 33 63 14 23 52 14 43 27 63 36 6 27 43 20 34
 * @param file The file path of the truth table
 * @return The parsed file as the struct Truth Table holding the information about the dimension and all the elements
 */
TruthTable *parseFile(char *file);

/**
 * Create a mapping of all the buckets between partition f and g, where f -> g.
 * @param f A partition of a function f
 * @param g A partition of a function g
 * @return A struct BucketsMap that holds all the possible mappings between f and g
 */
BucketsMap *mapBuckets(struct Partition *f, struct Partition *g);

void createBucketsMap(BucketsMap *bucketsMap, Node **domains, Partition *partition);

void selectRecursive(size_t i, bool *chosen, Node **domains, Partition *pG, BucketsMap *bucketsMap,
                     size_t *currentDomain);

void addDomain(BucketsMap *bucketsMap, size_t domainSize, size_t *domain);

size_t factorial(size_t value);

void calculateMultiplicities(TruthTable *f, size_t *multiplicities);

TruthTable *randomAffinePermutation(size_t dimension);

TruthTable *randomAffineFunction(size_t dimension);

TruthTable * compose(TruthTable *f, TruthTable *g);

void add(TruthTable *dest, TruthTable *src);

TruthTable *createTruthTable(TruthTable *f);

TtNode * outerPermutation(Partition *f, Partition *g, size_t dimension, size_t *basis, size_t *domain);

void
guessValuesOfL(size_t k, size_t *basis, size_t *images, Partition *f, Partition *g, size_t dimension, size_t *generated,
               bool *generatedImages, TtNode *a1, size_t *fClass, size_t *gClass, size_t *domain);

size_t *createClassRepresentation(Partition *partition, size_t dimension);

TruthTable * inverse(TruthTable *f);

bool *computeSetOfTs(TruthTable *f, size_t x);

Node *computeDomain(TruthTable *f, const bool *map);

bool innerPermutation(TruthTable *f, TruthTable *g, const size_t *basis, TruthTable *a2, FILE *fp);

bool dfs(Node **domains, size_t k, size_t *values, TruthTable *f, TruthTable *g, TruthTable *a2, const size_t *basis);

#endif //AFFINE_AFFINE_H
