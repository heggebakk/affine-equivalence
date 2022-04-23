#ifndef AFFINE_AFFINE_H
#define AFFINE_AFFINE_H

/**
 * Parse file containing the elements of a truth table. The first line is the n of the truth table. The
 * second line contains all the elements of the truth table:
 * For the GF(6):
 * 6
 * 0 1 8 15 27 14 35 48 53 39 43 63 47 41 1 1 41 15 15 47 52 6 34 22 20 33 36 23 8 41 8 47 36 52 35 53 35 39 20 22 33 34 48 53 39 48 6 23 22 33 63 14 23 52 14 43 27 63 36 6 27 43 20 34
 * @param file The file path of the truth table
 * @return The parsed file as the struct Truth Table holding the information about the n and all the elements
 */
TruthTable *parseFile(char *file);

/**
 * Create a mapping of all the buckets between partition F and G, where F -> G.
 * If the partitions have several buckets of the same size, there is several possible mappings
 * Ex. F has the buckets with sizes [1,1,20,42] and G has [1,42,20,1] then we have the mappings
 * F[0] -> G[0,3], F[1] -> G[0,3], F[2] -> G[2], F[3] -> G[1] s.t we get the bucketsMap: [0, 3, 2, 1] and [3, 0, 2, 1].
 * @param F A partition of a function F
 * @param G A partition of a function G
 * @return A struct BucketsMap that holds all the possible mappings between F and G
 */
BucketsMap *mapBuckets(struct Partition *F, struct Partition *G);

/**
 * A recursive function that creates all the possible mappings between F and G
 * @param i Recursive step
 * @param chosen Keep track of the usage of the buckets for the construction of the maps
 * @param maps A list of all the mappings
 * @param partition The Partition of function F
 * @param bucketsMap The BucketsMap to store the mappings to
 * @param currentBucket A list of buckets, that keeps track if a bucket is already in use in the current recursive step.
 */
void mapBucketsRecursively(size_t i, bool *chosen, Node **maps, Partition *partition, BucketsMap *bucketsMap,
                           size_t *currentBucket);

/**
 * Add a new mapping to the BucketsMap
 * @param bucketsMap The BucketsMap to add to
 * @param numBuckets Number of buckets in the partition
 * @param map The map to add
 */
void addMapping(BucketsMap *bucketsMap, size_t numBuckets, size_t *map);

size_t factorial(size_t value);

void countMultiplicities(TruthTable *f, size_t *multiplicities);

TtNode * outerPermutation(Partition *f, Partition *g, size_t dimension, size_t *basis, size_t *domain);

void
guessValuesOfL(size_t k, size_t *basis, size_t *images, Partition *f, Partition *g, size_t dimension, size_t *generated,
               bool *generatedImages, TtNode *a1, size_t *fClass, size_t *gClass, size_t *domain);

size_t *createClassRepresentation(Partition *partition, size_t dimension);

bool *computeSetOfTs(TruthTable *f, size_t x);

Node *computeDomain(TruthTable *f, const bool *map);

bool innerPermutation(TruthTable *f, TruthTable *g, const size_t *basis, TruthTable *a2, FILE *fp);

bool dfs(Node **domains, size_t k, size_t *values, TruthTable *f, TruthTable *g, TruthTable *a2, const size_t *basis);

#endif //AFFINE_AFFINE_H
