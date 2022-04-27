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
BucketsMap *mapBuckets(Partition *F, Partition *G);

/**
 * A recursive function that creates all the possible mappings between F and G.
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

/**
 * A product of all positive descending integers. Factorial of n is denoted by n!.
 * @param n
 * @return n!
 */
size_t factorial(size_t n);

/**
 * Count all the occurrences of the elements in F and store them in a list, occurrences.
 * @param F A function F holding the elements to count
 * @param occurrences A list that holds the information about the occurrences of each element
 */
void countElements(TruthTable *F, size_t *occurrences);

/**
 * Reconstructing all linear permutations L1, respecting the partitions induced by function F and G
 * @param F Partition of function F
 * @param G Partition of function G
 * @param n Dimension
 * @param basis A basis = {b_1, ..., b_n}
 * @param map Tells how F -> G
 * @return All linear permutations L1
 */
TtNode *outerPermutation(Partition *F, Partition *G, size_t n, size_t *basis, size_t *map);

/**
 * Recursive function for reconstruction of all linear permutations L1
 * @param k Recursive step
 * @param basis A basis {b_1,...,b_n}
 * @param images Images of the basis elements under L
 * @param F Partition of function F
 * @param G Partition of function G
 * @param n Dimension
 * @param generated A partial truth table for L
 * @param generatedImages List, same size as images, holds the information if the images has been generated or not
 * @param L1 List that should hold all the linear permutations L1
 * @param fBucket Map of the buckets of function F
 * @param gBucket Map of the buckets of function G
 * @param map Tells how F -> G
 */
void
guessValuesOfL(size_t k, size_t *basis, size_t *images, Partition *F, Partition *G, size_t n, size_t *generated,
               bool *generatedImages, TtNode *L1, size_t *fBucket, size_t *gBucket, size_t *map);

/**
 * Create a list that tells in which bucket each element belongs to.
 * @param F Partition of a function F
 * @param n Dimension
 * @return A list holding the information of which bucket each element belongs to.
 */
size_t *createBucketRepresentation(Partition *F, size_t n);

/**
 * Compute the set of t's where t = F[x] + F[y] + F[x + y]
 * @param F Function containing the elements to compute the t's over
 * @param x A fixed value to compute with
 * @return A set of the T's from the computation
 */
bool *computeSetOfTs(TruthTable *F, size_t x);

/**
 * Compute the restricted domain for the given list of T's
 * @param F Function F
 * @param map A set of T's that we want to compute the restricted domain over
 * @return The restricted domain represented as a linked list
 */
Node *computeRestrictedDomains(TruthTable *F, const bool *map);

/**
 * Reconstruction of the inner permutation L2
 * @param F The truth table of function F
 * @param G The truth table of function G
 * @param basis A basis {b_1, ..., b_n}
 * @param L2 The inner permutation
 * @return Returns True if reconstruction of L2 was successful, False otherwise
 */
bool innerPermutation(TruthTable *F, TruthTable *G, const size_t *basis, TruthTable *L2);

/**
 * A dept first search to reconstruct the inner permutation L2.
 * @param domains The restricted domains over the list of T's computed previously
 * @param k The recursive step
 * @param values Initially a empty list, in use for guessing basis elements that maps to elements
 * @param F The truth table of a function F
 * @param G The truth table of a function G
 * @param L2 The innerpermutation
 * @param basis A basis {b_1,...,b_n}
 * @return
 */
bool dfs(Node **domains, size_t k, size_t *values, TruthTable *F, TruthTable *G, TruthTable *L2, const size_t *basis);

#endif //AFFINE_AFFINE_H
