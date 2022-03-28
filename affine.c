#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include "structures.h"
#include <memory.h>
#include "affine.h"

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
    for (size_t i = 0; i < 1L << dimension; ++i) {
        fscanf(fp, "%zu", &tt->elements[i]);
    }
    fclose(fp);
    return tt;
}

BucketsMap *mapBuckets(struct Partition *f, struct Partition *g) {
    BucketsMap *bucketsMap = initBucketsMap();
    if (f->numBuckets != g->numBuckets) {
        return bucketsMap;
    }

    // Find all domains for the different buckets.
    struct Node **domains = malloc(sizeof(Node) * f->numBuckets);
    for (size_t i = 0; i < f->numBuckets; ++i) domains[i] = initNode();

    for (size_t i = 0; i < f->numBuckets; ++i) {
        bool sameSize = false;
        for (size_t j = 0; j < g->numBuckets; ++j) {
            if (f->bucketSizes[i] == g->bucketSizes[j]) {
                sameSize = true;
                addNode(domains[i], j);
            }
        }
        if (!sameSize) {
            return NULL;
        }
    }

    size_t numOfMappings = 1;
    bool *isCalculated = malloc(sizeof(bool) * f->numBuckets);
    memset(isCalculated, 0, sizeof(bool) * f->numBuckets);
    for (size_t i = 0; i < f->numBuckets; ++i) {
        Node *current = domains[i]->next;
        if (!isCalculated[current->data]) {
            isCalculated[current->data] = true;
            numOfMappings *= factorial(countNodes(domains[i]));
        }
    }
    free(isCalculated);

    bucketsMap->domains = malloc(sizeof(size_t *) * numOfMappings);
    createBucketsMap(bucketsMap, domains, f);

    // Free memory for domains
    for (size_t i = 0; i < f->numBuckets; ++i) {
        destroyNodes(domains[i]);
    }
    free(domains);

    return bucketsMap;
}

void createBucketsMap(BucketsMap *bucketsMap, Node **domains, Partition *partition) {
    bool *chosen = malloc(sizeof (bool) * partition->numBuckets);
    size_t *currentDomain = malloc(sizeof(size_t) * partition->numBuckets);
    memset(chosen, 0, sizeof(bool) * partition->numBuckets);
    selectRecursive(0, chosen, domains, partition, bucketsMap, currentDomain);

    free(chosen);
    free(currentDomain);
}

void selectRecursive(size_t i, bool *chosen, Node **domains, Partition *pG, BucketsMap *bucketsMap,
                     size_t *currentDomain) {
    if (i == pG->numBuckets) {
        size_t domainSize = pG->numBuckets;
        addDomain(bucketsMap, domainSize, currentDomain);
        return;
    }
    Node *current = domains[i]->next;
    while (current != NULL) {
        if (!chosen[current->data]) { // Check if we already have used the current data for the construction
            currentDomain[i] = current->data;
            chosen[current->data] = true;
            selectRecursive(i + 1, chosen, domains, pG, bucketsMap, currentDomain);
            chosen[current->data] = false;
        }
        current = current->next;
    }
}

void addDomain(BucketsMap *bucketsMap, size_t domainSize, size_t *domain) {
    size_t size = bucketsMap->numOfMappings;
    bucketsMap->domains[size] = malloc(sizeof(size_t *) * domainSize);
    memcpy(bucketsMap->domains[size], domain, sizeof(size_t) * domainSize);
    bucketsMap->numOfMappings += 1;
}

void calculateMultiplicities(TruthTable *truthTable, size_t *multiplicities) {
    size_t dimension = truthTable->dimension;
    for (size_t x = 0; x < 1L << dimension; ++x) {
        size_t y = truthTable->elements[x];
        multiplicities[y] += 1;
    }
}

size_t factorial(size_t value) {
    size_t factorial = 1;
    for (size_t i = 1; i < value + 1; ++i) {
        factorial *= i;
    }
    return factorial;
}

void add(TruthTable *dest, TruthTable *src) {
    for (size_t i = 0; i < 1L << dest->dimension; ++i) {
        dest->elements[i] ^= src->elements[i];
    }
}

TruthTable *compose(TruthTable *dest, TruthTable *src) {
    size_t dimension = dest->dimension;
    TruthTable *result = initTruthTable(dimension);
    for (size_t x = 0; x < 1L << dimension; ++x) {
        result->elements[x] = dest->elements[src->elements[x]];
    }
    return result;
}

TruthTable *randomLinearFunction(size_t dimension) {
    size_t entries = 1L << dimension;
    size_t basisImages[dimension];
    size_t listGenerated[entries];
    listGenerated[0] = 0;
    for (size_t i = 0; i < dimension; ++i) {
        size_t j = rand() % entries;
        basisImages[i] = j;
        for (size_t k = 0; k < 1L << i; ++k) {
            listGenerated[(1L << i) + k] = listGenerated[k] ^ j;
        }
    }
    TruthTable *result = initTruthTable(dimension);
    memcpy(result->elements, listGenerated, sizeof(size_t) * entries);
    return result;
}

TruthTable *randomLinearPermutation(size_t dimension) {
    size_t entries = 1L << dimension;
    size_t listGenerated[entries];
    size_t basisImages[dimension];
    bool *generated = malloc(sizeof(bool) * entries);
    memset(generated, 0, sizeof (bool) * entries);
    generated[0] = true;
    listGenerated[0] = 0;

    for (size_t i = 0; i < dimension; ++i) {
        size_t j = rand() % entries;
        while (generated[j]) {
            j = (j + 1) % entries;
        }
        basisImages[i] = j;
        for (size_t k = 0; k < 1L << i; ++k) {
            listGenerated[1L << i ^ k] = listGenerated[k] ^ j;
            generated[listGenerated[k] ^ j] = true;
        }
    }

    size_t c = rand() % entries;
    for (size_t i = 0; i < entries; ++i) {
        listGenerated[i] += (listGenerated[i] + c) % entries;
    }
    TruthTable *result = initTruthTable(dimension);
    memcpy(result->elements, listGenerated, sizeof(size_t) * entries);
    free(generated);
    return result;
}

TruthTable *createTruthTable(TruthTable *tt) {
    size_t dimension = tt->dimension;
    TruthTable *a1 = randomLinearPermutation(dimension);
    TruthTable *a2 = randomLinearPermutation(dimension);
    size_t c = rand() % 1L << dimension;
    for (size_t i = 0; i < 1L << dimension; ++i) {
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

void outerPermutation(Partition *f, Partition *g, size_t dimension, size_t *basis, TtNode *l1, size_t *domain) {
    size_t *images = malloc(sizeof(size_t) * dimension); // The images of the basis elements under l
    size_t *generated = malloc(sizeof(size_t) * 1L << dimension); // A partial truth table for l
    bool *generatedImages = malloc(sizeof(bool) * 1L << dimension);
    memset(generated, 0, sizeof(size_t) * 1L << dimension);
    memset(generatedImages, 0, sizeof(bool) * 1L << dimension);
    //generatedImages[0] = true;

    /**
     * Create dictionaries indexing buckets by elements
     * For instance, fClassPosition[i] would be the index of the bucket w.r.t. f containing the element i.
     */
    size_t *fClass = createClassRepresentation(f, dimension);
    size_t *gClass = createClassRepresentation(g, dimension);

    // Recursively guess the values of l on the basis (essentially, a dfs with backtracking upon contradiction
    guessValuesOfL(0, basis, images, f, g, dimension, generated, generatedImages, l1, fClass, gClass, domain);

    free(images);
    free(generated);
    free(generatedImages);
    free(fClass);
    free(gClass);
}

void
guessValuesOfL(size_t k, size_t *basis, size_t *images, Partition *f, Partition *g, size_t dimension, size_t *generated,
               bool *generatedImages, TtNode *l1, size_t *fClass, size_t *gClass, size_t *domain) {
    /**
     * If all basis elements have been assigned an image, and no contradictions have occurs, then we have found a
     * linear permutation preserving the partition. We reconstruct its truth table, and add it to the linked list
     * containing all permutations found by the search.
     */
    if (k == dimension) {
        TruthTable *new = initTruthTable(dimension);
        memcpy(new->elements, generated, sizeof(size_t) * 1L << dimension);
        addTtNode(l1, new);
        destroyTruthTable(new);
        return;
    }

    /**
     * We then take the bucket of the same size from the partition with respect to G. We know that the image of the
     * basis element must belong to that bucket.
     */
    size_t posBucketG = domain[fClass[basis[k]]];

    // We now go through all possible choices from the bucket
    for (size_t ick = 0; ick < g->bucketSizes[posBucketG]; ++ick) {
        size_t ck = g->buckets[posBucketG][ick];

        /**
         * Since we want the function to be a permutation, the image of the basis element should not be one of the
         * images that we already generated.
         */
        if (generatedImages[ck]) continue;

        /**
         * A contradiction can occur if assigning this value to the basis element causes some other element to map to
         * the wrong bucket by linearity. The following variable will be set to true if such a contradiction is
         * encountered.
         */
        bool problem = false;

        // This is done to handle the special case of k = 0, since otherwise we get (1L << (k - 1)) == MAX_INT.
        size_t LIMIT = k ? 1L << k : 1;

        /**
         * We now go trough all linear combinations of the basis elements that have been previously assigned.
         * Adding the newly guessed basis element to such combination allows us to derive one more value of the
         * function; if one of these values maps to the wrong bucket, we set "problem" = false, to indicate a
         * contradiction and backtrack.
         */
        for (size_t linearCombination = 0; linearCombination < LIMIT; ++linearCombination) {
            size_t x = linearCombination ^ basis[k];
            size_t y = ck;

            /**
             * The following loop simply XOR's all images corresponding to the linear combination, so that we get its
             * image over linearity.
             */
             if (k) {
                 for (size_t i = 0; i < k; ++i) {
                     if (1L << i & linearCombination) {
                         y ^= generated[basis[i]];
                     }
                 }
             }

             // Check for contradiction as described above
             if (f->bucketSizes[fClass[x]] != g->bucketSizes[gClass[y]]) {
                 problem = true;
		  //printf("Right now ck is %lu\n", ck);
		  //printf("Contradiction is due to %lu -> %lu\n",x ,y );
                 break;
             }

            // Add the new preimage-image pair to the partial truth table of the function
            generated[x] = y;

            // We also indicate that the image belongs to the set of generated images
            generatedImages[y] = true;
        }
        // If no contradiction is encountered, we go to the next basis element
        if (!problem) {
            images[k] = ck;
            guessValuesOfL(k + 1, basis, images, f, g, dimension, generated, generatedImages, l1, fClass, gClass,
                           domain);
        }

        // When backtracking, we need to reset the generated image indicators
        for (size_t linearCombinations = 0; linearCombinations < LIMIT; ++linearCombinations) {
            size_t y = ck;
            if (k) {
                for (size_t i = 0; i < k; ++i) {
                    if (1L << i & linearCombinations) {
                        y ^= generated[basis[i]];
                    }
                }
            }
            generatedImages[y] = false;
        }
    }
}

size_t *createClassRepresentation(Partition *partition, size_t dimension) {
    // Loop over each bucket and set the bucket pos for each value
    size_t *class = malloc(sizeof(size_t) * 1L << dimension);
    for (size_t i = 0; i < partition->numBuckets; ++i) {
        for (size_t j = 0; j < partition->bucketSizes[i]; ++j) {
            size_t value = partition->buckets[i][j];
            class[value] = i; // Set the bucket number in the position of the value
        }
    }
    return class;
}

TruthTable *inverse(TruthTable *truthTable) {
    size_t dimension = truthTable->dimension;
    TruthTable *result = initTruthTable(dimension);
    for (int x = 0; x < 1L << dimension; ++x) {
        size_t y = truthTable->elements[x];
        result->elements[y] = x;
    }
    return result;
}