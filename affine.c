#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include "structures.h"
#include <memory.h>
#include "affine.h"

TruthTable *parseFile(char *file) {
    size_t n; // Dimension of the truth table
    FILE *fp = fopen(file, "r");

    // Check if file is found
    if (fp == NULL) {
        printf("File, %s, not found", file);
        fclose(fp);
        exit(1);
    }

    // If the file is found, we start parsing the file:
    fscanf(fp, "%zu", &n); // Set the n of the truth table
    TruthTable *f = initTruthTable(n);

    // Set all the elements in the truth table. Should be 2^n elements.
    for (size_t i = 0; i < 1L << n; ++i) {
        fscanf(fp, "%zu", &f->elements[i]);
    }

    fclose(fp);
    return f;
}

BucketsMap *mapBuckets(struct Partition *f, struct Partition *g) {
    BucketsMap *bucketsMap = initBucketsMap();
    // Check for contradictions between f and g
    if (f->numBuckets != g->numBuckets) {
        return bucketsMap;
    }

    // Find all domains for the different buckets.
    struct Node **domains = malloc(sizeof(Node) * f->numBuckets);
    for (size_t i = 0; i < f->numBuckets; ++i) domains[i] = initNode();

    for (size_t i = 0; i < f->numBuckets; ++i) {
        bool sameSize = false; // For checking contradiction.
        for (size_t j = 0; j < g->numBuckets; ++j) {
            if (f->bucketSizes[i] == g->bucketSizes[j]) {
                sameSize = true;
                addNode(domains[i], j);
            }
        }
        // If the partitions don't have the same sizes of the buckets, we have a contradiction
        if (!sameSize) return NULL;
    }

    size_t numOfMappings = 1; // If no contradictions, we know that there is at least one valid mapping of f and g.

    bool *isCalculated = calloc(sizeof(bool), f->numBuckets); // Boolean map for keeping track if we have used a bucket
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
    bool *chosen = calloc(sizeof (bool), partition->numBuckets);
    size_t *currentDomain = malloc(sizeof(size_t) * partition->numBuckets);
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

void calculateMultiplicities(TruthTable *f, size_t *multiplicities) {
    size_t dimension = f->dimension;
    for (size_t x = 0; x < 1L << dimension; ++x) {
        size_t y = f->elements[x];
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

TruthTable *compose(TruthTable *f, TruthTable *g) {
    size_t dimension = f->dimension;
    TruthTable *result = initTruthTable(dimension);
    for (size_t x = 0; x < 1L << dimension; ++x) {
        result->elements[x] = f->elements[g->elements[x]];
    }
    return result;
}

TruthTable *randomAffineFunction(size_t dimension) {
    size_t entries = 1L << dimension;
    size_t listGenerated[entries];
    listGenerated[0] = 0;
    size_t basisImages[dimension];
    for (size_t i = 0; i < dimension; ++i) {
        size_t j = rand() % entries;
        basisImages[i] = j;
        for (int k = 0; k < 1L << i; ++k) {
            listGenerated[(1L << i) + k] = listGenerated[k] ^ j;
        }
    }
    TruthTable *result = initTruthTable(dimension);
    memcpy(result->elements, listGenerated, sizeof(size_t) * entries);
    size_t randConstant = rand() % entries;
    for (int i = 0; i < entries; ++i) {
        result->elements[i] ^= randConstant;
    }
    return result;
}

TruthTable *randomAffinePermutation(size_t dimension) {
    size_t entries = 1L << dimension;
    bool generated[entries];
    size_t listGenerated[entries];
    generated[0] = true;
    for (size_t i = 1; i < entries; ++i) {
        generated[i] = false;
    }
    listGenerated[0] = 0;

    size_t basisImages[dimension];
    for (int i = 0; i < dimension; ++i) {
        size_t j = rand() % entries;
        while (generated[j]) {
            j = (j + 1) % entries;
        }
        basisImages[i] = j;
        for (int k = 0; k < 1L << i; ++k) {
            listGenerated[1L << i ^ k] = listGenerated[k] ^ j;
            generated[listGenerated[k] ^ j] = true;
        }
    }
    TruthTable *result = initTruthTable(dimension);
    memcpy(result->elements, listGenerated, sizeof(size_t) * entries);
    size_t randConstant = rand() % entries;
    for (int i = 0; i < entries; ++i) {
        result->elements[i] ^= randConstant;
    }
    return result;
}

TruthTable *createTruthTable(TruthTable *f) {
    size_t dimension = f->dimension;
    TruthTable *a1 = randomAffinePermutation(dimension);
    TruthTable *a2 = randomAffinePermutation(dimension);
    TruthTable *a = randomAffineFunction(dimension);

    TruthTable *temp = compose(f, a2);
    TruthTable *g = compose(a1, temp);
    add(g, a);

    destroyTruthTable(a1);
    destroyTruthTable(a2);
    destroyTruthTable(a);
    destroyTruthTable(temp);
    return g;
}

TtNode * outerPermutation(Partition *f, Partition *g, size_t dimension, size_t *basis, size_t *domain) {
    TtNode *a1 = initTtNode();
    size_t *images = malloc(sizeof(size_t) * dimension); // The images of the basis elements under l
    size_t *generated = calloc(sizeof(size_t), 1L << dimension); // A partial truth table for l
    bool *generatedImages = calloc(sizeof(bool), 1L << dimension);

    /**
     * Create dictionaries indexing buckets by elements
     * For instance, fClassPosition[i] would be the index of the bucket w.r.t. f containing the element i.
     */
    size_t *fClass = createClassRepresentation(f, dimension);
    size_t *gClass = createClassRepresentation(g, dimension);

    // Recursively guess the values of l on the basis (essentially, a dfs with backtracking upon contradiction
    guessValuesOfL(0, basis, images, f, g, dimension, generated, generatedImages, a1, fClass, gClass, domain);

    free(images);
    free(generated);
    free(generatedImages);
    free(fClass);
    free(gClass);
    return a1;
}

void
guessValuesOfL(size_t k, size_t *basis, size_t *images, Partition *f, Partition *g, size_t dimension, size_t *generated,
               bool *generatedImages, TtNode *a1, size_t *fClass, size_t *gClass, size_t *domain) {
    /**
     * If all basis elements have been assigned an image, and no contradictions have occurs, then we have found a
     * linear permutation preserving the partition. We reconstruct its truth table, and add it to the linked list
     * containing all permutations found by the search.
     */
    if (k == dimension) {
        TruthTable *new = initTruthTable(dimension);
        memcpy(new->elements, generated, sizeof(size_t) * 1L << dimension);
        addTtNode(a1, new);
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
            guessValuesOfL(k + 1, basis, images, f, g, dimension, generated, generatedImages, a1, fClass, gClass,
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

TruthTable *inverse(TruthTable *f) {
    size_t dimension = f->dimension;
    TruthTable *inverse = initTruthTable(dimension);
    for (size_t x = 0; x < 1L << dimension; ++x) {
        size_t y = f->elements[x];
        inverse->elements[y] = x;
    }
    return inverse;
}

bool *computeSetOfTs(TruthTable *f, const size_t x) {
    size_t dimension = f->dimension;
    bool *map = calloc(sizeof(bool), 1L << dimension);
    for (size_t y = 0; y < 1L << dimension; ++y) {
        size_t t = f->elements[x] ^ f->elements[y] ^ f->elements[x ^ y];
        map[t] = true;
    }
    return map;
}

Node *computeDomain(TruthTable *f, const bool *map) {
    size_t dimension = f->dimension;
    bool *domain = calloc(sizeof(bool), 1L << dimension);
    for (size_t i = 0; i < 1L << dimension; ++i) {
        domain[i] = true;
    }
    for (size_t t = 0; t < 1L << dimension; ++t) {
        if (map[t]) {
            bool *tempSet = calloc(sizeof(bool), 1L << dimension);
            for (size_t x = 0; x < 1L << dimension; ++x) {
                for (size_t y = 0; y < 1L << dimension; ++y) {
                    if (t == (f->elements[x] ^ f->elements[y] ^ f->elements[x ^ y])) {
                        tempSet[x] = true;
                        tempSet[y] = true;
                        tempSet[x ^ y] = true;
                    }
                }
            }
            for (size_t i = 0; i < 1L << dimension; ++i) {
                domain[i] &= tempSet[i];
            }
            free(tempSet);
        }
    }
    Node *domainResult = initNode();
    for (size_t i = 0; i < 1L << dimension; ++i) {
        if (domain[i]) {
            addNode(domainResult, i);
        }
    }
    free(domain);
    return domainResult;
}

bool innerPermutation(TruthTable *f, TruthTable *g, const size_t *basis, TruthTable *a2, FILE *fp) {
    size_t dimension = f->dimension;
    Node **restrictedDomains = malloc(sizeof(Node **) * (dimension + 1));
    bool result;

    for (size_t i = 0; i < dimension; ++i) {
        bool *map = computeSetOfTs(g, basis[i]);
        restrictedDomains[i] = computeDomain(f, map);
        free(map);
    }

    size_t *values = malloc(sizeof(size_t) * dimension);

    size_t constant_term = g->elements[0];
    /* Guess of constant term of a2 */

    for (size_t c2 = 0; c2 < 1L << dimension; ++c2) {
        /* Only consider preimages of g(0) */
        if (f->elements[c2] != constant_term) {
            continue;
        }

        TruthTable *newG = initTruthTable(dimension);
        for (size_t x = 0; x < 1L << dimension; ++x) {
            newG->elements[x ^ c2] = g->elements[x];
        }

        result = dfs(restrictedDomains, 0, values, f, newG, a2, basis);
        if (result) {
            /* If we get a result, we have to add the constant to the linear function that we found in dfs, and check if
             * f * l2 + c = g */
            if (c2 != 0) {
                for (int x = 0; x < 1L << dimension; ++x) {
                    a2->elements[x] ^= c2;
                }
            }

            /* If everything went smoothly, we should have aPrime == g */
            TruthTable *aPrime = compose(f, a2);
            printf("Constant c2: %zu\n", c2);
            fprintf(fp, "Constant c2: %zu\n", c2);

            destroyTruthTable(aPrime);
            destroyTruthTable(newG);
            break;
        }
        destroyTruthTable(newG);
    }

    free(values);

    for (size_t i = 0; i < dimension; ++i) {
        destroyNodes(restrictedDomains[i]);
    }
    free(restrictedDomains);
    return result;
}

bool dfs(Node **domains, size_t k, size_t *values, TruthTable *f, TruthTable *g, TruthTable *a2, const size_t *basis) {
    size_t dimension = f->dimension;
    if (k == dimension) return true;

    Node *current = domains[k]->next;
    while (current != NULL) {
        /* Guess that basis element #k maps to current->data */
        values[k] = current->data;
        /* Fill up part of the truth table (on the span of the guessed elements) */
        _Bool problem = false;
        for (size_t linear_combination = 0; linear_combination < (1L << k); ++linear_combination) {
            /* We are using the standard basis, and therefore the linear combination is the same
             * as the vector describing it
             */
            size_t new_input = linear_combination ^ (1L << k);
            size_t new_value = a2->elements[linear_combination] ^ current->data;
            a2->elements[new_input] = new_value;
            /* Check for a violation of f * a2 = g */
            if (f->elements[new_value] != g->elements[new_input]) {
                /* Something is wrong, backtrack */
                problem = true;
                break;
            }
        }
        if (!problem) {
            if (dfs(domains, k + 1, values, f, g, a2, basis)) return true;
        }
        current = current->next;
    }
    return false;
}
