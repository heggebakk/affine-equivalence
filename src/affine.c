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
        printf("File, %s, not found\n", file);
        fclose(fp);
        exit(1);
    }

    // If the file is found, we start parsing the file:
    // First line should contain one integer, the dimension of the function F.
    fscanf(fp, "%zu", &n); // Set the dimension, n, for the truth table
    TruthTable *f = initTruthTable(n);

    // Set all the elements in the truth table. Should be 2^n elements.
    for (size_t i = 0; i < 1L << n; ++i) {
        fscanf(fp, "%zu", &f->elements[i]);
    }

    fclose(fp);
    return f;
}

size_t *mapPreImages(Partition *F, Partition *G) {
    // First we check for contradiction.
    if (F->numBuckets != G->numBuckets) return NULL;

    size_t *map = malloc(sizeof(size_t) * F->numBuckets); // The new map where F -> G

    // Loop over all the buckets and map F -> G with their pre-images.
    for (int i = 0; i < F->numBuckets; ++i) {
        for (int j = 0; j < G->numBuckets; ++j) {
            if (F->multiplicities[i] == G->multiplicities[j]) {
                // If the multiplicity in F[i] == G[j], we have found the multiplicities that maps to another
                map[i] = j;
                break;
            }
        }
    }
    return map;
}

void countElements(TruthTable *F, size_t *occurrences) {
    size_t dimension = F->n;
    for (size_t x = 0; x < 1L << dimension; ++x) {
        size_t y = F->elements[x];
        occurrences[y] += 1;
    }
}

bool outerPermutation(Partition *F, Partition *G, size_t n, size_t *basis, size_t *map, TruthTable *functionF,
                      TruthTable *functionG, bool affineSearch) {
    size_t *images = malloc(sizeof(size_t) * n); // The images of the basis elements under l
    size_t *generated = calloc(sizeof(size_t), 1L << n); // A partial truth table for l
    bool *generatedImages = calloc(sizeof(bool), 1L << n);
    bool foundSolution = false;

    /**
     * Create dictionaries indexing buckets by elements
     * For instance, fClassPosition[i] would be the index of the bucket w.r.t. F containing the element i.
     */
    size_t *fClass = createBucketRepresentation(F, n);
    size_t *gClass = createBucketRepresentation(G, n);

    // Recursively guess the values of l on the basis (essentially, a dfs with backtracking upon contradiction
    guessValuesOfL(0, basis, images, F, G, n, generated, generatedImages, fClass, gClass, map, &foundSolution,
                   functionF,
                   functionG, affineSearch);

    free(images);
    free(generated);
    free(generatedImages);
    free(fClass);
    free(gClass);
    return foundSolution;
}

void
guessValuesOfL(size_t k, size_t *basis, size_t *images, Partition *partitionF, Partition *partitionG, size_t n,
               size_t *generated, bool *generatedImages, size_t *fBucket, size_t *gBucket, size_t *map,
               bool *foundSolution, TruthTable *functionF, TruthTable *functionG, bool affineSearch) {
    if (*foundSolution) return;
    /**
     * If all basis elements have been assigned an image, and no contradictions have occurs, then we have found a
     * linear permutation preserving the partition. We reconstruct its truth table, and add it to the linked list
     * containing all permutations found by the search.
     */
    if (k == n) {
        TruthTable *currentL1 = initTruthTable(n);
        memcpy(currentL1->elements, generated, sizeof(size_t) * 1L << n);
        TruthTable *L1Inverse = inverse(currentL1); // L1^{-1}
        TruthTable *GPrime = compose(L1Inverse, functionG); // L1^{-1} * G = G'
        TruthTable *L2 = initTruthTable(n);
        L2->elements[0] = 0; // We know that the function is linear => L[0] -> 0

        if (innerPermutation(functionF, GPrime, basis, L2, affineSearch)) {
            /* At this point, we know (L1,L2) linear s.t. L1 * orthoderivativeF * L2 = orthoderivativeG */
            *foundSolution = true;
            printf("L1:\n");
            printTruthTable(currentL1);
            printf("L2:\n");
            printTruthTable(L2);
            destroyTruthTable(L1Inverse);
            destroyTruthTable(GPrime);
            destroyTruthTable(L2);
            destroyTruthTable(currentL1);
            return;
        }
        destroyTruthTable(currentL1);
        destroyTruthTable(L1Inverse);
        destroyTruthTable(GPrime);
        destroyTruthTable(L2);
        return;
    }
        /**
         * We then take the bucket of the same size from the partition with respect to G. We know that the image of the
         * basis element must belong to that bucket.
         */
    size_t posBucketG = map[fBucket[basis[k]]];

    // We now go through all possible choices from the bucket
    for (size_t ick = 0; ick < partitionG->bucketSizes[posBucketG]; ++ick) {
        size_t ck = partitionG->buckets[posBucketG][ick];

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
             if (partitionF->bucketSizes[fBucket[x]] != partitionG->bucketSizes[gBucket[y]]) {
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
            guessValuesOfL(k + 1, basis, images, partitionF, partitionG, n, generated, generatedImages, fBucket,
                           gBucket,
                           map, foundSolution, functionF, functionG, affineSearch);
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

size_t *createBucketRepresentation(Partition *F, size_t n) {
    // Loop over each bucket and set the bucket pos for each value
    size_t *class = malloc(sizeof(size_t) * 1L << n);
    for (size_t i = 0; i < F->numBuckets; ++i) {
        for (size_t j = 0; j < F->bucketSizes[i]; ++j) {
            size_t value = F->buckets[i][j];
            class[value] = i; // Set the bucket number in the position of the value
        }
    }
    return class;
}

bool *computeSetOfTs(TruthTable *F, const size_t x) {
    size_t dimension = F->n;
    bool *map = calloc(sizeof(bool), 1L << dimension);
    for (size_t y = 0; y < 1L << dimension; ++y) {
        size_t t = F->elements[x] ^ F->elements[y] ^ F->elements[x ^ y];
        map[t] = true;
    }
    return map;
}

Node *computeRestrictedDomains(TruthTable *F, const bool *map) {
    size_t dimension = F->n;
    bool *domain = calloc(sizeof(bool), 1L << dimension);
    for (size_t i = 0; i < 1L << dimension; ++i) {
        domain[i] = true;
    }
    for (size_t t = 0; t < 1L << dimension; ++t) {
        if (map[t]) {
            bool *tempSet = calloc(sizeof(bool), 1L << dimension);
            for (size_t x = 0; x < 1L << dimension; ++x) {
                for (size_t y = 0; y < 1L << dimension; ++y) {
                    if (t == (F->elements[x] ^ F->elements[y] ^ F->elements[x ^ y])) {
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

bool innerPermutation(TruthTable *F, TruthTable *G, const size_t *basis, TruthTable *L2, bool affineSearch) {
    size_t dimension = F->n;
    Node **restrictedDomains = malloc(sizeof(Node **) * (dimension + 1));
    bool result;

    for (size_t i = 0; i < dimension; ++i) {
        bool *map = computeSetOfTs(G, basis[i]);
        restrictedDomains[i] = computeRestrictedDomains(F, map);
        free(map);
    }

    size_t *values = malloc(sizeof(size_t) * dimension);

    size_t constant_term = G->elements[0];
    /* Guess of constant term of L2 */

    if (affineSearch) {
        for (size_t c2 = 0; c2 < 1L << dimension; ++c2) {
            /* Only consider preimages of G(0) */
            if (F->elements[c2] != constant_term) {
                continue;
            }

            TruthTable *newG = initTruthTable(dimension);
            for (size_t x = 0; x < 1L << dimension; ++x) {
                newG->elements[x ^ c2] = G->elements[x];
            }

            result = dfs(restrictedDomains, 0, values, F, newG, L2, basis);
            if (result) {
                /* If we get a result, we have to add the constant to the linear function that we found in dfs, and check if
             * F * l2 + c = G */
                if (c2 != 0) {
                    for (int x = 0; x < 1L << dimension; ++x) {
                        L2->elements[x] ^= c2;
                    }
                }

                /* If everything went smoothly, we should have aPrime == G */
                TruthTable *aPrime = compose(F, L2);

                destroyTruthTable(aPrime);
                destroyTruthTable(newG);
                break;
            }
            destroyTruthTable(newG);
        }
    } else {
        result = dfs(restrictedDomains, 0, values, F, G, L2, basis);
        if (result) {
            /* If everything went smoothly, we should have aPrime == G */
            TruthTable *aPrime = compose(F, L2);

            destroyTruthTable(aPrime);
        }
    }

    free(values);

    for (size_t i = 0; i < dimension; ++i) {
        destroyNodes(restrictedDomains[i]);
    }
    free(restrictedDomains);
    return result;
}

bool dfs(Node **domains, size_t k, size_t *values, TruthTable *F, TruthTable *G, TruthTable *L2, const size_t *basis) {
    size_t dimension = F->n;
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
            size_t new_value = L2->elements[linear_combination] ^ current->data;
            L2->elements[new_input] = new_value;
            /* Check for a violation of F * L2 = G */
            if (F->elements[new_value] != G->elements[new_input]) {
                /* Something is wrong, backtrack */
                problem = true;
                break;
            }
        }
        if (!problem) {
            if (dfs(domains, k + 1, values, F, G, L2, basis)) return true;
        }
        current = current->next;
    }
    return false;
}
