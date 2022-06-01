#include <stdbool.h>
#include <memory.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "structures.h"
#include "equivalence.h"

TruthTable *initTruthTable(size_t n) {
    TruthTable *tt = malloc(sizeof(TruthTable));
    tt->n = n; // Dimension of the function
    tt->elements = malloc(sizeof(size_t) * 1L << n); // Allocate memory to fit all the elements from the function
    return tt;
}

void add(TruthTable *dest, TruthTable *src) {
    for (size_t x = 0; x < 1L << dest->n; ++x) {
        dest->elements[x] ^= src->elements[x]; // F[x] = F[x] + G[x]
    }
}

TruthTable *compose(TruthTable *f, TruthTable *g) {
    size_t dimension = f->n;
    TruthTable *result = initTruthTable(dimension);
    for (size_t x = 0; x < 1L << dimension; ++x) {
        result->elements[x] = f->elements[g->elements[x]]; // F[G[x]]
    }
    return result;
}

TruthTable *inverse(TruthTable *f) {
    size_t dimension = f->n;
    TruthTable *inverse = initTruthTable(dimension);
    for (size_t x = 0; x < 1L << dimension; ++x) {
        size_t y = f->elements[x];
        inverse->elements[y] = x;
    }
    return inverse;
}

TruthTable *randomLinearFunction(size_t n) {
    size_t entries = 1L << n;
    size_t listGenerated[entries];
    listGenerated[0] = 0;
    size_t basisImages[n];
    for (size_t i = 0; i < n; ++i) {
        size_t j = rand() % entries;
        basisImages[i] = j;
        for (int k = 0; k < 1L << i; ++k) {
            listGenerated[(1L << i) + k] = listGenerated[k] ^ j;
        }
    }
    TruthTable *newFunction = initTruthTable(n);
    memcpy(newFunction->elements, listGenerated, sizeof(size_t) * entries);
    return newFunction;
}

TruthTable *randomLinearPermutation(size_t n) {
    size_t entries = 1L << n;
    bool generated[entries];
    size_t listGenerated[entries];
    generated[0] = true;
    for (size_t i = 1; i < entries; ++i) {
        generated[i] = false;
    }
    listGenerated[0] = 0;

    size_t basisImages[n];
    for (int i = 0; i < n; ++i) {
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
    // Add a random constant to the new function, to create an affine function
    TruthTable *newFunction = initTruthTable(n);
    memcpy(newFunction->elements, listGenerated, sizeof(size_t) * entries);
    return newFunction;
}

TruthTable *createAffineTruthTable(TruthTable *F) {
    size_t n = F->n;
    size_t entries = 1L << n;

    // A1 * F * A2 + A = G
    TruthTable *A1 = randomLinearPermutation(n);
    TruthTable *A2 = randomLinearPermutation(n);
    TruthTable *A = randomLinearFunction(n);
    // Copy A1 A2 to L1 L2 before we add random constants to A1 A2 A, to create affine functions.

    // A random constant c, where c is in 2^n
    size_t constant1 = rand() % entries;
    size_t constant2 = rand() % entries;
    size_t constant3 = rand() % entries;
    for (int i = 0; i < entries; ++i) {
        // Add the constant
        A1->elements[i] ^= constant1;
        A2->elements[i] ^= constant2;
        A->elements[i] ^= constant3;
    }

    TruthTable *FComposeA2 = compose(F, A2);
    TruthTable *G = compose(A1, FComposeA2);
    add(G, A);

    destroyTruthTable(A1);
    destroyTruthTable(A2);
    destroyTruthTable(A);
    destroyTruthTable(FComposeA2);

    return G;
}

TruthTable *createLinearFunction(TruthTable *F) {
    size_t n = F->n;
    TruthTable *L1 = randomLinearPermutation(n);
    TruthTable *L2 = randomLinearPermutation(n);

    TruthTable *temp = compose(F, L2);
    TruthTable *G = compose(L1, temp);

    destroyTruthTable(L1);
    destroyTruthTable(L2);
    destroyTruthTable(temp);

    return G;
}

void printTruthTable(TruthTable *tt) {
    for (int i = 0; i < 1L << tt->n; ++i) {
        if (i < (1L << tt->n) - 1) {
            printf("%zu ", tt->elements[i]);
        } else {
            printf("%zu\n", tt->elements[i]);
        }
    }
    printf("\n");
}

void destroyTruthTable(TruthTable *tt) {
    free(tt->elements);
    free(tt);
}

Partition *initPartition(size_t n) {
    Partition *partition = malloc(sizeof(Partition));
    partition->multiplicities = malloc(sizeof(size_t) * n); // Malloc n lists
    partition->bucketSizes = malloc(sizeof(size_t) * n); // Malloc n lists
    partition->buckets = malloc(sizeof(size_t **) * n); // Malloc n lists of bucket lists.
    partition->numBuckets = 0;
}

void printPartitionBuckets(Partition *partition) {
    for (int i = 0; i < partition->numBuckets; ++i) {
        for (int j = 0; j < partition->bucketSizes[i]; ++j) {
            if (j == partition->bucketSizes[i] - 1) {
                printf("%zu\n", partition->buckets[i][j]);
            } else {
                printf("%zu ", partition->buckets[i][j]);
            }
        }
    }
    printf("\n");
}

void printPartition(Partition *F) {
    for (int i = 0; i < F->numBuckets; ++i) {
        printf("Bucket #%d = %zu\n", i, F->multiplicities[i]);
    }
}

Partition *partitionTt(TruthTable *tt) {
    size_t dimension = tt->n;
    size_t *multiplicities = calloc(sizeof(size_t), 1L << dimension);
    Partition *partition = initPartition(dimension);
    countElements(tt, multiplicities);

    for (int i = 0; i < 1L << dimension; ++i) {
        size_t numBuckets = partition->numBuckets; // init value = 0
        size_t multiplicity = multiplicities[i]; // Get i'th element in list of multiplicities

        // Check if multiplicity is in bucket, if false; add multiplicity to bucket
        bool multiplicityIsInBucket = false;
        for (int b = 0; b < numBuckets; ++b) {
            if (partition->multiplicities[b] == multiplicity) {
                multiplicityIsInBucket = true;
                size_t thisBucket = partition->bucketSizes[b];
                partition->buckets[b][thisBucket] = i; // Add the element to the correct bucket
                partition->bucketSizes[b] += 1;
                break;
            }
        }
        if (!multiplicityIsInBucket) {
            // The current multiplicity is not in the lists; add a new bucket to the buckets list
            partition->buckets[numBuckets] = malloc(sizeof(size_t) * 1L << dimension); // Allocate memory for a new bucket list
            partition->bucketSizes[numBuckets] = 1; // Increase the number of buckets by 1
            partition->multiplicities[numBuckets] = multiplicity; // Add the multiplicity of the new bucket
            partition->buckets[numBuckets][0] = i; // Add the element that belongs to this bucket
            partition->numBuckets += 1;
        }
    }
    free(multiplicities);
    return partition;
}

void destroyPartition(Partition *partition) {
    for (int i = 0; i < partition->numBuckets; ++i) {
        free(partition->buckets[i]);
    }
    free(partition->bucketSizes);
    free(partition->multiplicities);
    free(partition->buckets);
    free(partition);
}

Node *initNode() {
    Node *newNode = malloc(sizeof(Node));
    newNode->data = 0;
    newNode->next = NULL;
    return newNode;
}

void addNode(Node *head, size_t data) {
    Node *newNode = initNode();
    newNode->data = data;
    newNode->next = head->next;
    head->next = newNode;
}
size_t getNode(Node *head, size_t index) {
    Node *current = head;
    for (size_t i = 0; i < index; ++i) {
        current = current->next;
    }
    return current->data;
}

size_t countNodes(Node *head) {
    size_t count = 0;
    Node *current = head->next;

    while (current != NULL) {
        count += 1;
        current = current->next;
    }
    return count;
}

void printNodes(Node *head) {
    if (head->next == NULL) {
        printf("Linked list is empty.\n");
        return;
    }
    Node *current = head->next;
    while (current != NULL) {
        printf("%zu ", current->data);
        current = current->next;
    }
    printf("\n");
}

void destroyNodes(Node *head) {
    Node *current = NULL;
    while (head != NULL) {
        current = head;
        head = head->next;
        free(current);
    }
}

TtNode *initTtNode() {
    TtNode *newNode = malloc(sizeof(TtNode));
    newNode->data = NULL;
    newNode->next = NULL;
    return newNode;
}

void addTtNode(TtNode *head, TruthTable *data) {
    size_t dimension = data->n;
    if (head->data == NULL) {
        head->data = initTruthTable(dimension);
        memcpy(head->data->elements, data->elements, sizeof(size_t) * 1L << dimension);
        return;
    }
    TtNode *newNode = malloc(sizeof(TtNode));
    newNode->data = initTruthTable(dimension);
    memcpy(newNode->data->elements, data->elements, sizeof(size_t) * 1L << dimension);
    newNode->next = head->next;
    head->next = newNode;
}

TruthTable *getTtNode(TtNode *head, size_t index) {
    TtNode *current = head;
    for (size_t i = 0; i < index; ++i) {
        current = current->next;
    }
    return current->data;
}

size_t countTtNodes(TtNode *head) {
    if (head->data == NULL) return 0;
    size_t count = 1;
    TtNode *current = head;
    while (current->next != NULL) {
        count += 1;
        current = current->next;
    }
    return count;
}

void destroyTtNode(TtNode *head) {
    if (head->data == NULL) {
        free(head);
        return;
    }
    while (head->data != NULL) {
        TtNode *current = head;
        head = head->next;
        destroyTruthTable(current->data);
        free(current);
        if (head == NULL) {
            free(head);
            break;
        }
    }
}

RunTimes *initRunTimes() {
    RunTimes *newTime = malloc(sizeof(RunTimes));
    newTime->total = 0.0;
    return newTime;
}

double stopTime(double runTime, clock_t startParsing) {
    runTime += (double) (clock() - startParsing) / CLOCKS_PER_SEC;
    return runTime;
}

void printTimes(RunTimes *runTimes) {
    printf("Total time spent: %f \n", runTimes->total);
}

void destroyRunTimes(RunTimes *runTimes) {
    free(runTimes);
}

void addConstant(TruthTable *F, size_t c) {
    size_t dimension = F->n;
    for (size_t i = 0; i < 1L << dimension; ++i) {
        F->elements[i] ^= c;
    }
}

size_t *createStandardBasis(size_t n) {
    size_t *basis = malloc(sizeof(size_t) * n);
    for (size_t i = 0; i < n; ++i) {
        basis[i] = 1L << i;
    }
    return basis;
}

