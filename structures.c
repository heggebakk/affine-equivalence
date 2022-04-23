#include <stdbool.h>
#include <memory.h>
#include <stdlib.h>
#include <stdio.h>
#include "structures.h"
#include "affine.h"

TruthTable *initTruthTable(size_t n) {
    TruthTable *tt = malloc(sizeof(TruthTable));
    tt->n = n; // Dimension of the function
    tt->elements = malloc(sizeof(size_t) * 1L << n); // Allocate memory to fit all the elements from the function
    return tt;
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

void writeTruthTable(TruthTable *tt, FILE *filepath) {
    for (int i = 0; i < 1L << tt->n; ++i) {
        if (i < (1L << tt->n) - 1) {
            fprintf(filepath, "%zu ", tt->elements[i]);
        } else {
            fprintf(filepath, "%zu\n", tt->elements[i]);
        }
    }
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

Partition *partitionTt(TruthTable *tt) {
    size_t dimension = tt->n;
    size_t *multiplicities = malloc(sizeof(size_t) * 1L << dimension);
    Partition *partition = initPartition(dimension);
    memset(multiplicities, 0, sizeof (size_t) * 1L << dimension);
    countMultiplicities(tt, multiplicities);

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

BucketsMap *initBucketsMap() {
    BucketsMap *new = malloc(sizeof(BucketsMap));
    new->numOfMappings = 0;
    return new;
}

void destroyBucketsMap(BucketsMap *bucketsMap) {
    for (int i = 0; i < bucketsMap->numOfMappings; ++i) {
        free(bucketsMap->mappings[i]);
    }
    if (bucketsMap->numOfMappings != 0) {
        free(bucketsMap->mappings);
    }
    free(bucketsMap);
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