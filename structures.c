#include <stdbool.h>
#include <memory.h>
#include <stdlib.h>
#include <stdio.h>
#include "structures.h"
#include "affine.h"

TruthTable *initTruthTable(size_t dimension) {
    TruthTable *newTt = malloc(sizeof(TruthTable));
    newTt->dimension = dimension;
    newTt->elements = malloc(sizeof(size_t) * 1L << dimension);
    return newTt;
}

void destroyTruthTable(TruthTable *truthTable) {
    free(truthTable->elements);
    free(truthTable);
}

Partition *initPartition(size_t dimension) {
    Partition *partition = malloc(sizeof(Partition));
    partition->multiplicities = malloc(sizeof(size_t) * dimension);
    partition->bucketSizes = malloc(sizeof(size_t) * dimension);
    partition->buckets = malloc(sizeof(size_t **) * dimension);
    partition->numBuckets = 0;
}

Partition *partitionTt(TruthTable *truthTable) {
    Partition *partition;
    size_t dimension = truthTable->dimension;
    size_t *multiplicities = malloc(sizeof(size_t) * 1L << dimension);
    memset(multiplicities, 0, sizeof (size_t) * 1L << dimension);
    calculateMultiplicities(truthTable, multiplicities);
    for (int i = 0; i < 1L << dimension; ++i) {
        printf("%zu ", multiplicities[i]);
    }
    printf("\n");
    partition = initPartition(1L << dimension);

    for (int i = 0; i < 1L << dimension; ++i) {
        size_t numBuckets = partition->numBuckets; // init value = 0
        size_t multiplicity = multiplicities[i]; // Get i'th element in list of multiplicities

        // Check if multiplicity is in bucket, if false; add multiplicity to bucket
        bool multiplicityIsInBucket = false;
        for (int b = 0; b < numBuckets; ++b) {
            if (partition->multiplicities[b] == multiplicity) {
                multiplicityIsInBucket = true;
                partition->buckets[b][partition->bucketSizes[b]] = i;
                partition->bucketSizes[b] += 1;
                break;
            }
        }
        if (!multiplicityIsInBucket) {
            // Add a new bucket to the buckets list
            partition->buckets[numBuckets] = malloc(sizeof(size_t) * 1L << dimension);
            partition->bucketSizes[numBuckets] = 1;
            partition->multiplicities[numBuckets] = multiplicity;
            partition->buckets[numBuckets][0] = i;
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

void addNode(Node *head, int data) {
    Node *newNode = initNode();
    newNode->data = data;
    newNode->next = head->next;
    head->next = newNode;
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
