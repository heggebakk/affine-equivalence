#ifndef AFFINE_STRUCTURES_H
#define AFFINE_STRUCTURES_H

#include <stdbool.h>
#include <memory.h>
#include <stdlib.h>
#include <stdio.h>

/**
 * A truth table which holds the information of the dimension of the truth table and all its elements.
 */
typedef struct TruthTable {
    size_t dimension;
    size_t *elements;
} TruthTable;

/**
 * Initialize a new TruthTable
 * @param dimension The dimension of the TruthTable
 * @return The new TruthTable
 */
TruthTable *initTruthTable(size_t dimension);

void printTruthTable(TruthTable *truthTable);

/**
 * Free memory for struct TruthTable
 * @param truthTable The truth table to be freed
 */
void destroyTruthTable(TruthTable *truthTable);

typedef struct Partition {
    size_t numBuckets;
    size_t *multiplicities;
    size_t *bucketSizes;
    size_t **buckets;
} Partition;

Partition *initPartition(size_t dimension);

void printPartition(Partition *partition);

Partition *partitionTt(TruthTable *truthTable);

void destroyPartition(Partition *partition);

typedef struct Node {
    size_t data; // The data for this node
    struct Node *next; // Pointer to the next node in the linked list
} Node;

Node *initNode();

void addNode(Node *head, int data);

size_t countNodes(Node *head);

void printNodes(Node *head);

typedef struct BucketsMap {
    size_t **domains;
    size_t **mappings;
    size_t numOfMappings;
} BucketsMap;

BucketsMap *initBucketsMap();

#endif //AFFINE_STRUCTURES_H
