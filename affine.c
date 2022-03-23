#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include "structures.h"
#include <memory.h>
#include "affine.h"

void createBucketsMap(BucketsMap *bucketsMap, struct Node **domains, struct Partition *pG, size_t dimension,
                      size_t numOfMappings);

void selectRecursive(int i, bool *chosen, struct Node **domains, struct Partition *pG, size_t dimension,
                     BucketsMap *bucketsMap);

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
    for (int i = 0; i < 1L << dimension; ++i) {
        fscanf(fp, "%zu", &tt->elements[i]);
    }
    fclose(fp);
    return tt;
}

BucketsMap *mapBuckets(struct Partition *f, struct Partition *g, size_t dimension) {
    BucketsMap *bucketsMap = initBucketsMap();
    if (f->numBuckets != g->numBuckets) {
        printf("F and G is not affine\n");
        exit(0);
    }

    // Find all domains for the different buckets.
    struct Node **domains = malloc(sizeof(struct Node) * f->numBuckets);
    for (int i = 0; i < f->numBuckets; ++i) domains[i] = initNode();

    for (int i = 0; i < f->numBuckets; ++i) {
        bool sameSize = false;
        for (int j = 0; j < g->numBuckets; ++j) {
            if (f->bucketSizes[i] == g->bucketSizes[j]) {
                sameSize = true;
                addNode(domains[i], j);
            }
        }
        if (!sameSize) {
            printf("F and G is not affine!\n");
            exit(0);
        }
    }

    size_t numOfMappings = 1;
    bool *isCalculated = malloc(sizeof(bool) * f->numBuckets);
    memset(isCalculated, 0, sizeof(bool) * f->numBuckets);
    for (int i = 0; i < f->numBuckets; ++i) {
        Node *current = domains[i]->next;
        if (!isCalculated[current->data]) {
            isCalculated[current->data] = true;
            numOfMappings *= factorial(countNodes(domains[i]));
        }
    }
    free(isCalculated);
    printf("Num mappings = %zu", numOfMappings);
}

void calculateMultiplicities(TruthTable *truthTable, size_t *multiplicities) {
    size_t dimension = truthTable->dimension;
    for (int x = 0; x < 1L << dimension; ++x) {
        size_t y = truthTable->elements[x];
        multiplicities[y] += 1;
    }
}

size_t factorial(size_t value) {
    size_t factorial = 1;
    for (int i = 1; i < value + 1; ++i) {
        factorial *= i;
    }
    return factorial;
}