#include "structures.h"
#include "affine.h"

size_t *createBasis(size_t dimension);

void addConstant(TruthTable *tt, size_t c);

int main() {
    char *filename = "resources/q_6_1.tt";
    size_t dimension;
    size_t *basis;
    TruthTable *functionF = parseFile(filename);
//    TruthTable *functionG = createTruthTable(functionF);
    TruthTable *functionG = parseFile("resources/g.tt");
    printTruthTable(functionF);
    printTruthTable(functionG);
    Partition *partitionF = partitionTt(functionF);
    dimension = functionF->dimension;
    basis = createBasis(dimension);

    // Need to test for all possible constants, 0..2^n - 1.
    for (int c = 0; c < 1L << dimension; ++c) {
        _Bool foundSolution = false; /* for breaking out of nested loops */
        TruthTable *gPrime = initTruthTable(dimension);
        memcpy(gPrime->elements, functionG->elements, sizeof(size_t) * 1L << dimension);
        addConstant(gPrime, c);
//        printf("G':\n");
//        printTruthTable(gPrime);
        Partition *partitionG = partitionTt(gPrime);
        BucketsMap *bucketsMap = mapBuckets(partitionF, partitionG);
//        printf("num of maps: %zu\n", bucketsMap->numOfMappings);

        for (size_t map = 0; map < bucketsMap->numOfMappings; ++map) {
            // Calculate outer permutation
            TtNode *a1 = outerPermutation(partitionF, partitionG, dimension, basis, bucketsMap->domains[map]);
            size_t numPermutations = countTtNodes(a1);

            for (size_t i = 0; i < numPermutations; ++i) {
                TruthTable *a1Prime = getTtNode(a1, i);
                TruthTable *a1Inverse = inverse(a1Prime);
                TruthTable *gDoublePrime = compose(a1Inverse, gPrime);
                TruthTable *a2 = initTruthTable(dimension);
                a2->elements[0] = 0;

                if (innerPermutation(functionF, gDoublePrime, basis, a2)) {
                    printf("Hello!\n");
                    printf("a2:\n");
                    printTruthTable(a2);
                    foundSolution = true;
                }
                destroyTruthTable(a1Inverse);
                destroyTruthTable(gDoublePrime);
                destroyTruthTable(a2);

                if (foundSolution) break;
            }
            destroyTtNode(a1);

            if (foundSolution) break;
        }
        destroyTruthTable(gPrime);
        destroyBucketsMap(bucketsMap);
        destroyPartition(partitionG);

        if (foundSolution) break;
    }

    destroyTruthTable(functionF);
    destroyTruthTable(functionG);
    destroyPartition(partitionF);
    free(basis);

    return 0;
}

void addConstant(TruthTable *tt, size_t c) {
    size_t dimension = tt->dimension;
    for (size_t i = 0; i < 1L << dimension; ++i) {
        tt->elements[i] ^= c;
    }
}

size_t *createBasis(size_t dimension) {
    size_t *basis = malloc(sizeof(size_t) * (dimension));
    for (size_t i = 0; i < dimension; ++i) {
        basis[i] = 1L << i;
    }
    return basis;
}
