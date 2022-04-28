#include "adjoint.h"

_Bool dot(size_t a, size_t b) {
    return __builtin_popcountl(a & b) % 2;
}

TruthTable *recursiveAdjoint(TruthTable *L, TruthTable *La, bool *assignedElements, size_t i, size_t n) {
    if (i >= n) {
        return La;
    }

    /* Guess the value of basis element number i, i.e. 2^i */
    size_t basisElement = 1L << i;
    for (size_t value = 1; value < 1L << n; ++value) {
        if (assignedElements[value]) {
            continue;
        }

        /* Check for contradictions */
        _Bool problem = false;
        for (size_t linearCombination = 0; linearCombination < 1L << i; ++linearCombination) {
            /* Now we know the value linearCombination XOR basisElement */
            size_t newInput = linearCombination ^ basisElement;
            size_t newOutput = La->elements[linearCombination] ^ value;
            for (size_t x = 0; x < 1L << n; ++x) {
                if (dot(L->elements[x], newInput) != dot(x, newOutput)) {
                    problem = true;
                    break;
                }
                La->elements[newInput] = newOutput;
            }
            if (problem) {
                break;
            }
        }

        if (!problem) {
            TruthTable *potentialResult = recursiveAdjoint(L, La, assignedElements, i + 1, n);
            if (potentialResult) {
                return potentialResult;
            }
        }

        /* Reset Boolean map */
        for (size_t linearCombination = 0; linearCombination < 1L << i; ++linearCombination) {
            size_t new_output = La->elements[linearCombination] ^ value;
            assignedElements[new_output] = false;
        }

    }
    return 0;
}

TruthTable *adjoint(TruthTable *L) {
    size_t n = L->n;
    TruthTable *LAdjoint = initTruthTable(n);
    LAdjoint->elements[0] = 0;
    _Bool *assignedElements = calloc(sizeof(_Bool), 1L << n);
    assignedElements[0] = true;
    TruthTable *adjointTt = recursiveAdjoint(L, LAdjoint, assignedElements, 0, n);
    free(assignedElements);
    return adjointTt;
}

