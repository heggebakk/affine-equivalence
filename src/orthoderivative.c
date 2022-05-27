#include "orthoderivative.h"

/**
 * @author Nikolay S. Kaleyski
 */

#define POPCOUNT_FUNCTION __builtin_popcountl

TruthTable *orthoderivative(TruthTable *F) {
    size_t dimension = F->n;
    size_t entries = 1L << dimension;
    TruthTable *od = initTruthTable(dimension);

    /* Compute each element of the orthoderivative manually: o(a) must be such that
     * the dot product o(a) * (F(x) + F(a+x) + F(a) + F(0)) is equal to 0 for all x. */
    od->elements[0] = 0;
    _Bool problem = false;
    size_t count = 0;

    for (size_t a = 1; a < entries; ++a) {
        for (size_t possible_value = 1; possible_value < entries; ++possible_value) {
            problem = false;
            for (size_t x = 0; x < entries; ++x) {
                size_t derivative = F->elements[0] ^ F->elements[a] ^ F->elements[x] ^ F->elements[x ^ a];
                if (POPCOUNT_FUNCTION(possible_value & derivative) % 2) {
                    problem = true;
                    break;
                }
            }
            if (!problem) {
                od->elements[a] = possible_value;
                count +=1 ;
                break;
            }
        }
    }
    if (count < entries - 1) {
        printf("Orhtoderivative not working for one of the functions:\n");
        printTruthTable(F);
        exit(0);
    }
    return od;
}
