#include "orthoderivative.h"

/**
 * @author Nikolay S. Kaleyski
 */

#define POPCOUNT_FUNCTION __builtin_popcountl

TruthTable *orthoderivative(TruthTable *f) {
    size_t dimension = f->n;
    size_t entries = 1L << dimension;
    TruthTable *od = initTruthTable(dimension);

    /* Compute each element of the orthoderivative manually: o(a) must be such that
     * the dot product o(a) * (F(x) + F(a+x) + F(a) + F(0)) is equal to 0 for all x. */
    od->elements[0] = 0;
    _Bool problem = false;

    for (size_t a = 1; a < entries; ++a) {
        for (size_t possible_value = 1; possible_value < entries; ++possible_value) {
            problem = false;
            for (size_t x = 0; x < entries; ++x) {
                size_t derivative = f->elements[0] ^ f->elements[a] ^ f->elements[x] ^ f->elements[x ^ a];
                if (POPCOUNT_FUNCTION(possible_value & derivative) % 2) {
                    problem = true;
                    break;
                }
            }
            if (!problem) {
                od->elements[a] = possible_value;
                break;
            }
        }
    }
    return od;
}
