#include "adjoint.h"

_Bool dot(size_t a, size_t b) {
    return __builtin_popcountl(a & b) % 2;
}

_Bool is_it_really_adjoint(TruthTable * L, TruthTable *La) {
    size_t n = L->n;
    for(size_t x = 0; x < (1L << n); ++x) {
        for(size_t y = 0; y < (1L << n); ++y) {
            if (dot(L->elements[x],y) != dot(x,La->elements[y])) {
                return false;
            }
        }
    }
    return true;
}

TruthTable * recursive_adjoint(TruthTable * L, TruthTable * La, _Bool * assigned_elements, size_t i, size_t n) {
//    printf("Entering stage number %lu\n", i);
//    printf("Previously assigned: %lu -> %lu\n", i-1, La->elements[(1L << (i-1))]);
    if (i >= n) {
        return La;
    }

    /* Guess the value of basis element number i, i.e. 2^i */
    size_t basis_element = (1L << i);
    for(size_t value = 1; value < (1L << n); ++value) {
        if(assigned_elements[value]) {
            continue;
        }

        /* Check for contradictions */
        _Bool problem = false;
        for(size_t lincomb = 0; lincomb < (1L << (i)); ++lincomb) {
            /* Now we know the value lincomb XOR basis_element */
            size_t new_input = lincomb ^ basis_element;
            size_t new_output = La->elements[lincomb] ^ value;
            for(size_t x = 0; x < (1L << n); ++x) {
                if (dot(L->elements[x],new_input) != dot(x,new_output)) {
                    problem = true;
                    break;
                }
                La->elements[new_input] = new_output;
                //assigned_elements[new_output] = true;
            }
            if(problem) {
                break;
            }
        }

        if(!problem) {
            TruthTable * potential_result = recursive_adjoint(L, La, assigned_elements, i+1, n);
            if(potential_result) {
                return potential_result;
            }
        }

        /* Reset Boolean map */
        for(size_t lincomb = 0; lincomb < (1L << (i)); ++lincomb) {
            size_t new_input = lincomb ^ basis_element;
            size_t new_output = La->elements[lincomb] ^ value;
            assigned_elements[new_output] = false;
        }

    }
    return 0;
}

TruthTable * adjoint(TruthTable * L) {
    size_t n = L->n;
    TruthTable * La = initTruthTable(n);
    La->elements[0] = 0;
    _Bool * assigned_elements = calloc( sizeof(_Bool), (1L << n) );
    assigned_elements[0] = true;
    return recursive_adjoint(L, La, assigned_elements, 0, n);
}

