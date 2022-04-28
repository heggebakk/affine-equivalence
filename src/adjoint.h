#ifndef ADJOINT_H
#define ADJOINT_H

#include <stdio.h>
#include "structures.h"
/**
 * The dot product of a and b
 */
_Bool dot(size_t a, size_t b);

/* Finds the adjoint operator of the linear permutation L */
/**
 * Finds the adjoint operator of the linear permutation L
 * @param L Linear permutation
 * @return The adjoint of a linear permutation L
 */
TruthTable *adjoint(TruthTable *L);

#endif
