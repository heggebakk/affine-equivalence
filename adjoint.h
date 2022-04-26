#ifndef ADJOINT_H
#define ADJOINT_H

#include <stdio.h>
#include "structures.h"

_Bool dot(size_t a, size_t b);

_Bool is_it_really_adjoint(TruthTable * L, TruthTable *La);

/* Finds the adjoint operator of the linear permutation L */
TruthTable * adjoint(TruthTable * L);

#endif 
