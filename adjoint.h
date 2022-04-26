#ifndef ADJOINT_H
#define ADJOINT_H

#include <stdio.h>
#include "structures.h"

_Bool dot(size_t a, size_t b);

/* Finds the adjoint operator of the linear permutation L */
TruthTable * adjoint(TruthTable * L);

#endif 
