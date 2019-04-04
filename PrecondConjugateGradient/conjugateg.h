#ifndef CONJUGATEG_H
#define CONJUGATEG_H

#include "problemsetup.h"
#include "vectorops.h"
#include "precondition.h"

void precond_conjugate_gradient(Grid *x, Grid *b);
void conjugate_gradient(Grid *x, Grid *b);

#endif
