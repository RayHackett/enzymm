// ==================================================================
// Annulus.h
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Declaration of Annulus construction and destruction.
// ==================================================================

#ifndef ANNULUS_H
#define ANNULUS_H

#include "Region.h"

// ==================================================================
// Methods for Annulus manipulation
// ==================================================================
// create(u,a,b,d)		Make region {x in R^d : a <= |x-u| <= b }.
// free(A)				Free the region given (or use R->free)
// ==================================================================

extern Region *Annulus_create(double*,double,double,int);
extern void Annulus_free(Region*);

// ==================================================================

#endif

