// ==================================================================
// Super.h
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Declaration of type Superposition and its methods
// ==================================================================

#ifndef SUPER_H
#define SUPER_H

// ==================================================================
// Forward declarations
// ==================================================================
// Superposition		Type respresenting a superposition
// ==================================================================

typedef struct _Superposition Superposition;

// ==================================================================
// Methods of type Superposition
// ==================================================================
// create()				Create an empty superposition object
// free(S)				Free superposition S and assocatited mem.
// associate(x,y)		Associate vectors x and y in the superpsn.
// count(S)				Return number of associated vector pairs
// rmsd(S)				Return the rmsd (computes if out-of-date)
// rmsd100(S)			As above but computes rmsd100 (see code)
// centroid(S,k)		Centroid of left or right set of points
// rotation(S)			Returns rotation matrix (as 9 doubles)
// ==================================================================

extern Superposition *Superposition_create(void);
extern void Superposition_free(Superposition*);
extern void Superposition_align(Superposition*,const double*,const double*);
extern int Superposition_count(const Superposition*);
extern double Superposition_rmsd(Superposition*);
extern double Superposition_rmsd100(Superposition*);
extern const double *Superposition_centroid(Superposition*,int);
extern const double *Superposition_rotation(Superposition*);

// ==================================================================

#endif

