// ==================================================================
// Scanner.h
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Declaration of type Scanner.
// ==================================================================

#ifndef SCANNER_H
#define SCANNER_H

#include "Template.h"
#include "Molecule.h"

// ==================================================================
// Forward declarations
// ==================================================================
// Scanner					A Jess scan on a molecule
// ==================================================================

typedef struct _Scanner Scanner;

// ==================================================================
// Methods of type Scanner
// ==================================================================
// create(M,T,r)			Create object to scan M with template T
// free(S)					Free memory associated with S
// next(S)					Next result (an array of Atoms)
// ==================================================================

extern Scanner *Scanner_create(Molecule*,Template*,double,double);
extern void Scanner_free(Scanner*);
extern Atom **Scanner_next(Scanner*, int);
extern double Scanner_rmsd(Scanner*);

// ==================================================================

#endif

