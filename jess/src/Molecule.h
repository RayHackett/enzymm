// ==================================================================
// Molecule.h
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Declaration of type Molecule (a PDB reader).
// ==================================================================

#ifndef MOLECULE_H
#define MOLECULE_H

#include "Atom.h"
#include <stdio.h>

// ==================================================================
// Forward declarations
// ==================================================================
// Molecule					A complete molecule
// ==================================================================

typedef struct _Molecule Molecule;

// ==================================================================
// Methods of type Molecule
// ==================================================================
// create(file)				Create molecule from PDB file
// free(M)					Free memory associated with molecule M
// count(M)					Count number of atoms in the molecule
// atom(M,k)				Return pointer to atom k (see Atom.h)
// id(M)					The PDB code (if found)
// ==================================================================

extern Molecule *Molecule_create(FILE*,int,float);
extern void Molecule_free(Molecule*);
extern int Molecule_count(const Molecule*);
extern const Atom *Molecule_atom(const Molecule*,int);
extern const char *Molecule_id(const Molecule*);

// ==================================================================


#endif

