// ==================================================================
// Jess.h
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Declaration of types Jess, JessQuery and their methods.
// ==================================================================

#ifndef JESS_H
#define JESS_H

#include "Super.h"
#include "Molecule.h"
#include "Template.h"
#include "Atom.h"

// ==================================================================
// Forward declarations
// ==================================================================
// Jess					The Jess module
// JessQuery			The Jess query type
// ==================================================================

typedef struct _Jess Jess;
typedef struct _JessQuery JessQuery;

// ==================================================================
// Methods of type Jess
// ==================================================================
// create()				Create a Jess module
// free(J)				Free Jess object J (AND all templates)
// addTemplate(J,T)		Add T to the list of templates for J
// query(J,M,t)			Start a query on M using J threshold t
// ==================================================================

extern Jess *Jess_create(void);
extern void Jess_free(Jess*);
extern void Jess_addTemplate(Jess*,Template*);
extern JessQuery *Jess_query(Jess*,Molecule*,double,double);

// ==================================================================
// Methods of type JessQuery
// ==================================================================
// free(Q)				Frees the query object (NOT the molecule)
// next(Q)				Finds next result (true if successful)
// template(Q)			Returns the template for the hit
// molecule(Q)			Returns the molecule in which hit was found
// atoms(Q)				Array of atoms for the hit
// superposition(Q)		The superposition 
// ==================================================================

extern void JessQuery_free(JessQuery*);
extern int JessQuery_next(JessQuery*, int);
extern Template *JessQuery_template(JessQuery*);
extern const Molecule *JessQuery_molecule(JessQuery*);
extern Atom **JessQuery_atoms(JessQuery*);
extern Superposition *JessQuery_superposition(JessQuery*);

// ==================================================================

#endif

