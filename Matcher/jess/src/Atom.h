// ==================================================================
// Atom.h
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Declaration of type Atom and its methods.
// ==================================================================

#ifndef ATOM_H
#define ATOM_H

// ==================================================================
// Forward declarations
// ==================================================================
// Atom					A PDB ATOM record
// ==================================================================

typedef struct _Atom Atom;

// ==================================================================
// type Atom
// ==================================================================
// The fields correspond more or less exactly to those found in the
// PDB ATOM record spec. Notable exceptions are the chainID field
// which is set to '0' if it is blank (cf CATH) and the coordinate
// fields which are merged into an array Atom::x[3]. Finally, all
// text fields have blanks replaced by underscores.
// ==================================================================

struct _Atom
{
	int serial;
	char name[5];
	char altLoc;
	char resName[4];
	//Riziotis edit
	char chainID1;
	char chainID2;
	//char chainID;
	int resSeq;
	char iCode;
	double x[3];
	double occupancy;
	double tempFactor;
	char segID[4];
	char element[3];
	int charge;
};

// ==================================================================
// Methods of type Atom
// ==================================================================
// parse(A,s)			Parse string s as a PDB ATOM; true=>success
// ==================================================================

extern int Atom_parse(Atom*,const char*);

// ==================================================================

#endif


