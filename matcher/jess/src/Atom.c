// ==================================================================
// Atom.c
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Implementation of methods of type Atom.
// ==================================================================

#include "Atom.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

// ==================================================================
// Methods of type Atom
// ==================================================================

static void Atom_copyToken(char *d,const char *s,int n)
{
	int i;

	// Copy string from s to d replacing blanks by the
	// underscore where appropriate.

	strncpy(d,s,n);
	d[n]=0;
	for(i=0; i<n; i++)
	{
		if(isspace(d[i])) d[i]='_';
	}
}

int Atom_parse(Atom *A, const char *s)
{
	int n;

	// Check the record type...

	if(strncmp(s,"ATOM",4) && strncmp(s,"HETATM",6)!=0)
	{
		return 0;
	}

	// Get string length and zero the atom structure

	n = strlen(s);
	memset(A,0,sizeof(Atom));

	// Get the extant fields...

	if(n>6)  A->serial = atoi(&s[6]); else return 1;
	if(n>12) Atom_copyToken(A->name,&s[12],4); else return 1;
	if(n>16) A->altLoc = s[16]; else return 1;
	if(n>17) Atom_copyToken(A->resName,&s[17],3); else return 1;
	//Riziotis edit
	if(n>20) A->chainID1 = s[20]; else return 1;
	if(n>21) A->chainID2 = isspace(s[21]) ? '0':s[21]; else return 1;
	//if(n>21) A->chainID = isspace(s[21]) ? '0':s[21]; else return 1;
	if(n>22) A->resSeq = atoi(&s[22]); else return 1;
	if(n>26) A->iCode = s[26]; else return 1;
	if(n>30) A->x[0] = atof(&s[30]); else return 1;
	if(n>38) A->x[1] = atof(&s[38]); else return 1;
	if(n>46) A->x[2] = atof(&s[46]); else return 1;
	if(n>54) A->occupancy = atof(&s[54]); else return 1;
	if(n>60) A->tempFactor = atof(&s[60]); else return 1;
	if(n>72) Atom_copyToken(A->segID,&s[72],4); else return 1;
	if(n>76) Atom_copyToken(A->element,&s[76],2); else return 1;
	if(n>78) A->charge = atoi(&s[78]); else return 1;

	// All's well :-)

	return 1;
}

// ==================================================================

