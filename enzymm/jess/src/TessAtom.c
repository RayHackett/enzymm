// ==================================================================
// TessAtom.c
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Implementation of type TessAtom.
// ==================================================================

#include "TessAtom.h"
#include "Atom.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

// ==================================================================
// Local type TessAtom
// ==================================================================
// code					The code describing the match mode
// resSeq				The residue sequence number
// nameCount			Number of atom names present
// resNameCount			Number of residue names present
// chainID				The chain field
// name[k]				kth atom name alternate
// resName[k]			kth residue name alternate
// pos[k]				kth coordinate of position of atom
// distWeight[k]			kth atom distance threshold modifier (weight)
// ==================================================================

struct _TessAtom
{
	int code;
	int resSeq;
	int nameCount;
	int resNameCount;
	char chainID1; //Riziotis edit
	char chainID2;
	//char chainID; 
	char **name;
	char **resName;
	double pos[3];
	double distWeight;
};

// ==================================================================
// Methods of type TessAtom (ARGGH!!!)
// ==================================================================

static const char *TessAtom_convertResidue(char p)
{
	// Convert from single-letter code to three-letter code.

	int k;
	struct _dummy {char one; const char *three; };
	static const struct _dummy table[21] =
	{
		{ 'A',"ALA" },
		{ 'C',"CYS" },
		{ 'D',"ASP" },
		{ 'E',"GLU" },
		{ 'F',"PHE" },
		{ 'G',"GLY" },
		{ 'H',"HIS" },
		{ 'I',"ILE" },
		{ 'K',"LYS" },
		{ 'L',"LEU" },
		{ 'M',"MET" },
		{ 'N',"ASN" },
		{ 'P',"PRO" },
		{ 'Q',"GLN" },
		{ 'R',"ARG" },
		{ 'S',"SER" },
		{ 'T',"THR" },
		{ 'V',"VAL" },
		{ 'W',"TRP" },
		{ 'Y',"TYR" },
		{ 'X',"XXX" }
	};

	for(k=0; k<21 && p!=table[k].one; k++);
	return k<21 ? table[k].three:NULL;
}

TessAtom *TessAtom_create(const char *s)
{
	TessAtom *A;
	Atom a;
	int i,k,m;
	int rc,ac;
	int rc1,ac1;
	int nest;
	int rq;
	const char *q;
	const char *tmp;
	void *p;


	// 0. Parse the record as a standard PDB atom. We must
	// at least have all fields up to the last coord field.

	if(strlen(s)<54 || !Atom_parse(&a,s))
	{
		return NULL;
	}

	// 1. Find the number of extra atom names or
	// residue names...

	q=&s[54]; // Was 54, 66 when we parse distance threshold from temperature field
	k=strlen(q);
	nest=0;
	rc=ac=1;
	for(m=0; m<k; m++)
	{
		if(q[m]=='(')
		{
			if(nest!=0) return NULL;
			nest=1;
			ac++;
		}

		else if(q[m]==')')
		{
			if(nest!=1) return NULL;
			nest=0;
		}

		else if(nest==0 && TessAtom_convertResidue(q[m]))
		{
			rc++;
		}
	}
	if(nest!=0) return NULL;

	// Now, allocate sufficient memory for all
	// the atom names and residue names plus
	// the TessAtom structure itself...

	rq = sizeof(TessAtom);
	rq += sizeof(char*)*(ac+rc);
	rq += sizeof(char)*(5*ac+4*rc);
	A = (TessAtom*)calloc(1,rq);

	// Set up basic fields...

	A->code = a.serial;
	A->resSeq = a.resSeq;
	A->pos[0] = a.x[0];
	A->pos[1] = a.x[1];
	A->pos[2] = a.x[2];
	A->distWeight = a.tempFactor;
	A->chainID1 = a.chainID1;
	A->chainID2 = a.chainID2;
	//A->chainID = a.chainID;
	A->nameCount=ac;
	A->resNameCount=rc;

	// Set up all the pointers to the
	// residue name and atom name fields

	p=&A[1];
	A->name=p;
	p+=sizeof(char*)*ac;
	for(m=0; m<ac; m++)
	{
		A->name[m]=p;
		p+=5;
	}

	A->resName=p;
	p+=sizeof(char*)*rc;
	for(m=0; m<rc; m++)
	{
		A->resName[m]=p;
		p+=4;
	}

	// Copy the name and resName fields into
	// the arrays at index 0.

	strncpy(A->name[0],a.name,4);
	strncpy(A->resName[0],a.resName,3);

	// Finally, get all the extra fields at the end
	// of the PDB record...

	q=&s[54];
	rc1=1;
	ac1=1;
	while(rc1<rc || ac1<ac)
	{
		if(isspace(*q))
		{
			q++;
			continue;
		}

		else if(*q=='(')
		{
			q++;

			// Find the length of the atom name field (up
			// to 4 chars long...)

			for(m=0; m<4 && q[m]!=')'; m++);

			if(q[m]!=')')
			{
				// This is an error - the atom field may
				// be up to 4 chars long only...

				free(A);
				return NULL;
			}
			

			// Copy the atom name and replace blanks with
			// underscores...

			strncpy(A->name[ac1],q,m);
			for(i=0; i<4; i++)
			{
				if(isspace(A->name[ac1][i]) || !A->name[ac1][i])
				{
					A->name[ac1][i]='_';
				}
			}
			ac1++;

			// Skip to just after the closing bracket.

			q=&q[m+1];
			
		}

		else
		{
			// It must be a single-letter residue code.
			// Check it out...
			
			for(m=0; m<=4; m++){
				if(isspace(q[m])){
					break;
				}
				tmp = TessAtom_convertResidue(q[m]);
				if(!tmp)
				{
					free(A);
					return NULL;
				}

				strcpy(A->resName[rc1++],tmp);
			}
		}
	}

	// Sorted!

	return A;
}

//Riziotis edit
char TessAtom_chainID1(const TessAtom *A)
{
	return A->chainID1;
}

char TessAtom_chainID2(const TessAtom *A)
{
	return A->chainID2;
}

double TessAtom_distWeight(const TessAtom *A)
{
	return A->distWeight;
}
  
//char TessAtom_chainID(const TessAtom *A)
//{
//	return A->chainID;
//}

const double *TessAtom_position(const TessAtom *A)
{
	return A->pos;
}

void TessAtom_free(TessAtom *A)
{
	if(A) free(A);
}

int TessAtom_resSeq(const TessAtom *A)
{
	return A->resSeq;
}

static int TessAtom_isCarbon(const Atom *A)
{
	return A->name[0]=='_' && A->name[1]=='C' ? 1:0;
}

//Riziotis 
static int TessAtom_isHydrogen(const Atom *A)
{
	return A->name[0]=='_' && A->name[1]=='H' ? 1:0;
}

static int TessAtom_isInSamePosition(const TessAtom *T, const Atom *A)
{
	int k;

	for(k=0; k<T->nameCount; k++)
	{
		if(A->name[0]=='_')
		{
			if(A->name[2]==T->name[k][2]) return 1;
		}
		else
		{
			if(A->name[1]==T->name[k][1] && A->name[2]==T->name[k][2]) return 1;
		}
	}

	return 0;
}
//End Riziotis

static int TessAtom_isMainChain(const Atom *A)
{
	if(strcasecmp(A->name,"_CA_")==0) return 1;
	if(strcasecmp(A->name,"_N__")==0) return 1;
	if(strcasecmp(A->name,"_O__")==0) return 1;

	return 0;
}

static int TessAtom_matchName(const TessAtom *T, const Atom *A)
{
	int k;

	for(k=0; k<T->nameCount; k++)
	{
		if(strcasecmp(A->name,T->name[k])==0) return 1;
	}

	return 0;
}

static int TessAtom_matchResName(const TessAtom *T, const Atom *A)
{
	int k;

	for(k=0; k<T->resNameCount; k++)
	{
		if(strcasecmp(A->resName,T->resName[k])==0) return 1;
	}

	return 0;
}

static int TessAtom_matchType(const TessAtom *T, const Atom *A)
{
	// There may still be some problems with this matching routine.
	// But I don't give a flying f**k.

	int k;

	for(k=0; k<T->nameCount; k++)
	{
		if(A->name[0]=='_')
		{
			if(A->name[1]==T->name[k][1]) return 1;
		}
		else
		{
			if(A->name[0]==T->name[k][0] && A->name[1]==T->name[k][1]) return 1;
		}
	}

	return 0;
}

int TessAtom_match(const TessAtom *T, const Atom *A)
{
	// Return true if A matches T...

	switch(T->code)
	{
	case -1:
	case 0:
		// An exact match on both atom name and residue name.

		if(!TessAtom_matchName(T,A)) return 0;
		return TessAtom_matchResName(T,A);

	case 1:
		// An exact match on residue name and any non-carbon
		// side-chain atom.

		if(TessAtom_isCarbon(A)) return 0;
		if(TessAtom_isHydrogen(A)) return 0; // Riziotis - Only match heavy atoms
		if(TessAtom_isMainChain(A)) return 0;
		return TessAtom_matchResName(T,A);

	case 2:
		// Any non-carbon atom in the list of residues given...

		if(TessAtom_isCarbon(A)) return 0;
		if(TessAtom_isHydrogen(A)) return 0; // Riziotis - Only match heavy atoms
		return TessAtom_matchResName(T,A);

	case 3:
		// Atom type specified is residue(s) given

		if(!TessAtom_matchType(T,A)) return 0;
		if(TessAtom_isHydrogen(A)) return 0; // Riziotis - Only match heavy atoms
		return TessAtom_matchResName(T,A);

	case 4:
		// Non-carbon main-chain in given residue(s)

		if(TessAtom_isCarbon(A)) return 0;
		if(TessAtom_isHydrogen(A)) return 0; // Riziotis - Only match heavy atoms
		if(!TessAtom_isMainChain(A)) return 0;
		return TessAtom_matchResName(T,A);

	case 5:
		// Any main-chain atom in the given residue(s)

		if(!TessAtom_isMainChain(A)) return 0;
		if(TessAtom_isHydrogen(A)) return 0; // Riziotis - Only match heavy atoms
		return TessAtom_matchResName(T,A);

	case 6:
		// Any side-chain atom in the given residue(s)

		if(TessAtom_isMainChain(A)) return 0;
		if(TessAtom_isHydrogen(A)) return 0; // Riziotis - Only match heavy atoms
		return TessAtom_matchResName(T,A);

	case 7:
		// Any atom (in the specified resdiue)

		if(TessAtom_isHydrogen(A)) return 0; // Riziotis - Only match heavy atoms
		return TessAtom_matchResName(T,A);

	//Riziotis options
	case 8:
		//Any atom in the same position in the given residue(s)
		
		if(TessAtom_isHydrogen(A)) return 0; // Riziotis - Only match heavy atoms
		if(!TessAtom_isInSamePosition(T,A)) return 0;
		return TessAtom_matchResName(T,A);


	// Gail's options follow (no residue specificity)

	case 100:
		// name match

		return TessAtom_matchName(T,A);

	case 101:
		// any non-carbon side-chain atom.

		if(TessAtom_isCarbon(A)) return 0;
		if(TessAtom_isMainChain(A)) return 0;
		return 1;

	case 102:
		// Any non-carbon atom...

		if(TessAtom_isCarbon(A)) return 0;
		return 1;

	case 103:
		// Atom type specified (no residue considered)

		if(!TessAtom_matchType(T,A)) return 0;
		return 1;

	case 104:
		// Non-carbon main-chain

		if(TessAtom_isCarbon(A)) return 0;
		if(!TessAtom_isMainChain(A)) return 0;

	case 105:
		// Any main-chain atom in the given residue(s)

		return TessAtom_isMainChain(A);

	case 106:
		// Any side-chain atom in the given residue(s)

		if(TessAtom_isMainChain(A)) return 0;
		return 1;

	case 107:
		// Any atom AT ALL!!

		return 1;

	default:
		fprintf(stderr,"TessAtom_match: unknown match code (%i)\n",T->code);
		exit(-1);
	}
}

// ==================================================================

