// ==================================================================
// TessTemplate.c
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Implementation of TessTemplate creation and oracles.
// ==================================================================

#include "TessTemplate.h"
#include "TessAtom.h"
#include "Annulus.h"
#include "Join.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

// ==================================================================
// Forward declaration of local types
// ==================================================================
// Node					A node in a linked list of TessAtoms
// TessTemplate			The implementation of type TessTemplate
// ==================================================================

typedef struct _Node Node;
typedef struct _TessTemplate TessTemplate;

// ==================================================================
// Local type Node
// ==================================================================

struct _Node
{
	TessAtom *atom;
	Node *succ;
};

// ==================================================================
// Local type TessTemplate
// ==================================================================
// count				Number of atoms in the template
// atom[k]				Ptr to TessAtom k
// distance[i][j]		Distance between atoms i and j
// symbol				The name of the template
// dim					The dimension of the template (# residues)
// ==================================================================

struct _TessTemplate
{
	int count;
	TessAtom **atom;
	double **distance;
	char *symbol;
	int dim;
};

// ==================================================================
// Oracles of type TessTemplate
// ==================================================================

static int TessTemplate_count(const Template *T)
{
	const TessTemplate *J = (const TessTemplate*)&T[1];
	return J->count;
}

static int TessTemplate_match(const Template *T,int k,const Atom *A)
{
	const TessTemplate *J = (const TessTemplate*)&T[1];
	return TessAtom_match(J->atom[k],A);
}

static int TessTemplate_range(const Template *T,int i,int j,double *a,double *b)
{
	const TessTemplate *J = (const TessTemplate*)&T[1];
	*a = *b = J->distance[i][j];
	return 1;
}

static const double *TessTemplate_position(const Template *T, int k)
{
	const TessTemplate *J=(const TessTemplate*)&T[1];
	return TessAtom_position(J->atom[k]);
}

static double TessTemplate_distWeight(const Template *T, int k)
{
	const TessTemplate *J=(const TessTemplate*)&T[1];
	return TessAtom_distWeight(J->atom[k]);
}

static int TessTemplate_check(const Template *T, Atom **A, int k, int ignore_chain)
{
	const TessTemplate *J = (const TessTemplate*)&T[1];
	int i;
	int a,b,c,d,e,f;

	for(i=0; i<k-1; i++)
	{
		// Compare chain ids

		//Riziotis edit
		c = A[i]->chainID2-A[k-1]->chainID2;
		d = TessAtom_chainID2(J->atom[i])-TessAtom_chainID2(J->atom[k-1]);
		//c = A[i]->chainID-A[k-1]->chainID;
		//d = TessAtom_chainID(J->atom[i])-TessAtom_chainID(J->atom[k-1]);

		if(ignore_chain==1)
		{
			c=0;
			d=0;
		}

		if(c==0 && d!=0) return 0;
		if(c!=0 && d==0) return 0;

		if(c!=0)
		{
			continue;
		}

		// Compare residue sequence numbers

		a = A[i]->resSeq-A[k-1]->resSeq;
		b = TessAtom_resSeq(J->atom[i])-TessAtom_resSeq(J->atom[k-1]);

		if(a==0 && b!=0) return 0;
		if(a!=0 && b==0) return 0;
	}

	return 1;
}

static const char *TessTemplate_name(const Template *T)
{
	const TessTemplate *J = (const TessTemplate*)&T[1];
	return J->symbol;
}

static double TessTemplate_logE(const Template *T,double rmsd, int n)
{
	const TessTemplate *J = (const TessTemplate*)&T[1];
	const double logA = -8.5;
	const double alpha = 2.5;
	const double beta = -0.7;

	// Approximate number of hits per molecule atom is:
	// A * rmsd^alpha * J->dim^beta (empirically). This gives
	// us log(expected number of hits in molecule of size n)

	return logA + alpha*rmsd + beta*(double)J->dim + log((double)n);
}

// ==================================================================
// Private methods of type TessTemplate
// ==================================================================

static void TessTemplate_free(Template *T)
{
	TessTemplate *J;
	int i;

	if(T)
	{
		J=(TessTemplate*)&T[1];
		if(J->symbol) free(J->symbol);

		for(i=0; i<J->count; i++)
		{
			TessAtom_free(J->atom[i]);
		}

		free(T);
	}
}

// ==================================================================
// Creation of type TessTemplate
// ==================================================================

Template *TessTemplate_create(FILE *file,const char *sym)
{
	Template *T;
	TessTemplate *J;
	TessAtom *A;
	Node *head,*n;
	const double *x;
	const double *y;
	double tmp;
	int rq,i,j,k,count;
	char buf[0x100];

	// head is the head of a list of TessAtoms which
	// we create by parsing file. count is the number
	// of atom templates encountered.

	head=NULL;
	count=0;

	// Loop through the entire stream, line by line.

	while(fgets(buf,0x100,file))
	{
		if(strncmp(buf,"ATOM",4)==0 || strncmp(buf,"HETATM",6)==0)
		{
			// Parse a single TESS/Jess atom record...

			A = TessAtom_create(buf);

			// If some errors occur parsing the template
			// free all memory and return NULL.

			if(!A)
			{
				while(head)
				{
					n=head->succ;
					TessAtom_free(head->atom);
					free(head);
					head=n;
				}

				return NULL;
			}

			// No errors thus far - add the atom template
			// to the list we're creating...

			n = (Node*)calloc(1,sizeof(Node));
			n->succ=head;
			head=n;
			n->atom=A;
			count++;
		}
	}

	if(count==0) return NULL;

	// Set up the memory for the template. This is a
	// Template, followed by a TessTemplate followed
	// by various arrays needed for the TessTemplate.
	// To allow easy destruction it makes sense to
	// allocate it all in one big chunk, then set
	// pointers later...

	rq = sizeof(Template)+sizeof(TessTemplate);
	rq += count*sizeof(TessAtom*);
	rq += count*count*sizeof(double);
	rq += count*sizeof(double*);

	T = (Template*)calloc(1,rq);
	J = (TessTemplate*)&T[1];
	J->atom=(TessAtom**)&J[1];
	J->distance=(double**)&J->atom[count];

	J->distance[0]=(double*)&J->distance[count];
	for(i=1; i<count; i++)
	{
		J->distance[i]=(double*)&J->distance[i-1][count];
	}

	// Set up the method pointers in the Template
	// structure at the head of the object...

	T->free=TessTemplate_free;
	T->match=TessTemplate_match;
	T->position=TessTemplate_position;
	T->count=TessTemplate_count;
	T->range=TessTemplate_range;
	T->check=TessTemplate_check;
	T->name=TessTemplate_name;
	T->logE=TessTemplate_logE;
	T->distWeight=TessTemplate_distWeight;

	// Set up the data fields

	J->symbol=strdup(sym);
	J->count=count;

	// Create pointers to the atoms found in the
	// template to allow access by atom index rather
	// than iteration through a list.

	for(i=0; i<count; i++)
	{
		n=head->succ;
		J->atom[count-i-1]=head->atom;
		free(head);
		head=n;
	}

	// Compute the distances in the template.

	for(i=0; i<count; i++)
	{
		J->distance[i][i]=0.0;
		x=TessAtom_position(J->atom[i]);

		for(j=i+1; j<count; j++)
		{
			J->distance[i][j]=0.0;
			y=TessAtom_position(J->atom[j]);

			for(k=0; k<3; k++)
			{
				tmp = x[k]-y[k];
				J->distance[i][j] += tmp*tmp;
			}

			J->distance[i][j]=sqrt(J->distance[i][j]);
			J->distance[j][i]=J->distance[i][j];
		}
	}

	// Compute the dimension of the template...

	for(J->dim=0,i=0; i<count; i++)
	{
		// Have we already seen this residue sequence number?

		for(j=0; j<i; j++)
		{
			if(TessAtom_resSeq(J->atom[i])==TessAtom_resSeq(J->atom[j]))
			{
				j=i+10;
			}
		}

		if(j<=i)
		{
			J->dim++;
		}
	}

	return T;
}

// ==================================================================
