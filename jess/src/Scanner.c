// ==================================================================
// Scanner.c
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Implementation of type Scanner (the main Jess query object).
// ==================================================================

#include "Scanner.h"
#include "KdTree.h"
#include "Region.h"
#include "Annulus.h"
#include "Join.h"
#include <stdlib.h>
#include <string.h>

// ==================================================================
// Local type CandidateSet
// ==================================================================
// count				Number of atoms in the set
// atom[k]				Points to ATOM record for kth candidate
// coord[k]				Points to coordinates for kth candidate
// ==================================================================

typedef struct _CandidateSet CandidateSet;

struct _CandidateSet
{
	int count;
	Atom **atom;
	double **coord;
};

// ==================================================================
// Declaration of methods of local type CandidateSet
// ==================================================================
// create(M,T,k)		Create from molecule M, atom k of T
// free(S)				Free candidate set
// ==================================================================

static CandidateSet *CandidateSet_create(Molecule*,Template*,int);
static void CandidateSet_free(CandidateSet*);

// ==================================================================
// type Scanner
//===================================================================
// template				The template object
// set[k]				Set of candidates for atom k
// tree[k]				Tree of candidate positions for atom k
// query[k]				Current query state for tree k
// index[k]				Index of result[k] in set[k].
// result[k]			kth atom of current result set
// region[i]			Temporary region pointer
// count				= template->count(template)
// threshold			The global distance cutoff  
// max_total_threshold		The maximum value the distance cutoff can take
// 				after adding the global and single-residue
// 				distance cutoff
//===================================================================

struct _Scanner
{
	Template *template;
	CandidateSet **set;
	KdTree **tree;
	KdTreeQuery **query;
	int *index;
	Atom **atom;
	Region **region;
	int count;
	double threshold;
	double max_total_threshold;
};

// ==================================================================
// Methods of type Scanner
// ==================================================================

Scanner *Scanner_create(Molecule *M, Template *T,double r, double s)
{
	Scanner *S;
	int k,n=T->count(T);
	int m;

	S=(Scanner*)calloc(1,sizeof(Scanner));
	S->set=(CandidateSet**)calloc(n,sizeof(CandidateSet*));
	S->tree=(KdTree**)calloc(n,sizeof(KdTree*));
	S->query=(KdTreeQuery**)calloc(n,sizeof(KdTreeQuery*));
	S->index=(int*)calloc(n,sizeof(int));
	S->atom=(Atom**)calloc(n,sizeof(Atom*));
	S->region=(Region**)calloc(n,sizeof(Region*));

	S->template=T;
	S->threshold=r;
	S->max_total_threshold=s;
	S->count=n;

	for(k=0; k<n; k++)
	{
		S->index[k]=-1;
		S->set[k]=CandidateSet_create(M,T,k);

		if(S->set[k]->count==0)
		{
			Scanner_free(S);
			return NULL;
		}

		S->tree[k]=KdTree_create(S->set[k]->coord,S->set[k]->count,3);
	}

	if(S->count>0 && S->set[0]->count>0)
	{
		S->index[0]=0;
		S->atom[0]=S->set[0]->atom[0];
	}

	return S;
}

void Scanner_free(Scanner *S)
{
	int k,n;

	if(S)
	{
		n = S->template->count(S->template);

		for(k=0; k<n; k++)
		{
			if(S->set && S->set[k]) CandidateSet_free(S->set[k]);
			if(S->tree && S->tree[k]) KdTree_free(S->tree[k]);
			if(S->query && S->query[k]) KdTreeQuery_free(S->query[k]);
		}

		if(S->set) free(S->set);
		if(S->query) free(S->query);
		if(S->tree) free(S->tree);
		if(S->atom) free(S->atom);
		if(S->index) free(S->index);
		if(S->region) free(S->region);

		free(S);
	}
}

Atom **Scanner_next(Scanner *S, int ignore_chain)
{
	int j,k;
	double min,max;
	Region *J;

	double dynamic_threshold = S->threshold;

	k=S->count-1;

	// Attempt to find the next query result.

	while(k>=0)
	{
		// If k==S->count, we have a hit!

		if(k==S->count) break;

		// If k==0 we must find the next left-most
		// atom of the query...

		if(k==0)
		{
			S->index[0]++;

			if(S->index[0]>=S->set[0]->count)
			{
				// End of query...

				k=-1;
			}
			else
			{
				S->atom[0]=S->set[0]->atom[S->index[0]];
				k++;
			}

			continue;
		}

		// So k>0. If there is an active query for this
		// set then query it now...

		if(S->query[k])
		{
			S->index[k]=KdTreeQuery_next(S->query[k]);
			if(S->index[k]<0)
			{
				// The query ended. So we need to destroy
				// this query, then drop down a level...

				KdTreeQuery_free(S->query[k]);
				S->query[k]=NULL;
				S->atom[k]=NULL;
				k--;
			}
			else
			{
				// The query was successful(?) Remember the
				// atom, check n-ary constraints and continue
				// up...

				S->atom[k]=S->set[k]->atom[S->index[k]];
				if(S->template->check(S->template,S->atom,k+1,ignore_chain))
				{
					k++;
				}
			}

			continue;
		}

		// There is no active query for set k. If there is
		// no active query result for set k-1 then we need
		// to drop down again...

		if(S->index[k-1]<0)
		{
			k--;
			continue;
		}

		// So, there is an active query result at k-1 and
		// no active query at k; create a new query at
		// index k and try again (with the same k)

		for(j=0; j<k; j++)
		{
			S->template->range(S->template,j,k,&min,&max);

			dynamic_threshold = S->threshold + S->template->distWeight(S->template, j) + S->template->distWeight(S->template, k);
			// Limit threshold to a hard cutoff so execution does not suffer
			if(dynamic_threshold > S->max_total_threshold){
				dynamic_threshold = S->max_total_threshold;
			}
			min -= dynamic_threshold;
			max += dynamic_threshold;
			if(min<0.5) min=0.5;

			S->region[j]=Annulus_create(S->atom[j]->x,min,max,3);
		}

		J = Join_create(S->region,k,innerJoin);
		S->query[k]=KdTree_query(S->tree[k],J);
	}

	// If k<0 there is no more!

	if(k<0) return NULL;

	// Otherwise, the atoms are listed in S->atom.

	return S->atom;
}

// ==================================================================
// Methods of local type CandidateSet
// ==================================================================

static CandidateSet *CandidateSet_create(Molecule *M, Template *T, int k)
{
	CandidateSet *S;
	Atom *A;
	int n = Molecule_count(M);
	int m;

	S = (CandidateSet*)calloc(1,sizeof(CandidateSet));
	S->atom=(Atom**)calloc(n,sizeof(Atom*));

	for(m=0; m<n; m++)
	{
		A = (Atom*)Molecule_atom(M,m);
		if(T->match(T,k,A))
		{
			S->atom[S->count]=A;
			S->count++;
		}
	}

	S->atom=(Atom**)realloc(S->atom,sizeof(Atom*)*S->count);
	S->coord=(double**)calloc(S->count,sizeof(double*));

	for(m=0; m<S->count; m++)
	{
		S->coord[m]=S->atom[m]->x;
	}

	return S;
}

static void CandidateSet_free(CandidateSet *S)
{
	if(S)
	{
		if(S->atom) free(S->atom);
		if(S->coord) free(S->coord);
		free(S);
	}
}

// ==================================================================

