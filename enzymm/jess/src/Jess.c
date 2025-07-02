// ================================================================== Jess.c
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Implementation of types Jess and JessQuery.
// ==================================================================

#include "Jess.h"
#include "Molecule.h"
#include "Scanner.h"
#include "TessTemplate.h"
#include "Super.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// ==================================================================
// Forward declarations of local types
// ==================================================================
// Node					A linked list of templates
// ==================================================================

typedef struct _Node Node;

// ==================================================================
// type Jess
// ==================================================================
// head					The head of the list of nodes
// ==================================================================

struct _Jess
{
	Node *head;
};

// ==================================================================
// type JessQuery
// ==================================================================
// node					The current node
// scanner				The current scanner
// super				The current superposition
// reverseQ				True if superposition is reversed
// molecule				The molecule being scanned
// atoms				Array of Atoms which are hit
// threshold			The distance threshold
// ==================================================================

struct _JessQuery
{
	Node *node;
	Scanner *scanner;
	Superposition *super;
	int reverseQ;
	Molecule *molecule;
	Atom **atoms;
	double threshold;
	double max_total_threshold;
};

// ==================================================================
// type Node
// ==================================================================
// template				The template at this node
// next					The next node in the list
// ==================================================================

struct _Node
{
	Template *template;
	Node *next;
};

// ==================================================================
// Methods of type Jess
// ==================================================================

Jess *Jess_create(void)
{
	return (Jess*)calloc(1,sizeof(Jess));
}

void Jess_free(Jess *J)
{
	Node *n;
	Template *T;

	if(J)
	{
		while(J->head)
		{
			n=J->head->next;
			T=J->head->template;
			if(T) T->free(T);
			J->head=n;
		}
	}
}

void Jess_addTemplate(Jess *J, Template *T)
{
	Node *n;

	n=(Node*)calloc(1,sizeof(Node));
	n->template=T;
	n->next=J->head;
	J->head=n;
}

JessQuery *Jess_query(Jess *J, Molecule *M,double t,double s)
{
	JessQuery *Q;

	Q = (JessQuery*)calloc(1,sizeof(JessQuery));
	Q->node=J->head;
	Q->molecule=M;
	Q->threshold=t;
	Q->max_total_threshold=s;

	return Q;
}

// ==================================================================
// Methods of type JessQuery
// ==================================================================

void JessQuery_free(JessQuery *Q)
{
	if(Q)
	{
		Scanner_free(Q->scanner);
		Superposition_free(Q->super);
		free(Q);
	}
}

Template *JessQuery_template(JessQuery *Q)
{
	if(!Q->node) return NULL;
	return Q->node->template;
}

const Molecule *JessQuery_molecule(JessQuery *Q)
{
	return Q->molecule;
}

Atom **JessQuery_atoms(JessQuery *Q)
{
	if(!Q->atoms) return NULL;
	return Q->atoms;
}

Superposition *JessQuery_superposition(JessQuery *Q)
{
	int i;
	int count;
	Template *T;
	Atom **A;

	if(Q->super) return Q->super;

	A = Q->atoms;
	T = Q->node->template;
	count = T->count(T);
	Q->super=Superposition_create();

	for(i=0; i<count; i++)
	{
		Superposition_align(Q->super,A[i]->x,T->position(T,i));
	}

	return Q->super;
}

int JessQuery_next(JessQuery *Q, int ignore_chain)
{
	Template *T;
	Atom **A;
	Superposition *S;

	while(Q->node)
	{
		Superposition_free(Q->super);
		Q->super=NULL;

		if(!Q->scanner)
		{
			Q->scanner=Scanner_create(
				Q->molecule,
				Q->node->template,
				Q->threshold,
				Q->max_total_threshold
				);

			if(!Q->scanner)
			{
				Q->node=Q->node->next;
				continue;
			}
		}

		if((A=Scanner_next(Q->scanner, ignore_chain)))
		{
			Q->atoms=A;
			T = Q->node->template;

			return 1;
		}

		Scanner_free(Q->scanner);
		Q->scanner=NULL;

		Superposition_free(Q->super);
		Q->super=NULL;

		Q->atoms=NULL;
		Q->node=Q->node->next;
	}

	// All done...

	return 0;
}

// ==================================================================

