// ==================================================================
// Molecule.c
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Implementation of type Molecule.
// ==================================================================

#include "Molecule.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

// ==================================================================
// Local type Node
// ==================================================================
// next					Next node in the list of all nodes
// atom					Atom record stored at this node
// ==================================================================

typedef struct _Node Node;

struct _Node
{
	Node *next;
	Atom *atom;
};

// ==================================================================
// type Molecule
// ==================================================================
// count				Number of atoms in the molecule
// id					The molecule PDB code (if found)
// atom[k]				Pointer to kth atom in the molecule
// ==================================================================

struct _Molecule
{
	int count;
	char id[5];
	Atom *atom[0];
};

// ==================================================================
// Methods of type Molecule;
// ==================================================================

Molecule *Molecule_create(FILE *file, int ignore_endmdl, float conservation_cutoff)
{
	Node *head=NULL;
	Node *N;
	Atom A;
	Atom *pA;
	Molecule *M;
	char buf[0x100];
	char pdb[5];
	int count=0;

	pdb[0]=0;

	// Loop through the file and read all of the ATOM
	// records. Discard all altLoc atoms and those which
	// occur after the end of the first MODEL.

	memset(buf,0,0x100);
	while(fgets(buf,0x100,file))
	{
		// If we want to include all models (for instance in biounit PDB structures)
		if(ignore_endmdl==0 && strncmp(buf,"ENDMDL",6)==0){
			break;
		}
    
		// Get the PDB code if possible
		if(strncmp(buf,"HEADER",6)==0)
		{
			strncpy(pdb,&buf[62],4);
			pdb[4]=0;
		}

		// Parse an atom record if possible...
		if(Atom_parse(&A,buf)) //Riziotis edit: we want the altLoc atoms
		//if(Atom_parse(&A,buf) && isspace(A.altLoc))
		{
			// We got one! Create a new node in
			// the list...

			pA = (Atom*)calloc(1,sizeof(Atom));
			memcpy(pA,&A,sizeof(Atom));

			// If there is a user-defined conservation score cutoff in 
			// the temperature factor field, ignore this atom if its
			// score is lower than the cutoff
			
			if(conservation_cutoff>0 && pA->tempFactor<conservation_cutoff)
			{
				continue;
			}

			N = (Node*)calloc(1,sizeof(Node));
			N->next=head;
			N->atom=pA;
			head=N;
			count++;
		}
		
		memset(buf,0,0x100);
	}

	// Right, if count>0 we got some atoms. Otherwise
	// return NULL now!

	if(count<=0) return NULL;

	// Create the molecule...

	M = (Molecule*)calloc(1,sizeof(Molecule)+count*sizeof(Atom*));
	M->count=count;
	strcpy(M->id,pdb);

	// Loop through the list and add all the atoms to it.
	// Remember that we added them all backwards!

	while(count-->0)
	{
		M->atom[count]=head->atom;
		N=head->next;
		free(head);
		head=N;
	}

	// Check memory leaks??

	return M;
}

void Molecule_free(Molecule *M)
{
	int k;

	if(M)
	{
		for(k=0; k<M->count; k++)
		{
			if(M->atom[k]) free(M->atom[k]);
		}

		free(M);
	}
}

int Molecule_count(const Molecule *M)
{
	return M->count;
}

const Atom *Molecule_atom(const Molecule *M, int k)
{
	return !M || k<0 || k>=M->count ? NULL:M->atom[k];
}

const char *Molecule_id(const Molecule *M)
{
	if(strlen(M->id)==4 && strcmp(M->id, "    ")!=0) return M->id;
	return NULL;
}

// ==================================================================

