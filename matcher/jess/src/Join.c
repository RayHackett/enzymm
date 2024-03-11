// ==================================================================
// Join.c
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Implementation of the Join region code.
// ==================================================================

#include "Join.h"
#include <stdlib.h>
#include <string.h>

// ==================================================================
// type Join
// ==================================================================
// count				The number of regions in the join
// R[k]					The kth region in the join
// ==================================================================

typedef struct _Join Join;

struct _Join
{
	int count;
	Region *R[0];
};

// ==================================================================
// The oracles
// ==================================================================

static int Join_oro(Region *R,double *min,double *max,int dim)
{
	Join *J=(Join*)&R[1];
	int k;

	for(k=0; k<J->count; k++)
	{
		if(J->R[k]->intersectionQ(J->R[k],min,max,dim))
		{
			return 1;
		}
	}

	return 0;
}

static int Join_iro(Region *R,double *min,double *max,int dim)
{
	Join *J=(Join*)&R[1];
	int k;

	for(k=0; k<J->count; k++)
	{
		if(!(J->R[k]->intersectionQ(J->R[k],min,max,dim)))
		{
			return 0;
		}
	}

	return 1;
}

static int Join_opo(Region *R,double *x,int dim)
{
	Join *J=(Join*)&R[1];
	int k;

	for(k=0; k<J->count; k++)
	{
		if(J->R[k]->inclusionQ(J->R[k],x,dim))
		{
			return 1;
		}
	}

	return 0;
}

static int Join_ipo(Region *R,double *x,int dim)
{
	Join *J=(Join*)&R[1];
	int k;

	for(k=0; k<J->count; k++)
	{
		if(!(J->R[k]->inclusionQ(J->R[k],x,dim)))
		{
			return 0;
		}
	}

	return 1;
}

// ==================================================================
// Construction and destruction
// ==================================================================

Region *Join_create(Region **S,int count,JoinType type)
{
	Region *R;
	Join *J;
	int rq;

	rq = sizeof(Region)+sizeof(Join)+sizeof(Region*)*count;
	R = (Region*)calloc(1,rq);
	J = (Join*)&R[1];
	memcpy(J->R,S,sizeof(Region*)*count);

	R->free=Join_free;
	J->count=count;

	if(type==innerJoin)
	{
		R->intersectionQ = Join_iro;
		R->inclusionQ = Join_ipo;
	}
	else // if type==outerJoin
	{
		R->intersectionQ = Join_oro;
		R->inclusionQ = Join_opo;
	}

	return R;
}

void Join_free(Region *R)
{
	Join *J;
	int k;

	if(R)
	{
		J = (Join*)&R[1];

		for(k=0; k<J->count; k++)
		{
			if(J->R[k]) J->R[k]->free(J->R[k]);
		}

		free(R);
	}
}

// ==================================================================


