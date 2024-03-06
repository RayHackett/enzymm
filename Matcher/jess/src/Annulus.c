// ==================================================================
// Annulus.c
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Implementation of the Annulus region and its oracles
// ==================================================================

#include "Annulus.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

// ==================================================================
// Local "functions"
// ==================================================================

#define min(x,y) (x<y ? x:y)
#define max(x,y) (x>y ? x:y)

// ==================================================================
// type Annulus
// ==================================================================
// centre				The center of the annulus
// min,max				Limits of the radii
// dim					The dimension of the space
// ==================================================================

typedef struct _Annulus
{
	double *centre;
	double min;
	double max;
	int dim;
}
Annulus;

// ==================================================================
// Oracles for type Annulus
// ==================================================================

static int Annulus_po(Region *vA, double *x, int d)
{
	Annulus *A=(Annulus*)&vA[1];
	double tmp,sum;
	int i;

	// Does x lie within annulus A?

	if(A->dim!=d) return 0;

	for(sum=0.0,i=0; i<d; i++)
	{
		tmp = A->centre[i]-x[i];
		sum += tmp*tmp;
	}

	return sum<A->min || sum>A->max ? 0:1;
}

static int Annulus_ro(Region *vA, double *minBox, double *maxBox, int d)
{
	Annulus *A=(Annulus*)&vA[1];
	double minSum;
	double maxSum;
 	double t1,t2;
	double t3,t4;
	int i;

	if(d!=A->dim) return 0;

	// Does the box region [minBox,maxBox] intersect the annulus A?

	minSum=0.0;
	maxSum=0.0;
	for(i=0; i<d; i++)
	{
		t1 = A->centre[i]-minBox[i];
		t2 = A->centre[i]-maxBox[i];
		t1 *= t1;
		t2 *= t2;

		if(minBox[i]>A->centre[i] || maxBox[i]<A->centre[i])
		{
			minSum += min(t1,t2);
		}

		maxSum += max(t1,t2);
	}

	return minSum>A->max || maxSum<A->min ? 0:1;
}

// ==================================================================
// Methods of for regions of type Annulus
// ==================================================================

Region *Annulus_create(double *u, double a, double b, int d)
{
	Annulus *A;
	Region *R;
	int rq;
	double tmp;

	rq = sizeof(Region)+sizeof(Annulus)+d*sizeof(double);
	R = (Region*)calloc(1,rq);
	A = (Annulus*)&R[1];
	R->intersectionQ=Annulus_ro;
	R->inclusionQ=Annulus_po;
	R->free=Annulus_free;

	if(b<a)
	{
		tmp=b;
		b=a;
		a=tmp;
	}

	if(a<0.0) a=0.0;
	if(b<0.0) b=0.0;

	A->centre=(double*)&A[1];
	memcpy(A->centre,u,sizeof(double)*d);
	A->min=a*a;
	A->max=b*b;
	A->dim=d;

	return R;
}

void Annulus_free(Region *R)
{
	if(R) free(R);
}

// ==================================================================

