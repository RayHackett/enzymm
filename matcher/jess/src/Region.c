// ==================================================================
// Region.c
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Implementation of methods of abstract type Region
// ==================================================================

#include "Region.h"
#include <math.h>

// ==================================================================
// volumeR: recursive implementation
// ==================================================================
// R				The region
// p				The precision
// a,b				Limits of box
// d				Dimension of space
// r				Current splitting axis
// v				Current volume of box
// ==================================================================

static double Region_volumeR(
	Region *R,
	double p,
	double *a,
	double *b,
	int d,
	int r,
	double v
	)
{
	double split;
	double tmp;
	double sum;

	// First, if the volume of the box is small
	// enough, then just test for intersection

	if(v<=p)
	{
		if(R->intersectionQ(R,a,b,d))
		{
			return v;
		}
		else
		{
			return 0.0;
		}
	}

	// So we need to split the box into 2 at
	// coordinate axis r...

	split = (a[r]+b[r])/2.0;
	tmp=a[r];
	a[r]=split;

	sum = 0.0;

	if(R->intersectionQ(R,a,b,d))
		sum += Region_volumeR(R,p,a,b,d,(r+1)%d,v/2.0);

	a[r]=tmp;
	tmp=b[r];
	b[r]=split;

	if(R->intersectionQ(R,a,b,d))
		sum += Region_volumeR(R,p,a,b,d,(r+1)%d,v/2.0);

	b[r]=tmp;
	return sum;
}

// ==================================================================
// volume: the main interface
// ==================================================================

double Region_volume(Region *R,double p,double *a,double *b,int d)
{
	double v;
	int i;

	v=1.0;
	for(i=0; i<d; i++)
	{
		v *= a[i]-b[i];
	}
	v = fabs(v);

	return Region_volumeR(R,p,a,b,d,0,v);
}

// ==================================================================
// Test code
// ==================================================================
/*
#include "Annulus.h"
#include "Join.h"

int main(void)
{
	Region *R;
	Region *A,*B;
	Region *L[2];
	double u[3];
	double v[3];
	double c[3];
	double d[3];

	u[0]=u[1]=u[2]=0.0;
	v[0]=v[1]=0.0;

	for(v[2]=-2.0; v[2]<2.0; v[2]+=0.1)
	{

		c[0]=c[1]=c[2]=-50.0;
		d[0]=d[1]=d[2]=50.0;

		A = Annulus_create(u,0.8,1.0,3);
		B = Annulus_create(v,0.8,1.0,3);

		L[0]=A;
		L[1]=B;
		R = Join_create(L,2,innerJoin);

		printf("%6.2f %8.3g\n",v[2],Region_volume(R,1e-5,c,d,3));

		R->free(R);
	}

	return 0;
}
*/
// ==================================================================
