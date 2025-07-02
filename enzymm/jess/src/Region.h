// ==================================================================
// Region.h
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Declaration of abstract type Region.
// ==================================================================

#ifndef REGION_H
#define REGION_H

// ==================================================================
// Forward declarations
// ==================================================================
// Region				An abstraction of a region
// ==================================================================

typedef struct _Region Region;

// ==================================================================
// type Region
// ==================================================================
// intersectionQ		The intersection oracle
// inclusionQ			The inclusion oracle
// free(R)				Destroy region (free memory)
// ==================================================================

struct _Region
{
	int (*intersectionQ)(Region*,double*,double*,int);
	int (*inclusionQ)(Region*,double*,int);
	void (*free)(Region*);
};

// ==================================================================
// Methods of type Region
// ==================================================================
// volume(R,p,a,b,d)	Compute volume of R to precision p 
// ==================================================================

double Region_volume(Region*,double,double*,double*,int);

// ==================================================================

#endif

