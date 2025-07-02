// ==================================================================
// Join.h
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Declaration of region type Join. This joins two regions as a union
// or intersection. Caveat: the intersection oracle is sub-optimal
// but I don't know how to fix it!
// ==================================================================

#ifndef JOIN_H
#define JOIN_H

#include "Region.h"

// ==================================================================
// Join types
// ==================================================================
// innerJoin				An intersection of two regions
// outerJoin				The union of two regions
// ==================================================================

typedef enum {innerJoin,outerJoin} JoinType;

// ==================================================================
// Construction methods for region type Join
// ==================================================================
// create(R,n,type)			Create inner/outer join on R[0...n-1]
// free(J)					Frees join AND nested regions (J->free)
// ==================================================================

extern Region *Join_create(Region**,int,JoinType);
extern void Join_free(Region*);

// ==================================================================

#endif

