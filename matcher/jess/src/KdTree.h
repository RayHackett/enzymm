// ==================================================================
// KdTree.h
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Declaration of type KdTree, related types and methods.
// ==================================================================

#ifndef KDTREE_H
#define KDTREE_H

#include "Region.h"

// ==================================================================
// Forward declarations
// ==================================================================
// KdTree					A kd-tree
// KdTreeQuery				Result of query on a KdTree
// ==================================================================

typedef struct _KdTree KdTree;
typedef struct _KdTreeQuery KdTreeQuery;

// ==================================================================
// Methods of type KdTree
// ==================================================================
// create(u,n,k)			Create kd-tree on u[0],...,u[n-1]
// free(K)					Free the kd-tree K
// query(K,R)				Initialise a query object (see code)
// ==================================================================

extern KdTree *KdTree_create(double**,int,int);
extern void KdTree_free(KdTree*);
extern KdTreeQuery *KdTree_query(KdTree*,Region*);

// ==================================================================
// Methods of type KdTreeQuery
// ==================================================================
// free(Q)					Free memory for query object AND region
// next(Q)					Return next result in Q (-1 if no more)
// ==================================================================

extern void KdTreeQuery_free(KdTreeQuery*);
extern int KdTreeQuery_next(KdTreeQuery*);

// ==================================================================

#endif

