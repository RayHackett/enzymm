// ==================================================================
// Super.c
// Copyright (c) Jonathan Barker, 2001
// ==================================================================
// Implementation of type Superposition
// ==================================================================

#include "Super.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

// ==================================================================
// Local type Node
// ==================================================================
// data[]			The 2 vectors stored at this node
// succ				The next pair.
// ==================================================================

typedef struct _Node Node;

struct _Node
{
	double data[2][3];
	Node *succ;
};

// ==================================================================
// Methods of local type Node
// ==================================================================

static Node *Node_create(const double*, const double*);
static void Node_free(Node*);

// ==================================================================
// type Superposition
// ==================================================================
// upToDate			True if rmsd up to date.
// rmsd				The current rmsd
// rmsd100			The current rmsd100
// count			Number of vector pairs added
// head				The head of the list of pairs
// centre[k]		The centroid of the kth set (k=0,1)
// rotation			The rotation matrix used
// ==================================================================

struct _Superposition
{
	int upToDate;
	double rmsd;
	double rmsd100;
	int count;
	double rotation[9];
	double centre[2][3];
	Node *head;
};

// ==================================================================
// Private methods of type Superposition
// ==================================================================
// compute(S)		Compute the superposition
// ==================================================================

static void Superposition_compute(Superposition*);

// ==================================================================
// Local procedures
// ==================================================================
// min(a,b)				Returns min(a,b)
// max(a,b)				Returns max(a,b)
// rotate()				Used by subroutine jacobi()
// jacobi(M,P,v)		Computes diag(v) = P^T M P (M in, P,v out)
// superpose(a,b,n,M)	Computes superposition of arrays
// ==================================================================

static double min(double,double);
static double max(double,double);
static void rotate(double*,double*,int,int);
static int jacobi(double*,double*,double*);
static double superpose(double*,double*,int,double*);

// ==================================================================
// Methods of type Superposition
// ==================================================================

Superposition *Superposition_create(void)
{
	Superposition *S;
	S = (Superposition*)calloc(1,sizeof(Superposition));
	return S;
}

void Superposition_free(Superposition *S)
{
	Node *tmp;
	if(S)
	{
		while(S->head)
		{
			tmp=S->head->succ;
			Node_free(S->head);
			S->head=tmp;
		}

		free(S);
	}
}

void Superposition_align(Superposition *S,const double *x,const double *y)
{
	Node *n;

	S->upToDate=0;
	n=Node_create(x,y);
	n->succ=S->head;
	S->head=n;
	S->count++;
}

int Superposition_count(const Superposition *S)
{
	return S ? S->count:0;
}

double Superposition_rmsd(Superposition *S)
{
	if(!S || S->count<=1) return 0.;
	if(!S->upToDate) Superposition_compute(S);
	return S->rmsd;
}

double Superposition_rmsd100(Superposition *S)
{
	if(!S || S->count<=1) return 0.;
	if(!S->upToDate) Superposition_compute(S);
	return S->rmsd100;
}

void Superposition_compute(Superposition *S)
{
	Node *n;
	double *r1;
	double *r2;
	double c1[3];
	double c2[3];
	int c,i,j,k;
	double v[3];
	double rmsd;

	double tmp,tst;

	if(S->count<=1) return;

	// Copy all the vectors into an array.
	// This is not strictly necessary, but the
	// superposition code works on arrays
	// and I can't be arsed to rewrite it.

	c=S->count;
	r1=(double*)calloc(c,sizeof(double)*3);
	r2=(double*)calloc(c,sizeof(double)*3);

	for(n=S->head,k=0; k<c; k++,n=n->succ)
	{
		for(j=0; j<3; j++)
		{
			r1[3*k+j]=n->data[0][j];
			r2[3*k+j]=n->data[1][j];
		}
	}

	// Find the centroids of both sets...

	memset(c1,0,sizeof(double)*3);
	memset(c2,0,sizeof(double)*3);

	for(i=0; i<c; i++)
	{
		for(j=0; j<3; j++)
		{
			c1[j] += r1[3*i+j];
			c2[j] += r2[3*i+j];
		}
	}

	for(j=0; j<3; j++)
	{
		c1[j]/=(double)c;
		c2[j]/=(double)c;
	}

	// Shift each set by its centroid to make each
	// centroid the origin. Forgot this at first and
	// results were FUBAR!

	for(i=0; i<c; i++)
	{
		for(j=0; j<3; j++)
		{
			r1[3*i+j] -= c1[j];
			r2[3*i+j] -= c2[j];
		}
	}

	rmsd = superpose(r1,r2,c,S->rotation);

	// Set the easy fields...

	S->rmsd=rmsd;
	S->rmsd100=rmsd/(1+0.5*log(c/100.));
	S->upToDate=1;

	memcpy(S->centre[0],c1,sizeof(double)*3);
	memcpy(S->centre[1],c2,sizeof(double)*3);

	// Free temporary memory.

	free(r1);
	free(r2);
}

const double *Superposition_centroid(Superposition *S,int k)
{
	if(k<0 || k>1) return NULL;
	if(!S->upToDate) Superposition_compute(S);
	return S->centre[k];
}

const double *Superposition_rotation(Superposition *S)
{
	if(!S->upToDate) Superposition_compute(S);
	return S->rotation;
}


// ==================================================================
// The superposition algorithm stuff
// ==================================================================

static double min(double a, double b)
{
	return a<b ? a:b;
}

static double max(double a, double b)
{
	return a>b ? a:b;
}

static const double PRECISION = 1e-12;

static void rotate(double *W, double *P, int ip, int iq)
{
	// Here we sacrifice some speed/precision in
	// favour of simplicity. I dislike heuristic
	// optimistations - and distrust them too.

	double c,s;
	double pp,qq,pq;
	double t;
	int k;

	double colWp[3];
	double colWq[3];
	double colPp[3];
	double colPq[3];

	// Here we must "rotate" element (ip,iq).
	// Rather than follow the formula given
	// in NR/C, I'll do my own - probably
	// very slow in comparison since it uses
	// a trig function. But I understand it
	// so that's all that matters...

	// Thus, W -> R^T W R and P -> P R where
	// R is an (ip,iq) rotation matrix.
	// Below, c and s are the cosine and sine
	// of the angle of rotation theta.

	// Calculate the cos and sin required
	// The new W(ip,iq) is given by...
	//
	//  W(ip,iq) -> (c*c-s*s)W(ip,iq)-sc(W(ip,ip)-W(iq,iq))
	//
	// so this is the thing we want to zero.
	// Put A = W(ip,iq), B = (W(ip,ip)-W(iq,iq)),
	// c = cos(t), s = sin(t) and we have
	//
	//  2cos(2t)A = sin(2t)B
	//
	// Now we can solve for 2t using atan2.

	t = atan2(2*W[3*ip+iq],W[3*iq+iq]-W[3*ip+ip])/(double)2;

	c = cos(t);
	s = sin(t);

	// So now we have c and s. We need to calculate
	// the new values of the iqth and ipth columns of
	// both W and of P.

	for(k=0; k<3; k++)
	{
		colWp[k] = c*W[3*k+ip] - s*W[3*k+iq];
		colWq[k] = c*W[3*k+iq] + s*W[3*k+ip];
		colPp[k] = c*P[3*k+ip] - s*P[3*k+iq];
		colPq[k] = s*P[3*k+ip] + c*P[3*k+iq];
	}

	// Some of the above are not correct. In particular
	// whenever k=ip or k=iq we need a correction for
	// the W matrix...
    // Element 	W(ip,ip)

	pp = c*c*W[3*ip+ip]+s*s*W[3*iq+iq]-(double)2*s*c*W[3*ip+iq];

	// Element W(iq,iq)

	qq = s*s*W[3*ip+ip]+c*c*W[3*iq+iq]+(double)2*s*c*W[3*ip+iq];

	// Element W(ip,iq)

	pq = (double)0; // = (c*c-s*s)*W(ip,iq) + s*c*(W(ip,ip)-W(iq,iq))  ;-)

	// Put the new values into W and into P...

	for(k=0; k<3; k++)
	{
		W[3*k+ip] = colWp[k];
		W[3*ip+k] = colWp[k];
    	W[3*k+iq] = colWq[k];
		W[3*iq+k] = colWq[k];

		P[3*k+ip] = colPp[k];
		P[3*k+iq] = colPq[k];
	}

	W[3*ip+ip] = pp;
	W[3*iq+iq] = qq;
	W[3*ip+iq] = pq;
	W[3*iq+ip] = pq;

	// We're done, and ready for the next pass
}

static int jacobi(double *M, double *P, double *v)
{
	double W[9];
	int iterationCount=0;
	int done=0;
	int i,j;
	double sum;

	// W is a copy of M to be diagonlised, P is initially the
	// identity matrix...

	memcpy(W,M,sizeof(W));
	memset(P,0,sizeof(W));
	P[0]=P[4]=P[8]=(double)1;

	// Loop until finished...

	while(!done)
	{
		iterationCount++;

		// Sum the absolute values of the off-diagonal
		// elements of the working matrix...

		sum=(double)0;

		for(i=0; i<2; i++)
		{
			for(j=i+1; j<3; j++)
			{
				sum += fabs(W[3*i+j]);
			}
		}

		// If the sum is small enough we are done. If
		// not, we need to rotate all the off-diagonal
		// elements...

		if(sum<PRECISION)
		{
			done=1;
		}
  		else
		{
			// "Rotate" all the off-diagonal elements
			// of the working matrix in conjunction with
			// P. This is the guts of the algorithm.

			for(i=0; i<2; i++)
			{
				for(j=i+1; j<3; j++)
				{
					rotate(W,P,i,j);
				}
			}
		}
	}

	// Finally, copy the diagonal values of W to the
	// vector of eigenvalues...

	for(i=0; i<3; i++)
	{
		v[i]=W[4*i];
	}


	return iterationCount;
}


static double superpose(double *a,double *b, int n, double *M)
{
	double sumA=(double)0;
	double sumB=(double)0;
	double sumE=(double)0;
	double detX;
	double X[9];
	double XX[9];
	double P[9];
	double Q[9];
	double T[9];
	double e[3];
	double rmsd;
	double factor;
	int i,j,k;
	int flag;

	// Compute sum-of squares

	for(i=0; i<3*n; i++)
	{
		sumA += a[i]*a[i];
		sumB += b[i]*b[i];
	}

	// Compute the "covariances" X

	memset(X,0,sizeof(X));
	memset(XX,0,sizeof(XX));

	for(i=0; i<n; i++)
	{
		for(j=0; j<3; j++)
		{
			for(k=0; k<3; k++)
			{
				X[3*j+k] += a[3*i+j]*b[3*i+k];
			}
		}
	}

	// And compute det(X)...

	detX =
		+ X[0]*(X[4]*X[8]-X[5]*X[7])
		- X[1]*(X[3]*X[8]-X[5]*X[6])
		+ X[2]*(X[3]*X[7]-X[4]*X[6]);

	// Compute X^T.X (XX)

	for(i=0; i<3; i++)
	{
		for(j=0; j<3; j++)
		{
			for(k=0; k<3; k++)
			{
				XX[3*i+j] += X[3*k+i]*X[3*k+j];
			}
		}
	}

	// Compute the diagonalisation of XX using
	// the Jacobi algorithm...

	jacobi(XX,P,e);

	e[0]=max(e[0],0);
	e[1]=max(e[1],0);
	e[2]=max(e[2],0);


	sumE = sqrt(e[0]) + sqrt(e[1]) + sqrt(e[2]);

	if(detX<1e-8)
	{
		sumE -= (double)2*sqrt(min(e[0],min(e[1],e[2])));
		flag=1;
	}
	else
	{
		flag=0;
	}

	rmsd = sumA + sumB - (double)2*sumE;
	rmsd = sqrt(max(rmsd,0)/(double)n);

	// Compute the transform and return the rmsd.
	// This will fail if XX has rank<3. I need to
	// fix this...

	for(i=0; i<3; i++)
	{
		for(j=0; j<3; j++)
		{
			factor = sqrt(e[j]);
			if(flag && j==0)
			{
				factor = -factor;
			}

			T[3*i+j]=0.0;
			for(k=0; k<3; k++)
			{
				T[3*i+j] += X[3*i+k]*P[3*k+j]/factor;
			}
		}
	}

	for(i=0; i<3; i++)
	{
		for(j=0; j<3; j++)
		{
			M[3*i+j]=0.0;
			for(k=0; k<3; k++)
			{
				M[3*i+j] += P[3*i+k]*T[3*j+k];
			}
		}
	}

	return rmsd;
}

// ==================================================================
// Methods of type Node
// ==================================================================

static Node *Node_create(const double *x, const double *y)
{
	Node *n;
	if(!x || !y) return NULL;

	n=(Node*)malloc(sizeof(Node));
	memcpy(n->data[0],x,sizeof(double)*3);
	memcpy(n->data[1],y,sizeof(double)*3);
	return n;
}

static void Node_free(Node *n)
{
	if(n) free(n);
}

// ==================================================================

