// Delaunay.h: interface for the Delaunay class.
//   based on Clarkson's hull
//////////////////////////////////////////////////////////////////////

#ifndef _DELAUNAY_H
#define _DELAUNAY_H

#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#include <math.h>
#include <vector>
#include "Eigen/Dense"
 
#ifndef ulong
typedef unsigned long ulong;
#endif


class SimplexList;


class Simplex {

		friend class SimplexList;

	protected:

		unsigned np;
		unsigned *points;
		double *boxmin, *boxmax;

	public:

		bool TheSame ( unsigned *pts );
		bool TheSame ( const Simplex *S );

		unsigned GetPoint ( unsigned i ) { return points[i]; }
		unsigned GetFace ( unsigned i, unsigned *pts );
		void GetFace1 ( unsigned i, unsigned *pts );
		double *GetBoxMin() { return boxmin; }
		double *GetBoxMax() { return boxmax; }
		unsigned *GetPoints() { return points; }
		unsigned GetNPoints() { return np; }
		unsigned GetPointToFace ( Simplex &face );
		bool ContainsThePoint ( unsigned pnt );

};


class SimplexList {

		unsigned n, d, dbox;
		unsigned curr;
		unsigned *buffer;
		Simplex *smplx;
		double *smboxmin, *smboxmax;	// may not be of the d dimensionality
						// if the simplices are in higher
						// dimensional space

	public:

		SimplexList();
		void Init ( unsigned nsmp, unsigned dim, unsigned dimbox );
		~SimplexList();
		Simplex *Add ( unsigned *pts, double *minbox, double *maxbox );
		Simplex *Revise ( unsigned *pts, double *minbox, double *maxbox, Simplex **S, ulong *m, char *checked );
		Simplex *Revise ( unsigned *pts, double *minbox, double *maxbox );
		Simplex *Revise ( Simplex &S );
		inline Simplex *GetFirst ( unsigned &it ) { if (!curr) return 0; it=0; return &smplx[it++]; }
		inline Simplex *GetNext ( unsigned &it ) { if (it==curr) return 0; return &smplx[it++]; }
		Simplex &operator[] (int i) { return smplx[i]; }
		Simplex &operator[] (unsigned i) { return smplx[i]; }
		unsigned GetNoObjects() { return curr; }
		bool IsIn ( Simplex *S );
		void ScaleDown ( double minX, double maxX );

};


#define MAXDIM		150
#define MAXPOINTS	MAXDIM+5

#define FACESMPLX_BOTH		0
#define FACESMPLX_NOLEFT	1
#define FACESMPLX_NORIGHT	2


class Delaunay  
{
	unsigned npoints; 
	unsigned d;			// dimensionality
	Eigen::Vector3d Min, Max, Centre;
	const std::vector<Eigen::Vector3d> &X;	// points

	int 	**PointSimplices;	// list of indices to "simplices" 
					//    (d-simplices) for each point

	char *NSimp;			// number of "simplices" for each point

	SimplexList	Faces;		// list of all "faces" ( (d-1)-simplices )

	int	*face_simplex1,		// pointers to two simplices to which the face  
		*face_simplex2;		//     belongs (for each face)
						
	int 	**SimplexFaces;		// indices of the faces of each simplex

	unsigned MaxNSimp;		// the biggest number of simplices per point
	unsigned NBoundary;		// number of boundary faces

public:
	Delaunay( const std::vector<Eigen::Vector3d> &x, bool ext=false, bool dbg=false ) :
		npoints(x.size()), d(x[0].size()),
		Min(1e32,1e32,1e32), Max(-1e32,-1e32,-1e32),
		Centre(0,0,0), X(x), NSimp(0) { Init(ext); }

	~Delaunay();

	Simplex *GetFirstSimplex(unsigned &it) { return Smplx.GetFirst(it); }
	Simplex *GetNextSimplex(unsigned &it) { return Smplx.GetNext(it); }
	Eigen::Vector3d GetPoint ( unsigned i ) { return X[i]; }
	Eigen::Vector3d GetBoxMax() { return Max; }
	Eigen::Vector3d GetBoxMin() { return Min; }
	int GetNSimplices() { return Smplx.GetNoObjects(); }
	int GetNFaces() { return Faces.GetNoObjects(); }
	unsigned GetDimensionality() { return d; }
	unsigned GetNPoints() { return npoints; }
	void AddNextPoint() { curr_point--; }
	Eigen::Vector3d GetCenter() { return Centre; }
	Simplex &operator[] (int i) { return Smplx[i]; }
	Simplex &operator[] (unsigned i) { return Smplx[i]; }
	void GetSimplexPoints ( int simplex, Eigen::Vector3d *x ); 
		// assumed that the space for d+1 d-dimensional points is allocated
	void GetFaceSimplices ( int face, Eigen::Vector3d *S1, Eigen::Vector3d *S2 );
	void GetFaceSimplices ( int face, int &S1, int &S2 );
	bool IsBoundary ( int face );
		// return value indicates which simplex exists
	int GetFaceSimplices ( int face, unsigned *pts1, unsigned *pts2 );
	unsigned GetMaxNSimp(); 
		// calculates the maximal number of simplices belonging to a point
	unsigned GetNBoundary() { return NBoundary; }
	void GetStar ( int pnt, int *star, int &ns, int *ring, int &nr );
	void GetStar ( int pnt, int *star, int &ns ) { 
		for (int i=0; i<NSimp[pnt]; i++) star[i]=PointSimplices[pnt][i]; 
		ns=NSimp[pnt]; }
	void GetRing ( int pnt, int *ring, int &nr );
	Simplex &GetFace ( int face ) { return Faces[face]; }
	Simplex &GetSimplex ( int i ) { return Smplx[i]; }
	void GetSimplexFaces ( int simplex, int *faces );
	double GetSimplexSize ( int simplex );
	double GetFaceSize ( int face );
	Simplex &GetBiggestFace ( int simplex );
	int GetSecondSimpOfFace ( int face, int simplex ) {
		if ( face_simplex1[face] == -1 ) return face_simplex1[face];
		if ( face_simplex2[face] == -1 ) return face_simplex2[face];
		return face_simplex1[face]==simplex ? face_simplex2[face] : face_simplex1[face]; }
	void GetFacesWithEdge ( int simplex, Simplex &edge, int &face1, int &face2 );
		// there are always two of them
	bool IsBoundaryPoint ( int pnt );
	int FindEdge ( int face, SimplexList &edges );
	double CalculateSimplexVolume ( Simplex &S );

	void GetSimplicesOfPoint ( int pnt, std::vector<int> &smpls, int &ns ) { 
		unsigned it; Simplex *S; 
		ns=0;
		for (S=Smplx.GetFirst(it); S; S=Smplx.GetNext(it)) 
			if ( S->ContainsThePoint(pnt) ) smpls[ns++]=it-1;
	}

private:

	void Init ( bool ext = false );
	SimplexList Smplx;
	int curr_point;

};


#endif

