// Delaunay.cpp: implementation of the Delaunay class.
//   based on Clarkson's hull
//////////////////////////////////////////////////////////////////////

#include "hull.h"
#include "Delaunay.h"
#include <stdio.h>
#include <math.h>

extern "C" {
int comp( const void *arg1, const void *arg2 ) {
	if (*(unsigned*)arg1 == *(unsigned*)arg2) return 0;
	if (*(unsigned*)arg1 > *(unsigned*)arg2) return 1;
	return -1;
}
}

void Delaunay::Init ( bool ext ) {

	// not enough points
	if (npoints < d+1) return;

	double minX=std::numeric_limits<double>::max();
	double maxX=std::numeric_limits<double>::min();
	for (unsigned i=0; i<X.size(); i++) 
		for (unsigned k=0; k<d; k++) {
			if ( X[i][k] < minX ) minX = X[i][k];
			if ( X[i][k] > maxX ) maxX = X[i][k];
		}

	std::vector<Eigen::Vector3d> scaledX(X.size());
	for (unsigned i=0; i<X.size(); i++)
		for (unsigned k=0; k<d; k++) scaledX[i][k] = (X[i][k]-minX)/(maxX-minX);

	ClarksonHull CH ( scaledX, npoints, Smplx );
	Smplx.Init(CH.GetNSimpices(),d,d);
	CH.Output();

	if (!ext) {

		Faces.Init(Smplx.GetNoObjects()*(d+1),d-1,d-1);
		Simplex *S;
		unsigned it;
		double minbox[d-1], maxbox[d-1];
		for (S=GetFirstSimplex(it); S; S=GetNextSimplex(it)) {
			unsigned pts[MAXDIM+1];
			for (unsigned i=0; i<d+1; i++) {
				S->GetFace(i,pts);
				Faces.Revise ( pts, minbox, maxbox );
			}
		}
	}

	Smplx.ScaleDown ( minX, maxX );
	Faces.ScaleDown ( minX, maxX );

}

Delaunay::~Delaunay() {

	if (NSimp) {
		delete[] PointSimplices;
		delete[] NSimp;
		delete[] SimplexFaces;
		delete[] face_simplex1;
		delete[] face_simplex2;
	}

}


SimplexList::SimplexList() : 
		buffer(0), smplx(0), smboxmin(0), smboxmax(0) {}


void SimplexList::Init ( unsigned nsmp, unsigned dim, unsigned dimbox ) {

	n = nsmp;
	d = dim;
	dbox = dimbox;
	curr = 0;
	buffer = new unsigned [n*(d+1)];
	smboxmin = new double[n*dbox];
	smboxmax = new double[n*dbox];
	smplx = new Simplex[n];
	for (unsigned i=0; i<n; i++) {
		smplx[i].np = d+1;
		smplx[i].points = buffer + i*(d+1);
		smplx[i].boxmin = &smboxmin[i*dbox];
		smplx[i].boxmax = &smboxmax[i*dbox];
	}

}


SimplexList::~SimplexList() {

	if (buffer) delete[] buffer;
	if (smplx) delete[] smplx;
	if (smboxmin) delete[] smboxmin;
	if (smboxmax) delete[] smboxmax;

}


Simplex *SimplexList::Add ( unsigned *pts, double *minbox, double *maxbox ) {

	if ((unsigned)curr == n) return 0;
	qsort( (void*)pts, (size_t)d+1, sizeof(unsigned), comp );
	for (unsigned i=0; i<d+1; i++)
		smplx[curr].points[i] = pts[i];
	memcpy( smplx[curr].boxmin, minbox, dbox*sizeof(double) );
	memcpy( smplx[curr].boxmax, maxbox, dbox*sizeof(double) );

	return &smplx[curr++];

}


Simplex *SimplexList::Revise ( unsigned *pts, double *minbox, double *maxbox, Simplex **SS, ulong *m, char *checked ) {

	unsigned it;
	int i;
	Simplex *S;
	for (S=GetFirst(it),i=0; S; S=GetNext(it),i++) {
		if ( checked[i] ) continue;
		if (S->TheSame(pts)) {
			checked[i] = 1;
			*SS = S;
			*m = it-1;
			return 0;
		}
	}

	Simplex *ret = Add ( pts, minbox, maxbox );
	checked[i]=0;
	*SS = ret;
	return ret;

}

Simplex *SimplexList::Revise ( unsigned *pts, double *minbox, double *maxbox ) {

	unsigned it;
	int i;
	Simplex *S;
	for (S=GetFirst(it),i=0; S; S=GetNext(it),i++) {
		if (S->TheSame(pts)) {
			return 0;
		}
	}

	Simplex *ret = Add ( pts, minbox, maxbox );
	return ret;

}


Simplex *SimplexList::Revise ( Simplex &S ) {

	return Revise ( S.GetPoints(), S.GetBoxMin(), S.GetBoxMax() );

}


bool SimplexList::IsIn ( Simplex *S ) {

	unsigned it;
	for (Simplex *s=GetFirst(it); s; s=GetNext(it))
		if (s->TheSame(S)) return true;

	return false;

}

void SimplexList::ScaleDown ( double minX, double maxX ) {

	unsigned it;
	for (Simplex *s=GetFirst(it); s; s=GetNext(it))
		for (unsigned k=0; k<d; k++) {
			s->boxmin[k] = minX + s->boxmin[k]*(maxX-minX);
			s->boxmax[k] = minX + s->boxmax[k]*(maxX-minX);
		}

}
 

unsigned Simplex::GetFace ( unsigned i, unsigned *pts ) {

	for (unsigned j=0, k=0; j<np; j++)
		if (i!=j) pts[k++] = points[j];

	return points[i];
}


void Simplex::GetFace1 ( unsigned i, unsigned *pts ) {

	for (unsigned j=0, k=0; j<np; j++)
		if (points[j]!=i) pts[k++] = points[j];

}


unsigned Simplex::GetPointToFace ( Simplex &face ) {

	int i,j, ret[MAXPOINTS];

	for (i=0; i<(int)np; i++) ret[i] = points[i];

	for (i=0; i<(int)np; i++)
		for (j=0; j<(int)np-1; j++)
			if ( face.GetPoint(j) == (unsigned)ret[i] ) 
				ret[i] = -1;

	for (i=0; i<(int)np; i++) if ( ret[i] != -1 ) return ret[i];

	return -1;
}


bool Simplex::TheSame ( unsigned *pts ) {

	for (unsigned i=0; i<np; i++)
		if (points[i] != pts[i]) return false;

	return true;

}


bool Simplex::TheSame ( const Simplex *S ) {

	for (unsigned i=0; i<np; i++)
		if (points[i] != S->points[i]) return false;

	return true;

}


bool Simplex::ContainsThePoint ( unsigned pnt ) {

	for (unsigned i=0; i<np; i++)
		if (points[i] == pnt) return true;

	return false;

}


// from here: https://www.sanfoundry.com/cpp-program-compute-determinant-matrix/
double det(int n, double mat[H_MAXDIM][H_MAXDIM])
{
    double d = 0;
    double submat[H_MAXDIM][H_MAXDIM];
    if (n == 2)
        return ((mat[0][0] * mat[1][1]) - (mat[1][0] * mat[0][1]));
    else
    {
        for (int c = 0; c < n; c++)
        {
            int subi = 0; //submatrix's i value
            for (int i = 1; i < n; i++)
            {
                int subj = 0;
                for (int j = 0; j < n; j++)
                {
                    if (j == c)
                        continue;
                    submat[subi][subj] = mat[i][j];
                    subj++;
                }
                subi++;
 
            }
            d = d + (pow(-1, c) * mat[0][c] * det(n - 1, submat));
        }
    }
    return d;
}


// from here: https://en.wikipedia.org/wiki/Simplex#Volume
double Delaunay::CalculateSimplexVolume ( Simplex &S ) {

	double mat[H_MAXDIM][H_MAXDIM];

	const Eigen::Vector3d &x0 = X[S.GetPoint(0)];
	for (unsigned i=1; i<d+1; i++) {
		const Eigen::Vector3d &xi = X[S.GetPoint(i)];
		for (unsigned j=0; j<d; j++) mat[i-1][j] = xi[j] - x0[j];
	}
	
	long factorial=1;
        for (unsigned i=1; i<=d; ++i) factorial *= i;

	return fabs(det(d, mat))/factorial;
		
}


void Delaunay::GetSimplexPoints ( int simplex, Eigen::Vector3d *x ) {

	for (unsigned i=0; i<d+1; i++)
		x[i] = X[Smplx[simplex].GetPoint(i)];

}


void Delaunay::GetFaceSimplices ( int face, Eigen::Vector3d *S1, Eigen::Vector3d *S2 ) {

	if ( face_simplex1[face] == -1 ) 
		S1[0] = Eigen::Vector3d(4294967295,4294967295,4294967295);
	else 
		for (unsigned i=0; i<d+1; i++)
			S1[i] = X[Smplx[face_simplex1[face]].GetPoint(i)];

	if ( face_simplex2[face] == -1 ) 
		S2[0] = Eigen::Vector3d(4294967295,4294967295,4294967295);
	else 
		for (unsigned i=0; i<d+1; i++)
			S2[i] = X[Smplx[face_simplex2[face]].GetPoint(i)];

}


void Delaunay::GetFaceSimplices ( int face, int &S1, int &S2 ) {

	S1 = face_simplex1[face];
	S2 = face_simplex2[face];

}


int Delaunay::GetFaceSimplices ( int face, unsigned *pts1, unsigned *pts2 ) {

	int ret = FACESMPLX_BOTH;

	if ( face_simplex1[face] == -1 ) 
		ret = FACESMPLX_NOLEFT;
	else 
		for (unsigned i=0; i<d+1; i++)
			pts1[i] = Smplx[face_simplex1[face]].GetPoint(i);

	if ( face_simplex2[face] == -1 ) 
		ret = FACESMPLX_NORIGHT;
	else 
		for (unsigned i=0; i<d+1; i++)
			pts2[i] = Smplx[face_simplex2[face]].GetPoint(i);

	return ret;

}


bool Delaunay::IsBoundary ( int face ) {

	return face_simplex1[face]==-1 || face_simplex2[face]==-1;

}


inline bool inlist ( int np, int *list, int item ) {
	for (int i=0; i<np; i++) if ( list[i] == item ) return true;
	return false;
}


void Delaunay::GetStar ( int pnt, int *star, int &ns, int *ring, int &nr ) {

	static unsigned pts[MAXDIM+1];

	nr = 0;
	ns = 0;

	for (int i=0; i<NSimp[pnt]; i++) {

		int simpnum = PointSimplices[pnt][i];
		Smplx[simpnum].GetFace1(pnt,pts);

		for (int j=0; j<(int)d+1; j++) {
			int facenum = SimplexFaces[simpnum][j];
			if ( !Faces[facenum].TheSame(pts) ) {
				if ( !inlist(ns,star,facenum) ) star[ns++] = facenum;
			}
			else ring[nr++] = facenum;
		}

	}

}


void Delaunay::GetRing ( int pnt, int *ring, int &nr ) {

	static unsigned pts[MAXDIM+1];

	nr = 0;

	for (int i=0; i<NSimp[pnt]; i++) {

		int simpnum = PointSimplices[pnt][i];
		Smplx[simpnum].GetFace1(pnt,pts);

		for (int j=0; j<(int)d+1; j++) {
			int facenum = SimplexFaces[simpnum][j];
			if ( Faces[facenum].TheSame(pts) ) 
				ring[nr++] = facenum;
		}

	}

}


void Delaunay::GetSimplexFaces ( int simplex, int *faces ) {

	for (int i=0; i<(int)d+1; i++) faces[i] = SimplexFaces[simplex][i];

}


double Delaunay::GetSimplexSize ( int simplex ) {

	// the volume of the sircumbscribing box

	double vol=1;

	for (int i=0; i<(int)d; i++)
		vol *= Smplx[simplex].GetBoxMax()[i] - Smplx[simplex].GetBoxMin()[i];

	return vol;

}


double Delaunay::GetFaceSize ( int face ) {

	// the volume of the sircumbscribing box

	double vol=1;

	for (int i=0; i<(int)d; i++)
		vol *= Faces[face].GetBoxMax()[i] - Faces[face].GetBoxMin()[i];

	return vol;

}


Simplex &Delaunay::GetBiggestFace ( int simplex ) {

	double max=-1e32, size;
	int imax;

	for (int i=0; i<(int)d+1; i++) {
		size = GetFaceSize(SimplexFaces[simplex][i]);
		if ( size>max ) { max = size; imax = SimplexFaces[simplex][i]; }
	}

	return GetFace(imax);

}


void Delaunay::GetFacesWithEdge ( int simplex, Simplex &edge, int &face1, int &face2 ) {

	int face, i,j,k,k1;

	face1 = face2 = -1;

	// for each face
	for (i=0; i<(int)d+1; i++) {

		face = SimplexFaces[simplex][i];

		// for each edge of the face
		for (j=0; j<(int)d; j++) {
			k1 = 0;
			unsigned pts[MAXDIM];
			for (k=0; k<(int)d; k++)
				if (k!=j) pts[k1++] = Faces[face].GetPoint(k);

			if ( edge.TheSame(pts) ) {
				if ( face1 == -1 ) face1 = face;
				else face2 = face;
			}
		}
				
	}

}


bool Delaunay::IsBoundaryPoint ( int pnt ) {

	int i,ns,nr;

	int *star = new int[d*GetMaxNSimp()];
	int *ring = new int[GetMaxNSimp()];

	GetStar ( pnt, star, ns, ring, nr );

	for (i=0; i<ns; i++)
		if ( IsBoundary(star[i]) ) {
			delete[] star;
			delete[] ring;
			return true;
		}

	delete[] star;
	delete[] ring;
	return false;

}


int Delaunay::FindEdge ( int face, SimplexList &edges ) {

	int i,j,k;
	Simplex &Face = GetFace(face);
	unsigned pts[MAXDIM+1];

	for (i=0; i<(int)d; i++) {
		k=0;
		for (j=0; j<(int)d; j++) 
			if ( j!=i ) pts[k++] = Face.GetPoint(j);
		for (j=0; j<(int)d; j++)
			if ( edges[j].TheSame(pts) ) return j;
	}

	return -1;

}

unsigned Delaunay::GetMaxNSimp() { 

	unsigned nmax=0;

	for (unsigned i=0; i<npoints; i++) {

		unsigned it; Simplex *S; 
		unsigned ns=0;
		for (S=Smplx.GetFirst(it); S; S=Smplx.GetNext(it)) 
			if ( S->ContainsThePoint(i) ) ns++;

		if ( nmax<ns ) nmax = ns;

	}

	return nmax;

}
