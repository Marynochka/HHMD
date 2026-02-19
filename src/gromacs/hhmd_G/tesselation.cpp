#include "tesselation.h"


#define NMaxSimplPerPoint 100


RealSimplex::RealSimplex ( 
	int point,
	Simplex &oS, 
	std::vector<Eigen::Vector3d> &oX, 
	Eigen::Vector3d sh ) : S(&oS), rX(&oX), shift(sh) {

	point_nu = point;
	int n=0;
	for (int i=0; i<4; i++)
		if (S->GetPoint(i)==point) continue;
		else other_points[n++] = S->GetPoint(i);

	std::vector<Eigen::Vector3d> &X = *rX;

	for (int i=0; i<3; i++) { Min[i]=1e32; Max[i]=-1e32; }
	Eigen::Vector3d x ( X[point_nu]+shift );
	for (int i=0; i<3; i++) { 
		if (x[i]<Min[i]) Min[i]=x[i];
		if (x[i]>Max[i]) Max[i]=x[i];
	}
	for (int k=0; k<3; k++) {
		x = X[other_points[k]]+shift;
		for (int i=0; i<3; i++) { 
			if (x[i]<Min[i]) Min[i]=x[i];
			if (x[i]>Max[i]) Max[i]=x[i];
		}
	}

}


void RealSimplex::PrintSimplex() {

	FILE* f = fopen("simplex.dat", "w");

	PrintSimplex(f);

	fclose(f);

}


void RealSimplex::PrintSimplex ( FILE* f ) {

	std::vector<Eigen::Vector3d> &X = *rX;
	Eigen::Vector3d r_nu, r_mu, r_gamma, r_sigma;
	r_nu 	= X[point_nu]+shift;
	r_mu 	= X[other_points[0]]+shift;
	r_gamma	= X[other_points[1]]+shift;
	r_sigma	= X[other_points[2]]+shift;

	fprintf ( f, "%.3f %.3f %.3f\n", r_nu[0],r_nu[1],r_nu[2] );
	fprintf ( f, "%.3f %.3f %.3f\n", r_mu[0],r_mu[1],r_mu[2] );
	fprintf ( f, "%.3f %.3f %.3f\n", r_gamma[0],r_gamma[1],r_gamma[2] );
	fprintf ( f, "%.3f %.3f %.3f\n", r_sigma[0],r_sigma[1],r_sigma[2] );
/*	fprintf ( f, "%.3f %.3f %.3f\n", r_nu[0],r_nu[1],r_nu[2] );
	fprintf ( f, "%.3f %.3f %.3f\n", r_sigma[0],r_sigma[1],r_sigma[2] );
	fprintf ( f, "%.3f %.3f %.3f\n", r_mu[0],r_mu[1],r_mu[2] );
	fprintf ( f, "%.3f %.3f %.3f\n", r_gamma[0],r_gamma[1],r_gamma[2] );
	fprintf ( f, "%.3f %.3f %.3f\n", r_sigma[0],r_sigma[1],r_sigma[2] );
	fprintf ( f, "\n\n\n" );
*/
}


void Star::Build ( int p, Delaunay &D, std::vector<Eigen::Vector3d> &X, Tesselation &tssl ) {

	point = p;

	nRS = 0;
	RS.resize(2*NMaxSimplPerPoint);

	int ns;
	std::vector<int> smpls(3*D.GetMaxNSimp());
	Eigen::Vector3d shift;

	D.GetSimplicesOfPoint ( point, smpls, ns );

	shift[0]=0; shift[1]=0; shift[2]=0;
	for (int i=0; i<ns; i++) 
		RS[nRS++] = RealSimplex ( point, D[smpls[i]], X, shift );
/*
	// sides    !!! tested only for nx=ny=nz
	if ( tssl.IsMinX(point) ) {
		D.GetSimplicesOfPoint ( tssl.PeriodicMinX(point), smpls, ns );
		shift[0]=-(tssl.Max[0]-tssl.Min[0]); shift[1]=0; shift[2]=0;
		for (int i=0; i<ns; i++) 
			RS[nRS++] = RealSimplex ( tssl.PeriodicMinX(point), D[smpls[i]], X, shift );
	}

	if ( tssl.IsMaxX(point) ) {
		D.GetSimplicesOfPoint ( tssl.PeriodicMaxX(point), smpls, ns );
		shift[0]=(tssl.Max[0]-tssl.Min[0]); shift[1]=0; shift[2]=0;
		for (int i=0; i<ns; i++) 
			RS[nRS++] = RealSimplex ( tssl.PeriodicMaxX(point), D[smpls[i]], X, shift );
	}

	if ( tssl.IsMinY(point) ) {
		D.GetSimplicesOfPoint ( tssl.PeriodicMinY(point), smpls, ns );
		shift[0]=0; shift[1]=-(tssl.Max[1]-tssl.Min[1]); shift[2]=0;
		for (int i=0; i<ns; i++) 
			RS[nRS++] = RealSimplex ( tssl.PeriodicMinY(point), D[smpls[i]], X, shift );
	}

	if ( tssl.IsMaxY(point) ) {
		D.GetSimplicesOfPoint ( tssl.PeriodicMaxY(point), smpls, ns );
		shift[0]=0; shift[1]=(tssl.Max[1]-tssl.Min[1]); shift[2]=0;
		for (int i=0; i<ns; i++) 
			RS[nRS++] = RealSimplex ( tssl.PeriodicMaxY(point), D[smpls[i]], X, shift );
	}

	if ( tssl.IsMinZ(point) ) {
		D.GetSimplicesOfPoint ( tssl.PeriodicMinZ(point), smpls, ns );
		shift[0]=0; shift[1]=0; shift[2]=-(tssl.Max[2]-tssl.Min[2]);
		for (int i=0; i<ns; i++) 
			RS[nRS++] = RealSimplex ( tssl.PeriodicMinZ(point), D[smpls[i]], X, shift );
	}

	if ( tssl.IsMaxZ(point) ) {
		D.GetSimplicesOfPoint ( tssl.PeriodicMaxZ(point), smpls, ns );
		shift[0]=0; shift[1]=0; shift[2]=(tssl.Max[2]-tssl.Min[2]);
		for (int i=0; i<ns; i++) 
			RS[nRS++] = RealSimplex ( tssl.PeriodicMaxZ(point), D[smpls[i]], X, shift );
	}

	// edges    !!! tested only for nx=ny=nz
	if ( tssl.IsMinX(point) && tssl.IsMinY(point) )  {
		D.GetSimplicesOfPoint ( tssl.PeriodicMinXMinY(point), smpls, ns );
		shift[0]=-(tssl.Max[0]-tssl.Min[0]); shift[1]=-(tssl.Max[1]-tssl.Min[1]); shift[2]=0;
		for (int i=0; i<ns; i++) 
			RS[nRS++] = RealSimplex ( tssl.PeriodicMinXMinY(point), D[smpls[i]], X, shift );
	}

	if ( tssl.IsMinX(point) && tssl.IsMaxY(point) )  {
		D.GetSimplicesOfPoint ( tssl.PeriodicMinXMaxY(point), smpls, ns );
		shift[0]=-(tssl.Max[0]-tssl.Min[0]); shift[1]=(tssl.Max[1]-tssl.Min[1]); shift[2]=0;
		for (int i=0; i<ns; i++) 
			RS[nRS++] = RealSimplex ( tssl.PeriodicMinXMaxY(point), D[smpls[i]], X, shift );
	}

	if ( tssl.IsMaxX(point) && tssl.IsMinY(point) )  {
		D.GetSimplicesOfPoint ( tssl.PeriodicMaxXMinY(point), smpls, ns );
		shift[0]=(tssl.Max[0]-tssl.Min[0]); shift[1]=-(tssl.Max[1]-tssl.Min[1]); shift[2]=0;
		for (int i=0; i<ns; i++) 
			RS[nRS++] = RealSimplex ( tssl.PeriodicMaxXMinY(point), D[smpls[i]], X, shift );
	}

	if ( tssl.IsMaxX(point) && tssl.IsMaxY(point) )  {
		D.GetSimplicesOfPoint ( tssl.PeriodicMaxXMaxY(point), smpls, ns );
		shift[0]=(tssl.Max[0]-tssl.Min[0]); shift[1]=(tssl.Max[1]-tssl.Min[1]); shift[2]=0;
		for (int i=0; i<ns; i++) 
			RS[nRS++] = RealSimplex ( tssl.PeriodicMaxXMaxY(point), D[smpls[i]], X, shift );
	}

	if ( tssl.IsMinX(point) && tssl.IsMinZ(point) )  {
		D.GetSimplicesOfPoint ( tssl.PeriodicMinXMinZ(point), smpls, ns );
		shift[0]=-(tssl.Max[0]-tssl.Min[0]); shift[1]=0; shift[2]=-(tssl.Max[2]-tssl.Min[2]);
		for (int i=0; i<ns; i++) 
			RS[nRS++] = RealSimplex ( tssl.PeriodicMinXMinZ(point), D[smpls[i]], X, shift );
	}

	if ( tssl.IsMinX(point) && tssl.IsMaxZ(point) )  {
		D.GetSimplicesOfPoint ( tssl.PeriodicMinXMaxZ(point), smpls, ns );
		shift[0]=-(tssl.Max[0]-tssl.Min[0]); shift[1]=0; shift[2]=(tssl.Max[2]-tssl.Min[2]);
		for (int i=0; i<ns; i++) 
			RS[nRS++] = RealSimplex ( tssl.PeriodicMinXMaxZ(point), D[smpls[i]], X, shift );
	}

	if ( tssl.IsMaxX(point) && tssl.IsMinZ(point) )  {
		D.GetSimplicesOfPoint ( tssl.PeriodicMaxXMinZ(point), smpls, ns );
		shift[0]=(tssl.Max[0]-tssl.Min[0]); shift[1]=0; shift[2]=-(tssl.Max[2]-tssl.Min[2]);
		for (int i=0; i<ns; i++) 
			RS[nRS++] = RealSimplex ( tssl.PeriodicMaxXMinZ(point), D[smpls[i]], X, shift );
	}

	if ( tssl.IsMaxX(point) && tssl.IsMaxZ(point) )  {
		D.GetSimplicesOfPoint ( tssl.PeriodicMaxXMaxZ(point), smpls, ns );
		shift[0]=(tssl.Max[0]-tssl.Min[0]); shift[1]=0; shift[2]=(tssl.Max[2]-tssl.Min[2]);
		for (int i=0; i<ns; i++) 
			RS[nRS++] = RealSimplex ( tssl.PeriodicMaxXMaxZ(point), D[smpls[i]], X, shift );
	}

	if ( tssl.IsMinY(point) && tssl.IsMinZ(point) )  {
		D.GetSimplicesOfPoint ( tssl.PeriodicMinYMinZ(point), smpls, ns );
		shift[0]=0; shift[1]=-(tssl.Max[1]-tssl.Min[1]); shift[2]=-(tssl.Max[2]-tssl.Min[2]);
		for (int i=0; i<ns; i++) 
			RS[nRS++] = RealSimplex ( tssl.PeriodicMinYMinZ(point), D[smpls[i]], X, shift );
	}

	if ( tssl.IsMinY(point) && tssl.IsMaxZ(point) )  {
		D.GetSimplicesOfPoint ( tssl.PeriodicMinYMaxZ(point), smpls, ns );
		shift[0]=0; shift[1]=-(tssl.Max[1]-tssl.Min[1]); shift[2]=(tssl.Max[2]-tssl.Min[2]);
		for (int i=0; i<ns; i++) 
			RS[nRS++] = RealSimplex ( tssl.PeriodicMinYMaxZ(point), D[smpls[i]], X, shift );
	}

	if ( tssl.IsMaxY(point) && tssl.IsMinZ(point) )  {
		D.GetSimplicesOfPoint ( tssl.PeriodicMaxYMinZ(point), smpls, ns );
		shift[0]=0; shift[1]=(tssl.Max[1]-tssl.Min[1]); shift[2]=-(tssl.Max[2]-tssl.Min[2]);
		for (int i=0; i<ns; i++) 
			RS[nRS++] = RealSimplex ( tssl.PeriodicMaxYMinZ(point), D[smpls[i]], X, shift );
	}

	if ( tssl.IsMaxY(point) && tssl.IsMaxZ(point) )  {
		D.GetSimplicesOfPoint ( tssl.PeriodicMaxYMaxZ(point), smpls, ns );
		shift[0]=0; shift[1]=(tssl.Max[1]-tssl.Min[1]); shift[2]=(tssl.Max[2]-tssl.Min[2]);
		for (int i=0; i<ns; i++) 
			RS[nRS++] = RealSimplex ( tssl.PeriodicMaxYMaxZ(point), D[smpls[i]], X, shift );
	}


	// vertices    !!! tested only for nx=ny=nz
	if ( tssl.IsMinX(point) && tssl.IsMinY(point) && tssl.IsMinZ(point) )  {
		D.GetSimplicesOfPoint ( tssl.PeriodicMinXMinYMinZ(point), smpls, ns );
		shift[0]=-(tssl.Max[0]-tssl.Min[0]); shift[1]=-(tssl.Max[1]-tssl.Min[1]); shift[2]=-(tssl.Max[2]-tssl.Min[2]);
		for (int i=0; i<ns; i++) 
			RS[nRS++] = RealSimplex ( tssl.PeriodicMinXMinYMinZ(point), D[smpls[i]], X, shift );
	}

	if ( tssl.IsMinX(point) && tssl.IsMinY(point) && tssl.IsMaxZ(point) )  {
		D.GetSimplicesOfPoint ( tssl.PeriodicMinXMinYMaxZ(point), smpls, ns );
		shift[0]=-(tssl.Max[0]-tssl.Min[0]); shift[1]=-(tssl.Max[1]-tssl.Min[1]); shift[2]=(tssl.Max[2]-tssl.Min[2]);
		for (int i=0; i<ns; i++) 
			RS[nRS++] = RealSimplex ( tssl.PeriodicMinXMinYMaxZ(point), D[smpls[i]], X, shift );
	}

	if ( tssl.IsMinX(point) && tssl.IsMaxY(point) && tssl.IsMinZ(point) )  {
		D.GetSimplicesOfPoint ( tssl.PeriodicMinXMaxYMinZ(point), smpls, ns );
		shift[0]=-(tssl.Max[0]-tssl.Min[0]); shift[1]=(tssl.Max[1]-tssl.Min[1]); shift[2]=-(tssl.Max[2]-tssl.Min[2]);
		for (int i=0; i<ns; i++) 
			RS[nRS++] = RealSimplex ( tssl.PeriodicMinXMaxYMinZ(point), D[smpls[i]], X, shift );
	}

	if ( tssl.IsMinX(point) && tssl.IsMaxY(point) && tssl.IsMaxZ(point) )  {
		D.GetSimplicesOfPoint ( tssl.PeriodicMinXMaxYMaxZ(point), smpls, ns );
		shift[0]=-(tssl.Max[0]-tssl.Min[0]); shift[1]=(tssl.Max[1]-tssl.Min[1]); shift[2]=(tssl.Max[2]-tssl.Min[2]);
		for (int i=0; i<ns; i++) 
			RS[nRS++] = RealSimplex ( tssl.PeriodicMinXMaxYMaxZ(point), D[smpls[i]], X, shift );
	}


	if ( tssl.IsMaxX(point) && tssl.IsMinY(point) && tssl.IsMinZ(point) )  {
		D.GetSimplicesOfPoint ( tssl.PeriodicMaxXMinYMinZ(point), smpls, ns );
		shift[0]=(tssl.Max[0]-tssl.Min[0]); shift[1]=-(tssl.Max[1]-tssl.Min[1]); shift[2]=-(tssl.Max[2]-tssl.Min[2]);
		for (int i=0; i<ns; i++) 
			RS[nRS++] = RealSimplex ( tssl.PeriodicMaxXMinYMinZ(point), D[smpls[i]], X, shift );
	}

	if ( tssl.IsMaxX(point) && tssl.IsMinY(point) && tssl.IsMaxZ(point) )  {
		D.GetSimplicesOfPoint ( tssl.PeriodicMaxXMinYMaxZ(point), smpls, ns );
		shift[0]=(tssl.Max[0]-tssl.Min[0]); shift[1]=-(tssl.Max[1]-tssl.Min[1]); shift[2]=(tssl.Max[2]-tssl.Min[2]);
		for (int i=0; i<ns; i++) 
			RS[nRS++] = RealSimplex ( tssl.PeriodicMaxXMinYMaxZ(point), D[smpls[i]], X, shift );
	}

	if ( tssl.IsMaxX(point) && tssl.IsMaxY(point) && tssl.IsMinZ(point) )  {
		D.GetSimplicesOfPoint ( tssl.PeriodicMaxXMaxYMinZ(point), smpls, ns );
		shift[0]=(tssl.Max[0]-tssl.Min[0]); shift[1]=(tssl.Max[1]-tssl.Min[1]); shift[2]=-(tssl.Max[2]-tssl.Min[2]);
		for (int i=0; i<ns; i++) 
			RS[nRS++] = RealSimplex ( tssl.PeriodicMaxXMaxYMinZ(point), D[smpls[i]], X, shift );
	}

	if ( tssl.IsMaxX(point) && tssl.IsMaxY(point) && tssl.IsMaxZ(point) )  {
		D.GetSimplicesOfPoint ( tssl.PeriodicMaxXMaxYMaxZ(point), smpls, ns );
		shift[0]=(tssl.Max[0]-tssl.Min[0]); shift[1]=(tssl.Max[1]-tssl.Min[1]); shift[2]=(tssl.Max[2]-tssl.Min[2]);
		for (int i=0; i<ns; i++) 
			RS[nRS++] = RealSimplex ( tssl.PeriodicMaxXMaxYMaxZ(point), D[smpls[i]], X, shift );
	}
*/
	for (int i=0; i<3; i++) { Min[i]=1e32; Max[i]=-1e32; }
	for (int k=0; k<nRS; k++) RS[k].UpdateMinMax(Min,Max);

} 


Tesselation::Tesselation ( const std::vector<Eigen::Vector3d> &oX, double mi[3], double ma[3], int N[3] ) :
	X(oX), D(X), nx(N[0]), ny(N[1]), nz(N[2]) {

	for (int i=0; i<3; i++) { Min[i] = mi[i]; Max[i] = ma[i]; }

	St.resize(X.size());

	for (int k=0; k<X.size(); k++) //{ 
		St[k].Build ( k, D, X, *this ); 
// St[k].PrintAllSimplexesVolumes();}

//	double _r[3] = { 7.87751533322496300000, 0.22858474917249835000, 6.62005182668404310000 };
//	r = Eigen::Vector3d(_r);
//if ( St[k].in_the_box(r) ) St[k].PrintAllSimplexes(); }

}


bool Tesselation::in_this_cell ( int k, double _r[3] ) { 

	r = Eigen::Vector3d(_r);
	if ( r[0] <= 0 ) r[0] = 0.0000000000000000001;
	if ( r[1] <= 0 ) r[1] = 0.0000000000000000001;
	if ( r[2] <= 0 ) r[2] = 0.0000000000000000001;

	if ( !St[k].in_the_box(r) ) return false;

	for (int i=0; i<St[k].GetNSimplices(); i++) {
		if ( !St[k].in_simplex_box(i,r) ) continue;
		if ( !St[k].theta(i,r) ) continue;
		curr_star = k; found_simplex = i;
		return true;
	}

	return false; 

}

