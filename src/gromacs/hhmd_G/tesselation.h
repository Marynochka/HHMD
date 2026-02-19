#ifndef _TESSELATION_H
#define _TESSELATION_H

#include "Delaunay.h"


class RealSimplex {

		Simplex *S;
		std::vector<Eigen::Vector3d> *rX;
		Eigen::Vector3d shift;
		Eigen::Vector3d Min, Max;

		int point_nu;
		int other_points[3];

	public:

		RealSimplex ( int point, Simplex &oS, std::vector<Eigen::Vector3d> &oX, Eigen::Vector3d sh );
		RealSimplex() {}
		RealSimplex &operator= ( const RealSimplex &RS ) { 
			point_nu = RS.point_nu;
			for (int i=0; i<3; i++) other_points[i] = RS.other_points[i];
			S=RS.S; rX=RS.rX; shift=RS.shift; Min=RS.Min; Max=RS.Max; 
			return *this; 
		}

		Eigen::Vector3d operator[] (int i) { return (*rX)[S->GetPoint(i)]+shift; }
		void GetPoints ( Eigen::Vector3d &r_nu, 
				 Eigen::Vector3d &r_mu, 
				 Eigen::Vector3d &r_gamma, 
				 Eigen::Vector3d &r_sigma  ) { 
			r_nu 	= (*rX)[point_nu]+shift;
			r_mu 	= (*rX)[other_points[0]]+shift;
			r_gamma	= (*rX)[other_points[1]]+shift;
			r_sigma	= (*rX)[other_points[2]]+shift;
		}

		void UpdateMinMax ( Eigen::Vector3d &mi, Eigen::Vector3d &ma ) {
			for (int i=0; i<3; i++) { 
				if (Min[i]<mi[i]) mi[i]=Min[i];
				if (Max[i]>ma[i]) ma[i]=Max[i];
			}
		}
		bool in_the_box ( Eigen::Vector3d &r ) {
			for (int i=0; i<3; i++)
				if ( r[i]<Min[i] ) return false;
				else if ( r[i]>Max[i] ) return false;
			return true;
		}

		void PrintSimplex();
		void PrintSimplex ( FILE *f );

};


class Tesselation;


inline double Heaviside ( double x ) { if (x<0) return 0; if (x>0) return 1; return 1.0/2.0; }


class Star {

		int point;
		
		int nRS;
		std::vector<RealSimplex> RS;

		Eigen::Vector3d r_nu, r_mu, r_gamma, r_sigma;
		Eigen::Vector3d n;

		Eigen::Vector3d Min, Max;

	public:

		Star() {}
		~Star() {}
		void Build ( int p, Delaunay &D, std::vector<Eigen::Vector3d> &X, Tesselation &tssl );

		int GetNSimplices() { return nRS; }
		RealSimplex &operator[] (int i) { return RS[i]; }

		bool in_the_box ( Eigen::Vector3d &r ) {
			for (int i=0; i<3; i++)
				if ( r[i]<Min[i] ) return false;
				else if ( r[i]>Max[i] ) return false;
			return true;
		}
		bool in_simplex_box ( int i, Eigen::Vector3d &r ) {
			return RS[i].in_the_box(r);
		}

		// theta_{nu_mu_gamma_delta}, nu - the centre of the cell
		int theta ( int i, Eigen::Vector3d &r ) {

			RealSimplex &S = RS[i];
			S.GetPoints ( r_nu, r_mu, r_gamma, r_sigma );

			if ( !theta_nu ( r, r_nu, r_mu, r_gamma, r_sigma ) ) return 0;
			if ( !theta_nu ( r, r_mu, r_nu, r_gamma, r_sigma ) ) return 0;
			if ( !theta_nu ( r, r_gamma, r_mu, r_nu, r_sigma ) ) return 0;
			if ( !theta_nu ( r, r_sigma, r_mu, r_gamma, r_nu ) ) return 0;

			return 1;

		}

		double t ( int i, Eigen::Vector3d &r ) {

			RealSimplex &S = RS[i];
			S.GetPoints ( r_nu, r_mu, r_gamma, r_sigma );

			double V_nu = 
				 (r_mu[0]-r_nu[0])*((r_gamma[1]-r_nu[1])*(r_sigma[2]-r_nu[2])-(r_sigma[1]-r_nu[1])*(r_gamma[2]-r_nu[2]))
				-(r_mu[1]-r_nu[1])*((r_gamma[0]-r_nu[0])*(r_sigma[2]-r_nu[2])-(r_sigma[0]-r_nu[0])*(r_gamma[2]-r_nu[2]))
				+(r_mu[2]-r_nu[2])*((r_gamma[0]-r_nu[0])*(r_sigma[1]-r_nu[1])-(r_sigma[0]-r_nu[0])*(r_gamma[1]-r_nu[1]));

			double V_r = 
				 (r_mu[0]-r[0])*((r_gamma[1]-r[1])*(r_sigma[2]-r[2])-(r_sigma[1]-r[1])*(r_gamma[2]-r[2]))
				-(r_mu[1]-r[1])*((r_gamma[0]-r[0])*(r_sigma[2]-r[2])-(r_sigma[0]-r[0])*(r_gamma[2]-r[2]))
				+(r_mu[2]-r[2])*((r_gamma[0]-r[0])*(r_sigma[1]-r[1])-(r_sigma[0]-r[0])*(r_gamma[1]-r[1]));

			return V_r/V_nu;

		}

		void b ( double *c_k, int i ) {

			RealSimplex &S = RS[i];
			S.GetPoints ( r_nu, r_mu, r_gamma, r_sigma );

			c_k[0] = (r_gamma[1]-r_mu[1])*(r_sigma[2]-r_mu[2])-(r_sigma[1]-r_mu[1])*(r_gamma[2]-r_mu[2]);
			c_k[1] = (r_sigma[0]-r_mu[0])*(r_gamma[2]-r_mu[2])-(r_gamma[0]-r_mu[0])*(r_sigma[2]-r_mu[2]);
			c_k[2] = (r_gamma[0]-r_mu[0])*(r_sigma[1]-r_mu[1])-(r_sigma[0]-r_mu[0])*(r_gamma[1]-r_mu[1]);

			c_k[0] /= 6.0*GetSimplexVolume(i);
			c_k[1] /= 6.0*GetSimplexVolume(i);
			c_k[2] /= 6.0*GetSimplexVolume(i);

			static Eigen::Vector3d g = r_nu - (r_mu+r_gamma+r_sigma)/3;
			static Eigen::Vector3d c(c_k);
			if ( g.dot(c) < 0 ) { c_k[0] *= -1; c_k[1] *= -1; c_k[2] *= -1; }

		}

		void PrintSimplex ( int i ) { RS[i].PrintSimplex(); }
		double GetSimplexVolume ( int i ) { 
			RS[i].GetPoints ( r_nu, r_mu, r_gamma, r_sigma ); 
			return 0.1666666666666*fabs( (r_mu[0]-r_nu[0])*((r_gamma[1]-r_nu[1])*(r_sigma[2]-r_nu[2])-(r_sigma[1]-r_nu[1])*(r_gamma[2]-r_nu[2]))
				-(r_mu[1]-r_nu[1])*((r_gamma[0]-r_nu[0])*(r_sigma[2]-r_nu[2])-(r_sigma[0]-r_nu[0])*(r_gamma[2]-r_nu[2]))
				+(r_mu[2]-r_nu[2])*((r_gamma[0]-r_nu[0])*(r_sigma[1]-r_nu[1])-(r_sigma[0]-r_nu[0])*(r_gamma[1]-r_nu[1])) );
		}
		void PrintAllSimplexes() { FILE* f = fopen("simplexes.dat", "a"); for (int i=0; i<nRS; i++) RS[i].PrintSimplex(f); fclose(f); }
		void PrintAllSimplexesVolumes() { FILE* f = fopen("simplexes_V.dat", "a"); for (int i=0; i<nRS; i++) fprintf ( f, "%f\n", GetSimplexVolume(i) ); fclose(f); }
		Eigen::Vector3d get_mu ( int i ) { RealSimplex &S = RS[i]; S.GetPoints ( r_nu, r_mu, r_gamma, r_sigma ); return r_mu; }
		Eigen::Vector3d get_gamma ( int i ) { RealSimplex &S = RS[i]; S.GetPoints ( r_nu, r_mu, r_gamma, r_sigma ); return r_gamma; }
		Eigen::Vector3d get_sigma ( int i ) { RealSimplex &S = RS[i]; S.GetPoints ( r_nu, r_mu, r_gamma, r_sigma ); return r_sigma; }

		double GetVolume() {
		
			double vol=0;
			for (int i=0; i<nRS; i++)
				vol += GetSimplexVolume(i);
				
			return vol;
		}

	private:

		int theta_nu ( 	Eigen::Vector3d &r, 
				Eigen::Vector3d &r_nu, 
				Eigen::Vector3d &r_mu, 
				Eigen::Vector3d &r_gamma, 
				Eigen::Vector3d &r_sigma	) {

			n[0] = (r_sigma[1]-r_mu[1])*(r_gamma[2]-r_mu[2]) - (r_gamma[1]-r_mu[1])*(r_sigma[2]-r_mu[2]);
			n[1] = (r_gamma[0]-r_mu[0])*(r_sigma[2]-r_mu[2]) - (r_sigma[0]-r_mu[0])*(r_gamma[2]-r_mu[2]);
			n[2] = (r_sigma[0]-r_mu[0])*(r_gamma[1]-r_mu[1]) - (r_gamma[0]-r_mu[0])*(r_sigma[1]-r_mu[1]);

			double lambda = (n.dot(r_mu-r_nu))/n.squaredNorm();

			return Heaviside ( (lambda*n).dot( r_nu + lambda*n - r ) );

		}

};


class Tesselation {

		std::vector<Eigen::Vector3d> X;
		Delaunay D;

		double Min[3], Max[3];
		int nx, ny, nz;

		std::vector<Star> St;
		friend class Star;

		int curr_star, found_simplex;
		Eigen::Vector3d r;

	public:

		Tesselation ( const std::vector<Eigen::Vector3d> &oX, double mi[3], double ma[3], int N[3] );
		~Tesselation(){}

		bool in_this_cell ( int k, double r[3] );
		double mu() {  // uses information from in_this_cell above, so should be called right after in_this_cell
			return St[curr_star].t(found_simplex,r); 
		}
//		void c ( double *c_k, Vector<double> &g, Vector<double> &centr ) {  // uses information from in_this_cell above, so should be called right after in_this_cell
		void c ( double *c_k ) {  // uses information from in_this_cell above, so should be called right after in_this_cell
			St[curr_star].b(c_k,found_simplex); 
		}

		double GetStarVolume(int cell) {
			return St[cell].GetVolume();
		}

		int GetNpoints() { return X.size(); }

		void PrintFoundSimplex() { St[curr_star].PrintSimplex(found_simplex); }
		double GetFoundSimplexVolume() { return St[curr_star].GetSimplexVolume(found_simplex); }
		Eigen::Vector3d get_mu() { return St[curr_star].get_mu(found_simplex); }
		Eigen::Vector3d get_gamma() { return St[curr_star].get_gamma(found_simplex); }
		Eigen::Vector3d get_sigma() { return St[curr_star].get_sigma(found_simplex); }

		inline bool IsMinX ( int point ) { return X[point][0]==Min[0]; }
		inline bool IsMaxX ( int point ) { return X[point][0]==Max[0]; }
		inline bool IsMinY ( int point ) { return X[point][1]==Min[1]; }
		inline bool IsMaxY ( int point ) { return X[point][1]==Max[1]; }
		inline bool IsMinZ ( int point ) { return X[point][2]==Min[2]; }
		inline bool IsMaxZ ( int point ) { return X[point][2]==Max[2]; }

		// sides    !!! tested only for nx=ny=nz
		inline int PeriodicMinX ( int point ) { return point+nx*ny*nz-nx*ny; }
		inline int PeriodicMaxX ( int point ) { return point-nx*ny*nz+ny*nz; }
		
		inline int PeriodicMinY ( int point ) { return point+nx*nz-ny; }
		inline int PeriodicMaxY ( int point ) { return point-nx*nz+ny; }

		inline int PeriodicMinZ ( int point ) { return point+nz-1; }
		inline int PeriodicMaxZ ( int point ) { return point-nz+1; }

		// edges    !!! tested only for nx=ny=nz
		inline int PeriodicMinXMinY ( int point ) { return point+nx*ny*nz-nz; }
		inline int PeriodicMinXMaxY ( int point ) { return point+nx*ny*nz-nx*ny-nx*ny+nz; }
		inline int PeriodicMaxXMinY ( int point ) { return point-nx*ny*nz+nx*ny+nx*ny-nz; }
		inline int PeriodicMaxXMaxY ( int point ) { return point-nx*ny*nz+nz; }
		
                inline int PeriodicMinXMinZ ( int point ) { return point+nx*ny*nz-nx*nz+ny-1; }
		inline int PeriodicMinXMaxZ ( int point ) { return point+nx*ny*nz-nx*nz-ny+1; }
		inline int PeriodicMaxXMinZ ( int point ) { return point-nx*ny*nz+nx*nz+ny-1; }
		inline int PeriodicMaxXMaxZ ( int point ) { return point-nx*ny*nz+nx*nz-ny+1; }
		
		inline int PeriodicMinYMinZ ( int point ) { return point+ny*nz-1; }
		inline int PeriodicMinYMaxZ ( int point ) { return point+ny*nz-nx-nx+1; }
                inline int PeriodicMaxYMinZ ( int point ) { return point-ny*nz+nx+nx-1; }
		inline int PeriodicMaxYMaxZ ( int point ) { return point-ny*nz+1; }
		
		// corners    !!! tested only for nx=ny=nz
                inline int PeriodicMinXMinYMinZ ( int point ) { return point+nx*ny*nz-1; }
		inline int PeriodicMinXMinYMaxZ ( int point ) { return point+nx*ny*nz-ny-nz+1; }
		inline int PeriodicMinXMaxYMinZ ( int point ) { return point+nx*ny*nz-nx*ny-ny*nz+ny+nz-1; }
		inline int PeriodicMinXMaxYMaxZ ( int point ) { return point+nx*ny*nz-nx*ny-ny*nz+1; }
		
		inline int PeriodicMaxXMinYMinZ ( int point ) { return point-nx*ny*nz+nx*ny+ny*nz-1; }
		inline int PeriodicMaxXMinYMaxZ ( int point ) { return point-nx*ny*nz+nx*ny+ny*nz-nx-ny+1; }
                inline int PeriodicMaxXMaxYMinZ ( int point ) { return point-nx*ny*nz+nx+ny-1; }
		inline int PeriodicMaxXMaxYMaxZ ( int point ) { return point-nx*ny*nz+1; }

};



#endif
