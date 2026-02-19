/* hull.h */
/*
 * Ken Clarkson wrote this.  Copyright (c) 1995 by AT&T..
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 */

#ifndef HDR
#define HDR 1

#include <float.h>
#include <assert.h>

typedef double h_Coord;
typedef h_Coord* h_point;

#include "stormacs.h"
#include "Delaunay.h"


#define H_MAXDIM 8
#define H_BLOCKSIZE 10000
#define H_MAXBLOCKS 1000


typedef h_point h_site;


typedef struct h_basis_s {
	struct h_basis_s *next; /* free list */
	int ref_count;	/* storage management */
	int lscale;    /* the log base 2 of total scaling of vector */
	h_Coord sqa, sqb; /* sums of squared norms of a part and b part */
	h_Coord vecs[1]; /* the actual vectors, extended by malloc'ing bigger */
} h_basis_s;

typedef struct h_neighbor {
	h_site vert; /* vertex of simplex */
	struct h_simplex *simp; /* neighbor sharing all vertices but vert */
	h_basis_s *basis; /* derived vectors */
} h_neighbor;

typedef struct h_simplex {
	struct h_simplex *next;	/* free list */
	long visit;		/* number of last site visiting this simplex */
/*	float Sb; */
	short mark;
	h_basis_s* normal;	/* normal vector pointing inward */
	h_neighbor peak;		/* if null, remaining vertices give facet */
	h_neighbor neigh[1];	/* neighbors of simplex */
} h_simplex;

typedef struct h_fg_node h_fg;
typedef struct h_tree_node h_Tree;
struct h_tree_node {
    h_Tree *left, *right;
    h_site key;
    int size;   /* maintained to be the number of nodes rooted here */
    h_fg *fgs;
    h_Tree *next; /* freelist */
};

typedef struct h_fg_node {
	h_Tree *facets;
	double dist, vol;	/* of Voronoi face dual to this */
	h_fg *next;  		/* freelist */
	short mark;
	int ref_count;
} h_fg_node;
	

#define H_MAXPOINTS 10000

/* from hull.c */

#define max_blocks 10000
#define Nobj 10000


#define MG_BLOCKSIZE	1000
#define MG_MAXBLOCKS	100

class MemGuard {

		int npointers;
		int nblocks;
		void **body[MG_MAXBLOCKS];
		void **curr;

		void Assert ( int b ) {
			if ( !b ) {
				printf ( "\nerror in memguard\n" );
				exit(1);
			}
		}

		void AddPointer ( void *p ) {
			
			Assert ( nblocks+1 < MG_MAXBLOCKS );
			int n = npointers++ % MG_BLOCKSIZE;
			if ( !n ) {
				Assert ( nblocks < MG_MAXBLOCKS );
				curr = body[nblocks++]=(void**)malloc(MG_BLOCKSIZE*sizeof(void*));
			} 

			curr[n] = p;
			
		}

		void ReplacePointer ( void *prev, void *pnew ) {

			for (int i=0; i<nblocks; i++) {
				void **p = body[i];
				int j;
				for (j=0; j<MG_BLOCKSIZE && i*MG_BLOCKSIZE+j < npointers; j++) 
					if ( p[j] == prev ) {
						p[j] = pnew;
						return;
					}
				if ( i*MG_BLOCKSIZE+j >= npointers ) break;
			}

		}

	public:

		MemGuard() : npointers(0), nblocks(0) {}

		void *operator() (size_t size) {
			void *p = malloc(size);
			if ( p ) AddPointer(p);
			return p;
		}

		void *realloc ( void *prev, size_t size ) {

			void *pnew = ::realloc ( prev, size );
			if ( pnew ) ReplacePointer ( prev, pnew );
			return pnew;

		}

		void Free() {

			int i,j;
			for (i=0; i<nblocks; i++) {
				void **p = body[i];
				for (j=0; j<MG_BLOCKSIZE && i*MG_BLOCKSIZE+j < npointers; j++) 
					free ( p[j] );
				if ( i*MG_BLOCKSIZE+j >= npointers ) break;
			}

			for (i=0; i<nblocks; i++)
				free(body[i]);

			npointers = nblocks = 0;

		}

};


extern MemGuard memory;
extern h_simplex **vi_st;
extern long vnum;
extern long vi_ss;

class ClarksonHull {

		double mult_up;
		const std::vector<Eigen::Vector3d> &X;
		int dim;
		int npoints;
		double Huge;
		h_basis_s tt_basis, *tt_basisp, *infinity_basis;
		int num_blocks;
		h_basis_s *h_basis_s_list;
		h_simplex *h_simplex_list;
		long s_num;
		int	sc_lscale;
		double	sc_max_scale, sc_ldetbound, sc_Sb;
		h_basis_s *see_b;
		h_simplex *ma_ns;
		h_simplex **sear_st;
		long sear_ss;
		SimplexList &Smplx;

		h_simplex *root;
		int	pdim;				/* point dimension */
		h_site p;				/* the current site */
		h_Coord infinity[10];	/* point at infinity for Delaunay triang */
		int	rdim,	/* region dimension: (max) number of sites specifying region */
			cdim,	/* number of sites currently specifying region */
			h_site_size, /* size of malloc needed for a site */
			h_point_size;  /* size of malloc needed for a point */
		h_Coord mins[H_MAXDIM], maxs[H_MAXDIM];
		h_point site_blocks[H_MAXBLOCKS];
		int exact_bits;
		float b_err_min, b_err_min_sq;
		int h_basis_vec_size;
		size_t h_basis_s_size;
		size_t h_simplex_size;
		size_t h_Tree_size;
		h_Tree *h_Tree_list;
		size_t h_fg_size;
		h_fg *h_fg_list;
		bool dooutput;
		unsigned nsimplices;
		h_simplex *ch_root;
		long pnum;
		h_neighbor out_p_neigh;



		long site_numm(h_site p);
		h_simplex *build_convex_hull(short dim);
		h_site read_next_site(long j);
		h_site get_next_site(void);
		h_site new_site (h_site p, long j);
		h_basis_s *new_block_h_basis_s(int make_blocks) {	
			int i;
			static	h_basis_s *h_basis_s_block_table[max_blocks];
			h_basis_s *xlm, *xbt;
			static int num_h_basis_s_blocks;
			/* long before;	*/
			if (make_blocks) {
				/* DEBEXP(-10, num_##X##_blocks) */
				assert(num_h_basis_s_blocks<max_blocks);
				/* before = _memfree(0);*/
				xbt = h_basis_s_block_table[num_h_basis_s_blocks++] =
					(h_basis_s*)memory(Nobj * h_basis_s_size);
 				memset(xbt,0,Nobj * h_basis_s_size);
				assert(xbt);

				xlm = (h_basis_s*)((char*)xbt + (Nobj) * h_basis_s_size);
				for (i=0;i<Nobj; i++) {
					xlm = (h_basis_s*)((char*)xlm + (-1) * h_basis_s_size);
					xlm->next = h_basis_s_list;
					h_basis_s_list = xlm;
				}

				return h_basis_s_list;
			};

//			for (i=0; i<num_h_basis_s_blocks; i++)
//				free(h_basis_s_block_table[i]);
			num_h_basis_s_blocks = 0;
			h_basis_s_list = 0;
			return 0;
		}
		h_simplex *new_block_h_simplex(int make_blocks) {	
			int i;
			static	h_simplex *h_simplex_block_table[max_blocks];
			h_simplex *xlm, *xbt;
			static int num_h_simplex_blocks;
			/* long before;	*/
			if (make_blocks) {
				/* DEBEXP(-10, num_##X##_blocks) */
				assert(num_h_simplex_blocks<max_blocks);
				/* before = _memfree(0);*/
				xbt = h_simplex_block_table[num_h_simplex_blocks++] =
					(h_simplex*)memory(Nobj * h_simplex_size);
 				memset(xbt,0,Nobj * h_simplex_size);
				assert(xbt);

				xlm = (h_simplex*)((char*)xbt + (Nobj) * h_simplex_size);
				for (i=0;i<Nobj; i++) {
					xlm = (h_simplex*)((char*)xlm + (-1) * h_simplex_size);
					xlm->next = h_simplex_list;
					h_simplex_list = xlm;
				}

				return h_simplex_list;
			};

//			for (i=0; i<num_h_simplex_blocks; i++)
//				free(h_simplex_block_table[i]);
			num_h_simplex_blocks = 0;
			h_simplex_list = 0;
			return 0;
		}
		h_Tree *new_block_h_Tree(int make_blocks) {	
			int i;
			static	h_Tree *h_Tree_block_table[max_blocks];
			h_Tree *xlm, *xbt;
			static int num_h_Tree_blocks;
			/* long before;	*/
			if (make_blocks) {
				/* DEBEXP(-10, num_##X##_blocks) */
				assert(num_h_Tree_blocks<max_blocks);
				/* before = _memfree(0);*/
				xbt = h_Tree_block_table[num_h_Tree_blocks++] =
					(h_Tree*)memory(Nobj * h_Tree_size);
 				memset(xbt,0,Nobj * h_Tree_size);
				assert(xbt);

				xlm = (h_Tree*)((char*)xbt + (Nobj) * h_Tree_size);
				for (i=0;i<Nobj; i++) {
					xlm = (h_Tree*)((char*)xlm + (-1) * h_Tree_size);
					xlm->next = h_Tree_list;
					h_Tree_list = xlm;
				}

				return h_Tree_list;
			};

//			for (i=0; i<num_h_Tree_blocks; i++)
//				free(h_Tree_block_table[i]);
			num_h_Tree_blocks = 0;
			h_Tree_list = 0;
			return 0;
		}
		h_fg *new_block_h_fg(int make_blocks) {	
			int i;
			static	h_fg *h_fg_block_table[max_blocks];
			h_fg *xlm, *xbt;
			static int num_h_fg_blocks;
			/* long before;	*/
			if (make_blocks) {
				/* DEBEXP(-10, num_##X##_blocks) */
				assert(num_h_fg_blocks<max_blocks);
				/* before = _memfree(0);*/
				xbt = h_fg_block_table[num_h_fg_blocks++] =
					(h_fg*)memory(Nobj * h_fg_size);
 				memset(xbt,0,Nobj * h_fg_size);
				assert(xbt);

				xlm = (h_fg*)((char*)xbt + (Nobj) * h_fg_size);
				for (i=0;i<Nobj; i++) {
					xlm = (h_fg*)((char*)xlm + (-1) * h_fg_size);
					xlm->next = h_fg_list;
					h_fg_list = xlm;
				}

				return h_fg_list;
			};

//			for (i=0; i<num_h_fg_blocks; i++)
//				free(h_fg_block_table[i]);
			num_h_fg_blocks = 0;
			h_fg_list = 0;
			return 0;
		}
		void free_h_basis_s_storage(void) {new_block_h_basis_s(0);}
		void free_h_simplex_storage(void) {new_block_h_simplex(0);}
		void free_h_Tree_storage(void) {new_block_h_Tree(0);}
		void free_h_fg_storage(void) {new_block_h_fg(0);}
		int out_of_flat(h_simplex *root, h_point p);
		void get_basis_sede(h_simplex *s);
		int reduce(h_basis_s **v, h_point p, h_simplex *s, int k);
		h_Coord Vec_dot(h_point x, h_point y);
		int reduce_inner(h_basis_s *v, h_simplex *s, int k);
		h_Coord Norm2(h_point x);
		void Ax_plus_y(h_Coord a, h_point x, h_point y);
		double sc(h_basis_s *v,h_simplex *s, int k, int j);
		void Vec_scale_test(int n, h_Coord a, h_Coord *x);
		void Ax_plus_y_test(h_Coord a, h_point x, h_point y);
		h_simplex *extend_simplices(h_simplex *s);
		void connect(h_simplex *s);
		h_neighbor *op_vert(h_simplex *a, h_site b);
		h_neighbor *op_simp(h_simplex *a, h_simplex *b);
		h_point get_another_site(void);
		h_simplex *make_facets(h_simplex *seen);
		void make_output(h_simplex *root);
		void free_hull_storage(void);
		h_Coord Vec_dot_pdim(h_point x, h_point y);
		void buildhull (h_simplex *root);
		int sees(h_site p, h_simplex *s);
		void get_normal_sede(h_simplex *s);
		h_simplex *search(h_simplex *root);

	public:

		ClarksonHull ( const std::vector<Eigen::Vector3d> &oX, int np, SimplexList &oSmplx ) :
			mult_up(1e6), X(oX), dim(0), npoints(np),
			Huge(0), tt_basisp(&tt_basis), num_blocks(0),
			h_basis_s_list(0), h_simplex_list(0),
			s_num(0), sc_lscale(0),
			sc_max_scale(0), sc_ldetbound(0), sc_Sb(0),
			see_b(NULL), ma_ns(0), sear_st(0), sear_ss(MAXDIM),
			Smplx(oSmplx) {

			vi_st = 0;
			vnum = -1;
			vi_ss = 2000;

			out_p_neigh.vert=0;
			out_p_neigh.simp=0;
			out_p_neigh.basis=0;

			tt_basis.next=0;
			tt_basis.ref_count=1;
			tt_basis.lscale=-1;
			tt_basis.sqa=0;
			tt_basis.sqb=0;
			tt_basis.vecs[0]=0;

			for (int i=0; i<H_MAXDIM; i++ ) {
				mins[i] = DBL_MAX;
				maxs[i] = -DBL_MAX;
			}

			read_next_site(-1);

			h_point_size = h_site_size = sizeof(h_Coord)*dim;

			root = build_convex_hull( dim );

			dooutput = false;
			nsimplices = 0;
			make_output(root);
			dooutput = true;
//			nsimplices = 100000;

		}

		void Output() {
			make_output(root);
		}

		~ClarksonHull() {
			memory.Free();
//			free_hull_storage();
		}

		int GetNSimpices() { return nsimplices; }

		void vlist_out(h_point *v, int vdim, int Fin);

};




#endif
