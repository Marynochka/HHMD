/* hull.c : "combinatorial" functions for hull computation */

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <assert.h>
#include <string.h>

#include "hull.h"


MemGuard memory;

#define push(x) *(sear_st+tms++) = x;
#define pop(x)  x = *(sear_st + --tms);

#define vi_push(x) *(vi_st+tms++) = x;
#define vi_pop(x)  x = *(vi_st + --tms);


#define lookup(a,b,what,whatt)						\
{									\
	int i;								\
	h_neighbor *x;							\
	for (i=0, x = a->neigh; (x->what != b) && (i<cdim) ; i++, x++);	\
	if (i<cdim)							\
		return x;						\
	else {								\
/*		fprintf(DFILE,"adjacency failure,op_" #what ":\n");	*/\
/*		DEBTR(-10)					*/	\
/*		print_simplex_f(a, DFILE, &print_neighbor_full);	*/\
/*		print_##whatt(b, DFILE);			*/	\
/*		fprintf(DFILE,"---------------------\n");	*/	\
/*		print_triang(a,DFILE, &print_neighbor_full);	*/	\
		exit(1);						\
		return 0;						\
	}								\
}									\


h_neighbor *ClarksonHull::op_simp(h_simplex *a, h_simplex *b) {lookup(a,b,simp,simplex)}
	/* the neighbor entry of a containing b */

h_neighbor *ClarksonHull::op_vert(h_simplex *a, h_site b) {lookup(a,b,vert,site)}
	/* the neighbor entry of a containing b */


void ClarksonHull::connect(h_simplex *s) {
/* make neighbor connections between newly created simplices incident to p */

	h_site xf,xb,xfi;
	h_simplex *sb, *sf, *seen;
	int i;
	h_neighbor *sn;

	if (!s) return;
	assert(!s->peak.vert
		&& s->peak.simp->peak.vert==p
		&& !op_vert(s,p)->simp->peak.vert);
	if (s->visit==pnum) return;
	s->visit = pnum;
	seen = s->peak.simp;
	xfi = op_simp(seen,s)->vert;
	for (i=0, sn = s->neigh; i<cdim; i++,sn++) {
		xb = sn->vert;
		if (p == xb) continue;
		sb = seen;
		sf = sn->simp;
		xf = xfi;
		if (!sf->peak.vert) {	/* are we done already? */
			sf = op_vert(seen,xb)->simp;
			if (sf->peak.vert) continue;				
		} else do {
			xb = xf;
			xf = op_simp(sf,sb)->vert;
			sb = sf;
			sf = op_vert(sb,xb)->simp;
		} while (sf->peak.vert);

		sn->simp = sf;
		op_vert(sf,xf)->simp = s;

		connect(sf);
	}

}


				
h_simplex *ClarksonHull::make_facets(h_simplex *seen) {
/*
 * visit simplices s with sees(p,s), and make a facet for every neighbor
 * of s not seen by p
 */

	h_simplex *n;
	h_neighbor *bn;
	int i;

	if (!seen) return NULL;
	seen->peak.vert = p;

	for (i=0,bn = seen->neigh; i<cdim; i++,bn++) {
		n = bn->simp;
		if (pnum != n->visit) {
			n->visit = pnum;
			if (sees(p,n)) make_facets(n);
		} 
		if (n->peak.vert) continue;
		copy_simp(ma_ns,seen);
		ma_ns->visit = 0;
		ma_ns->peak.vert = 0;
		ma_ns->normal = 0;
		ma_ns->peak.simp = seen;
/*		ns->Sb -= ns->neigh[i].basis->sqb; */
		NULLIFY(h_basis_s,ma_ns->neigh[i].basis);
		ma_ns->neigh[i].vert = p;
		bn->simp = op_simp(n,seen)->simp = ma_ns;
	}
	return ma_ns;
}



h_simplex *ClarksonHull::extend_simplices(h_simplex *s) {
/*
 * p lies outside flat containing previous sites;
 * make p a vertex of every current simplex, and create some new simplices
 */

	int	i,
		ocdim=cdim-1;
	h_simplex *ns;
	h_neighbor *nsn;

	if (s->visit == pnum) return s->peak.vert ? s->neigh[ocdim].simp : s;
	s->visit = pnum;
	s->neigh[ocdim].vert = p;
	NULLIFY(h_basis_s,s->normal);
	NULLIFY(h_basis_s,s->neigh[0].basis);
	if (!s->peak.vert) {
		s->neigh[ocdim].simp = extend_simplices(s->peak.simp);
		return s;
	} else {
		copy_simp(ns,s);
		s->neigh[ocdim].simp = ns;
		ns->peak.vert = NULL;
		ns->peak.simp = s;
		ns->neigh[ocdim] = s->peak;
		inc_ref(h_basis_s,s->peak.basis);
		for (i=0,nsn=ns->neigh;i<cdim;i++,nsn++)
			nsn->simp = extend_simplices(nsn->simp);
	}
	return ns;
}


h_simplex *ClarksonHull::search(h_simplex *root) {
// return a simplex s that corresponds to a facet of the 
// current hull, and sees(p, s) 

	h_simplex *s;
	h_neighbor *sn;
	int i;
	long tms = 0;

//	if (!st) st = (h_simplex **)malloc((ss+MAXDIM+1)*sizeof(h_simplex*));
	if (!sear_st) sear_st = (h_simplex **)memory((sear_ss+MAXDIM+1)*sizeof(h_simplex*));
	push(root->peak.simp);
	root->visit = pnum;
	if (!sees(p,root))
		for (i=0,sn=root->neigh;i<cdim;i++,sn++) push(sn->simp);
	while (tms) {
		if (tms>sear_ss) {
			sear_st=(h_simplex**)memory.realloc(sear_st,
				((sear_ss+=sear_ss)+MAXDIM+1)*sizeof(h_simplex*));
			assert(sear_st);
		}
		pop(s);
		if (s->visit == pnum) continue;
		s->visit = pnum;
		if (!sees(p,s)) continue;
		if (!s->peak.vert) return s;
		for (i=0, sn=s->neigh; i<cdim; i++,sn++) push(sn->simp);
	}
	return NULL;
}


h_point ClarksonHull::get_another_site(void) {

	h_point pnext;

	pnext = get_next_site();
	if (!pnext) return NULL;
	pnum = site_numm(pnext)+2;
	return pnext;
}



void ClarksonHull::buildhull (h_simplex *root) {

	while (cdim < rdim) {
		p = get_another_site();
		if (!p) return;
		if (out_of_flat(root,p))
			extend_simplices(root);
		else
			connect(make_facets(search(root)));
	}
	while ((p = get_another_site()))
		connect(make_facets(search(root)));
}


h_site ClarksonHull::new_site (h_site p, long j) {

	assert(num_blocks+1<H_MAXBLOCKS);
	if (0==(j%H_BLOCKSIZE)) {
		assert(num_blocks < H_MAXBLOCKS);
		return(site_blocks[num_blocks++]=(h_site)memory(H_BLOCKSIZE*h_site_size));
	} else
		return p+dim;
}


h_site ClarksonHull::read_next_site(long j){

	int i=0;//, k=0;

	if (!dim) dim = X[0].size();
	if (j==-1) return 0;

	if (j >= npoints) return 0;

	p = new_site(p,j);

	for (i=0; i<dim; i++) {
		p[i] = double(X[j][i]);
		p[i] = floor(mult_up*p[i]+0.5);   
		mins[i] = (mins[i]<p[i]) ? mins[i] : p[i];
		maxs[i] = (maxs[i]>p[i]) ? maxs[i] : p[i];
	}

	return p;
}


h_site ClarksonHull::get_next_site(void) {
	return read_next_site(s_num++);
}


long ClarksonHull::site_numm(h_site p) {

	long i,j;

	if (p==infinity) return -1;
	if (!p) return -2;
	for (i=0; i<num_blocks; i++)
		if ((j=p-site_blocks[i])>=0 && j < H_BLOCKSIZE*dim) 
			return j/dim+H_BLOCKSIZE*i;
	return -3;
}


typedef void out_func(h_point *, int, int);
typedef h_simplex* visit_func(h_simplex *, int, out_func);
typedef int test_func(h_simplex *, int, void *);

typedef int test_func(h_simplex *, int, void *);

int hullt(h_simplex *s, int i, void *dummy) {return i>-1;}
int truet(h_simplex *s, int i, void *dum) {return 1;}
void *facet_test(h_simplex *s, void *dummy) {return (!s->peak.vert) ? s : NULL;}

void ClarksonHull::vlist_out(h_point *v, int vdim, int Fin) {

//	static int F;
	int sn;
	int j;

	if (Fin) {/*F=Fin;*/ if (!v) return;}

	if ( !dooutput ) {
		for (j=0;j<vdim;j++)  {
			sn = site_numm(v[j]);
			if ( sn == -1 ) return;
		}
		nsimplices++;
		return;
	}

	unsigned pts[MAXDIM];
	double minbox[MAXDIM], maxbox[MAXDIM];

	for (j=0; j<vdim; j++) { minbox[j] = 1e32; maxbox[j] = -1e32; }
	for (j=0;j<vdim;j++)  {
		pts[j] = sn = site_numm(v[j]);
		if ( sn == -1 ) return;
		for (int k=0; k<dim; k++) {
			if ( minbox[k] > X[sn][k] ) minbox[k] = (double)X[sn][k];
			if ( maxbox[k] < X[sn][k] ) maxbox[k] = (double)X[sn][k];
		}
	}
	
	Smplx.Add(pts,minbox,maxbox);

}


ClarksonHull *ClHu;

void vlist_out(h_point *v, int vdim, int Fin) {
	ClHu->vlist_out(v, vdim, Fin);
}


h_simplex *facets_print(h_simplex *s, int cdim, out_func *p) {

	static out_func *out_func_here;
	h_point v[MAXDIM];
	int j;

	if (p) {out_func_here = p; if (!s) return NULL;}

	for (j=0;j<cdim;j++) v[j] = s->neigh[j].vert;

	out_func_here(v,cdim,0);

	return NULL;
}

h_simplex *facet_test(h_simplex *s, int dum, out_func *dum1) {return (!s->peak.vert) ? s : NULL;}

h_simplex **vi_st;
long vnum = -1;
long vi_ss = 2000;

h_simplex *visit_triang_gen(h_simplex *s, visit_func *visit, test_func *test, int cdim) {
/* 
 * starting at s, visit simplices t such that test(s,i,0) is true,
 * and t is the i'th neighbor of s;
 * apply visit function to all visited simplices;
 * when visit returns nonNULL, exit and return its value
 */
	h_neighbor *sn;
	h_simplex *v;
	h_simplex *t;
	int i;
	long tms = 0;

	vnum--;
	if (!vi_st) { vi_st=(h_simplex**)memory((vi_ss+H_MAXDIM+1)*sizeof(h_simplex*)); assert(vi_st); }
	if (s) vi_push(s);
	while (tms) {

		if (tms>vi_ss) {//DEBEXP(-1,tms);
			vi_st=(h_simplex**)memory.realloc(vi_st,
				((vi_ss+=vi_ss)+H_MAXDIM+1)*sizeof(h_simplex*));
			assert(vi_st);
		}
		vi_pop(t);
		if (!t || t->visit == vnum) continue;
		t->visit = vnum;
		if ((v=(*visit)(t,cdim,0))) {return v;}
		for (i=-1,sn = t->neigh-1;i<cdim;i++,sn++)
			if ((sn->simp->visit != vnum) && sn->simp && test(t,i,0))
				vi_push(sn->simp);
	}
	return NULL;
}

h_simplex *visit_triang(h_simplex *root, visit_func *visit, int cdim)
/* visit the whole triangulation */
	{return visit_triang_gen(root, visit, truet, cdim);}


void ClarksonHull::make_output(h_simplex *root){

	ClHu = this;
	facets_print(0, cdim, ::vlist_out);
	visit_triang_gen(visit_triang(root, facet_test, cdim), facets_print, hullt, cdim);

}


