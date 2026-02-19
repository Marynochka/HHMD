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

#define max_blocks 10000
#define Nobj 10000


#define NEWL(X,p)						\
{								\
 	p = X##_list ? X##_list : new_block_##X(1);		\
	assert(p);						\
 	X##_list = p->next;					\
}								



#define NEWLRC(X,p)						\
{								\
	p = X##_list ? X##_list : new_block_##X(1);		\
	assert(p);						\
	X##_list = p->next;					\
	p->ref_count = 1;					\
}								\


#define FREEL(X,p)						\
{								\
	memset((p),0,X##_size);					\
	(p)->next = X##_list;					\
	X##_list = p;						\
}								


#define dec_ref(X,v)	{if ((v) && --(v)->ref_count == 0) FREEL(X,(v));}
#define inc_ref(X,v)	{if (v) v->ref_count++;}
#define NULLIFY(X,v)	{if ((v) && --(v)->ref_count == 0) FREEL(X,(v)); v = NULL;}



#define mod_refs(op,s)					\
{							\
	int i;						\
	h_neighbor *mrsn;					\
							\
	for (i=-1,mrsn=s->neigh-1;i<cdim;i++,mrsn++)	\
		op##_ref(h_basis_s, mrsn->basis);		\
}


#define copy_simp(h_new,s)			\
{	NEWL(h_simplex,h_new);			\
	memcpy(h_new,s,h_simplex_size);		\
	mod_refs(inc,s);			\
}						\



