#ifndef BINOMIALS_H_INCLUDED
#define BINOMIALS_H_INCLUDED

extern gdouble _BINOMIALS[] ;
#define _binomial(_m,_k) (_BINOMIALS[(_m)*((_m)+1)/2+(_k)])

#define _binomial_list(_m) &(_BINOMIALS[(_m)*((_m)+1)/2])

#endif /*BINOMIALS_H_INCLUDED*/
