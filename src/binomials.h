#ifndef _BINOMIALS_H_INCLUDED_
#define _BINOMIALS_H_INCLUDED_

extern gdouble _BINOMIALS[] ;
#define _binomial(_m,_k) (_BINOMIALS[(_m)*((_m)+1)/2+(_k)])

#define _binomial_list(_m) &(_BINOMIALS[(_m)*((_m)+1)/2])

#endif
