#include <R.h>

void risk_set_etm(int *n, int *lt, int *dim_nev, double *times,
		  int *from, int *to, double *entry, double *exit,
		  int *nrisk, int *ncens, int *nev, double *dna) {
	
    const int ltimes = *lt;
    const int dim_trans = dim_nev[1];
    const int nb = *n;
    int i, j, t;
    
    /* Computation of the risk set and transition matrix */
    
    for (j=0; j < nb; ++j) {
	for (i=0; i < ltimes; ++i) {
	    if (entry[j] < times[i] && exit[j] >= times[i]) {
		nrisk[i + *lt * (from[j] - 1)] += 1;
	    }
	    if (exit[j] == times[i] && to[j] != 0) {
		nev[dim_nev[1] * dim_nev[1]*i + from[j] - 1 + dim_nev[1] * (to[j] - 1)] += 1;
		break;
	    }
	}
    }
    
    for (i = 0; i < dim_trans; ++i) {
	nrisk[i * (*lt)] = nrisk[i * (*lt) + 1];
    }
	
    /* Nelson-Aalen increments */
    for (t = 0; t < ltimes; ++t) {
	for (j = 0; j < dim_trans; ++j) {
	    for (i = 0; i < dim_trans; ++i) {
		if (nrisk[i * (*lt) + t] != 0) {
			dna[dim_nev[1] * dim_nev[1] * t + j * dim_nev[1] + i] =
			    (double) nev[dim_nev[1] * dim_nev[1] * t + j * dim_nev[1]+i] / (double) nrisk[i*(*lt)+t];
		}
	    }
	}
    }
}

