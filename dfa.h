/* file: dfa.c	  J. Mietus, C-K Peng, and G. Moody	8 February 2001
		  Last revised:			        25 January 2005  v4.9

-------------------------------------------------------------------------------
dfa: Detrended Fluctuation Analysis (translated from C-K Peng's Fortran code)
Copyright (C) 2001-2005 Joe Mietus, C-K Peng, and George B. Moody

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program; if not, write to the Free Software Foundation, Inc., 59 Temple
Place - Suite 330, Boston, MA 02111-1307, USA.

You may contact the authors by e-mail (peng@physionet.org) or postal mail
(Beth Israel Deaconess Medical Center, Room KS-B26, 330 Brookline Ave., Boston,
MA 02215 USA).  For updates to this software, please visit PhysioNet
(http://www.physionet.org/).
_______________________________________________________________________________

This method was first proposed in:
  Peng C-K, Buldyrev SV, Havlin S, Simons M, Stanley HE, Goldberger AL. Mosaic
  organization of DNA nucleotides. Phys Rev E 1994;49:1685-1689.  [Available
  on-line at http://prola.aps.org/abstract/PRE/v49/i2/p1685_1]

A detailed description of the algorithm and its application to physiologic
signals can be found in:
  Peng C-K, Havlin S, Stanley HE, Goldberger AL. Quantification of scaling
  exponents and crossover phenomena in nonstationary heartbeat time series.
  Chaos 1995;5:82-87. [Abstract online at http://www.ncbi.nlm.nih.gov/entrez/-
   query.fcgi?cmd=Retrieve&db=PubMed&list_uids=11538314&dopt=Abstract]

If you use this program in support of published research, please include a
citation of at least one of the two references above, as well as the standard
citation for PhysioNet:
  Goldberger AL, Amaral LAN, Glass L, Hausdorff JM, Ivanov PCh, Mark RG,
  Mietus JE, Moody GB, Peng CK, Stanley HE.  PhysioBank, PhysioToolkit, and
  Physionet: Components of a New Research Resource for Complex Physiologic
  Signals. Circulation 101(23):e215-e220 [Circulation Electronic Pages;
  http://circ.ahajournals.org/cgi/content/full/101/23/e215]; 2000 (June 13). 
  
  Adapted by Thalita B. Veronese 2011 (March 26)
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fitab.h"

#define SWAP(a,b) {temp = (a); (a) = (b); (b) = temp;}

/* Function prototypes. */
int rscale(VecDoub &rs, long minbox, long maxbox, double boxratio);
Doub dfa_call(VecDoub_I &A);
void dfa(VecDoub seq, long npts, int nfit, VecDoub &rs, int nr, int sw, VecDoub &mse, MatDoub &x, VecDoub &betha, MatDoub &covar, MatDoub &covar0, VecInt &ipiv, VecInt &indxc, VecInt &indxr);
void setup(int nfit, VecDoub &rs, int nr, VecDoub &mse, MatDoub &x);
void cleanup(void);
void help(void);
double polyfit(MatDoub &x, double *y, long ndat, int nfit, VecDoub &betha, MatDoub &covar, MatDoub &covar0, VecInt &ipiv, VecInt &indxc, VecInt &indxr);
double media(VecDoub_I &A);
void transp(MatDoub_I &a, MatDoub_O &at);
void diag(VecDoub_I &w, MatDoub_O &dw);
void mmultiply(MatDoub_I &a, MatDoub_I &b, MatDoub_O &c);

Doub dfa_call(VecDoub_I &A)
{
    
    VecDoub seq, rs, mse;
    int nr, nfit = 2;
    Doub alpha = 0;
    
    int sw = 0;
    long minbox = 0L, maxbox = 0L, npts, temp;
    seq.resize(A.size()+1);
    seq[0] = 0;
//    cout << media(A) << endl;
    for (int i=1; i<A.size(); i++) 
//        seq[i] = A[i-1];
//        seq[i] = A[i-1] - media(A);
          seq[i] = seq[i-1] + A[i-1] - media(A);
//          seq[i] = seq[i-1] + A[i-1];
    npts = A.size();
//    cout << "npts = " << npts << endl;
    

    /* Set minimum and maximum box sizes. */
    if (minbox < 2*nfit) minbox = 2*nfit;
    if (maxbox == 0 || maxbox > npts/4) maxbox = npts/4;
//    if (maxbox == 0 || maxbox > npts/4) maxbox = npts;
    if (minbox > maxbox) {
	SWAP(minbox, maxbox);
	if (minbox < 2*nfit) minbox = 2*nfit;
    }
//    cout << "minbox " << minbox << endl;
//    cout << "maxbox " << maxbox << endl;
    
    /* Allocate and fill the box size array rs[].  rscale's third argument
       specifies that the ratio between successive box sizes is 2^(1/8). */
    nr = rscale(rs, minbox, maxbox, pow(2.0, 1.0/8.0));
//    cout << "nr " << nr << endl;
    
    MatDoub x(rs[nr]+1, nfit+1, 0.0);	/* matrix of abscissas and their powers, for polyfit(). */
    /* workspace for polyfit() */
    VecDoub betha(nfit+1, 0.0);
    MatDoub covar(nfit+1, nfit+1, 0.0), covar0(nfit+1, nfit+1, 0.0);
    VecInt indxc(nfit+1, 0), indxr(nfit+1, 0), ipiv(nfit+1, 0);

    mse.assign(nr+1, 0);
//    cout << "mse.size()" << mse.size() << endl;
    
    /* Allocate memory for dfa() and the functions it calls. */
    int j;    
//    x.assign(rs[nr]+1, nfit+1, 0.0);
    for (long i = 1; i <= rs[nr]; i++) {
	    x[i][1] = 1.0;
//      printf("x[%ld][1] = %g\n",i, x[i][1]);
	    x[i][2] = i;
//	    printf("x[%ld][2] = %g\n", i, x[i][2]);
	    for (int j = 3; j<=nfit; j++) {
	        x[i][j] = x[i][j-1] * i;
//	        printf("x[%ld][%d] = %g\n", i, j, x[i][j]); 
        }
    }    

    /* Measure the fluctuations of the detrended input data at each box size
       using the DFA algorithm; fill mse[] with these results. */
    dfa(seq, npts, nfit, rs, nr, sw, mse, x, betha, covar, covar0, ipiv, indxc, indxr); 

    VecDoub lrs(nr), lmse(nr);
    /* Output the results. */
    for (int i = 0; i < nr; i++) {
        lrs[i] = log10((double)(rs[i+1]));
        lmse[i] = log10(mse[i+1])/2.0;
//	    printf("%g %g\n", lrs[i], lmse[i]);
	} 
	
	MatDoub M(nr, 2, 1), Mlmse(nr, 1);
    for(int i=0;i<nr;i++)
	{
		M[i][1] = lrs[i];
		Mlmse[i][0] = lmse[i];
	}

	SVD svd(M);
//	cout << svd.u[0][0] << endl;
	MatDoub ut(svd.u.ncols(),svd.u.nrows()), dw(2,2), Minv(2,nr), tmp(svd.v.nrows(),dw.ncols()), a(2,1);
	transp(svd.u,ut);
	diag(svd.w,dw);
	mmultiply(svd.v,dw,tmp);
    mmultiply(tmp,ut,Minv);
   	mmultiply(Minv,Mlmse,a);
    alpha = a[1][0];
//    cout << alpha << endl;
			
/*	Fitab myfit(lrs, lmse);
    alpha = myfit.chi2; */
        
    return alpha;
}

/* rscale() allocates and fills rs[], the array of box sizes used by dfa()
   below.  The box sizes range from (exactly) minbox to (approximately) maxbox,
   and are arranged in a geometric series such that the ratio between
   consecutive box sizes is (approximately) boxratio.  The return value is
   the number of box sizes in rs[].
*/
int rscale(VecDoub &rs, long minbox, long maxbox, double boxratio)
{
    int ir, n;
    long rw;
    int rslen;	/* length of rs[] */

    /* Determine how many scales are needed. */
    rslen = log10(maxbox / (double)minbox) / log10(boxratio) + 1.5;
    //cout << "rslen " << rslen << endl;
    /* Thanks to Peter Domitrovich for pointing out that a previous version
       of the above calculation undercounted the number of scales in some
       situations. */
//    VecDoub auxrs(rslen);
    rs.assign(rslen+1, 0.0);
//    cout << rs.size() << endl;
    for (ir = 1, n = 2, rs[1] = minbox; n <= rslen && rs[n-1] < maxbox; ir++)
        if ((rw = minbox * pow(boxratio, ir) + 0.5) > rs[n-1])  {
            rs[n++] = rw;
//            cout << "rs[" << n << "] = " << rw << endl; 
        }
    if (rs[--n] > maxbox) --n;
//    cout << n << endl;
    return (n);
}



/* Detrended fluctuation analysis
    seq:	input data array
    npts:	number of input points
    nfit:	order of detrending (2: linear, 3: quadratic, etc.)
    rs:		array of box sizes (uniformly distributed on log scale)
    nr:		number of entries in rs[] and mse[]
    sw:		mode (0: non-overlapping windows, 1: sliding window)
   This function returns the mean squared fluctuations in mse[].
*/
void dfa(VecDoub seq, long npts, int nfit, VecDoub &rs, int nr, int sw, VecDoub &mse, MatDoub &x, VecDoub &betha, MatDoub &covar, MatDoub &covar0, VecInt &ipiv, VecInt &indxc, VecInt &indxr)
{

    long i, boxsize, inc, j;
    double stat;
//    cout << "nr = " << nr << endl;
    for (i = 1; i <= nr; i++) {
    
        boxsize = rs[i];
//        cout << "boxsize = " << boxsize << endl;
        if (sw) { inc = 1; stat = (int)(npts - boxsize + 1) * boxsize; }
    	else { inc = boxsize; stat = (int)(npts / boxsize) * boxsize; }
//    	cout << "inc = " << inc << ", stat " << stat << endl;
        for (mse[i] = 0.0, j = 0; j <= npts - boxsize; j += inc) { 
//            cout << "seq[" << j << "] = " << seq[j] << endl;
            mse[i] += polyfit(x, &seq[j], boxsize, nfit, betha, covar, covar0, ipiv, indxc, indxr); 
        }
        mse[i] /= stat;
//        cout << "mse[" << i << "] = " << mse[i] << endl;        
    }
}

/* polyfit() is based on lfit() and gaussj() from Numerical Recipes in C
   (Press, Teukolsky, Vetterling, and Flannery; Cambridge U. Press, 1992).  It
   fits a polynomial of degree (nfit-1) to a set of boxsize points given by
   x[1...boxsize][2] and y[1...boxsize].  The return value is the sum of the
   squared errors (chisq) between the (x,y) pairs and the fitted polynomial.
*/
double polyfit(MatDoub &x, double *y, long boxsize, int nfit, VecDoub &betha, MatDoub &covar, MatDoub &covar0, VecInt &ipiv, VecInt &indxc, VecInt &indxr)
{
    int icol, irow, j, k;
    double big, chisq, pivinv, temp;
    long i;
    static long pboxsize = 0L;

    /* This block sets up the covariance matrix.  Provided that boxsize
       never decreases (which is true in this case), covar0 can be calculated
       incrementally from the previous value. */
    if (pboxsize != boxsize) {	/* this will be false most of the time */
	if (pboxsize > boxsize)	/* this should never happen */
	    pboxsize = 0L;
	if (pboxsize == 0L)	/* this should be true the first time only */
	    for (j = 1; j <= nfit; j++)
		for (k = 1; k <= nfit; k++)
		    covar0[j][k] = 0.0;
	for (i = pboxsize+1; i <= boxsize; i++)
	    for (j = 1; j <= nfit; j++)
		for (k = 1, temp = x[i][j]; k <= j; k++) {
		    //cout << "temp = " << temp << endl;
		    covar0[j][k] += temp * x[i][k];
//		    cout << "covar0[" << j << "][" << k << "] = " << covar0[j][k] << endl; 
		}
	for (j = 2; j <= nfit; j++)
	    for (k = 1; k < j; k++) { 
    		covar0[k][j] = covar0[j][k];
//    		cout << "covar0[" << k << "][" << j << "] = " << covar0[k][j] << endl; 
        }
	pboxsize = boxsize;
    }
    
    for (j = 1; j <= nfit; j++) {
	betha[j] = ipiv[j] = 0;
	for (k = 1; k <= nfit; k++)
	    covar[j][k] = covar0[j][k];
    }
    for (i = 1; i <= boxsize; i++) {
	betha[1] += (temp = y[i]);
	betha[2] += temp * i;
    }
    if (nfit > 2)
	for (i = 1; i <= boxsize; i++)
	    for (j = 3, temp = y[i]; j <= nfit; j++)
		betha[j] += temp * x[i][j];
    for (i = 1; i <= nfit; i++) {
	big = 0.0;
	for (j = 1; j <= nfit; j++)
	    if (ipiv[j] != 1)
		for (k = 1; k <= nfit; k++) {
		    if (ipiv[k] == 0) {
			if ((temp = covar[j][k]) >= big ||
			    (temp = -temp) >= big) {
			    big = temp;
			    irow = j;
			    icol = k;
			}
		    }
		    else if (ipiv[k] > 1) {
    			printf("singular matrix\n");
			    exit(1);
			}
		}
	++(ipiv[icol]);
	if (irow != icol) {
	    for (j = 1; j <= nfit; j++) SWAP(covar[irow][j], covar[icol][j]);
	    SWAP(betha[irow], betha[icol]);
	}
	indxr[i] = irow;
	indxc[i] = icol;
	if (covar[icol][icol] == 0.0) {
	    printf("singular matrix()\n");
	    exit(1);
	}
	pivinv = 1.0 / covar[icol][icol];
	covar[icol][icol] = 1.0;
	for (j = 1; j <= nfit; j++) covar[icol][j] *= pivinv;
	betha[icol] *= pivinv;
	for (j = 1; j <= nfit; j++)
	    if (j != icol) {
		temp = covar[j][icol];
		covar[j][icol] = 0.0;
		for (k = 1; k <= nfit; k++) covar[j][k] -= covar[icol][k]*temp;
		betha[j] -= betha[icol] * temp;
	    }
    }
    chisq = 0.0;
    if (nfit <= 2)
	for (i = 1; i <= boxsize; i++) {
	    temp = betha[1] + betha[2] * i - y[i];
	    chisq += temp * temp;
	}
    else
	for (i = 1; i <= boxsize; i++) {
	    temp = betha[1] + betha[2] * i - y[i];
	    for (j = 3; j <= nfit; j++) temp += betha[j] * x[i][j];
	    chisq += temp * temp;
	}
//	cout << "chisq = " << chisq << endl;
    return (chisq);
}

double media(VecDoub_I &A)
{
    double soma = 0;
    for (int i=0; i<A.size(); i++) 
    {
        soma = soma + A[i];
    }
    soma = soma/A.size();
    return soma;
}
