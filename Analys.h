using namespace std;

/*The input '.gnl' file is formatted as following:
- 1st line: extension coefficients in ascending order: lambda_3 lambda_2 lambda_1 lambda_0
Ex: - time series composed of 20 points -> 1st line: 0 0 0 20
    - time series composed of 20 arrays of size 64x64 -> 1st line 0 64 64 20
    - time series composed of 20 hypercubes of size 64x64x64 -> 1st line 64 64 64 20
- from 2nd line: time series values
Reference: Veronese, T. B. "Generalized Numerical Lattice: a new concept for analytical
visualization and representation of time series systems", PhD Thesis, INPE, 2011. */

//This function returns the average value in a vector
void med(VecDoub_I &x, Doub &md)
{
   int n = x.size();
   md = 0;
   for (int i=0; i<n; i++)
   {
      md += x[i];
   }
   md = md/n;
}

void Analys(int k, VecInt &lb, VecDoub &elements, VecDoub &mi) {
       
   int nMi = 0, cont = 0; // analytical measurements counter
   k++; // variational degree
                     
   if (k==2) { // A(t)
            
      Int n = lb[3]; // n is the extension coefficient for the fundamental domain (most usual: time) 
      Doub alpha; // alpha is the scale exponent to be calculated by DFA technique
      VecDoub val = elements; // val is the vector to keep the time series measurements 	

      /* Call here the functions to analyse A(t) */

      alpha = dfa_call(val); // we calculate the scale exponent via DFA and store the result in alpha

      mi[0] = roundf(alpha * 1000) / 1000; // alpha is the first analytical measure
	    /* we truncate the value of mi[i] to better display the resulting GNG */
      /* Obs.: additional analytical measures should be stored in mi[0],...,mi[nmi]; nmi is declared as constant in gnlsys.cpp */
           
   } // k==2
            
   else if (k==3) { // A(t,x)  		

      cont = 0;
      MatDoub val(lb[3],lb[2]);
    	for (int i=0;i<lb[3];i++) {
    	   for (int j=0;j<lb[2];j++) {
    	      val[i][j] = elements[cont];    	           
		        cont++;
    	   }
      }	

   } // k==3
	    
   else if (k==4) { // A(t,x,y)
	    
      cont = 0;
      Doub issquare = (lb[2]!=lb[1]); // some methods can be applied only to square matrices
	    	  
      vector<MatDoub> val(lb[3]); // the 3D matrix is declared
      MatDoub matrix(lb[2],lb[1]);

   	  for (int i=0;i<lb[3];i++) { // i counts the time variation
		     val[i].resize(lb[2],lb[1]);
    	   for (int j=0;j<lb[2];j++) {
	    	    for (int k=0;k<lb[1];k++) {
	    	       val[i][j][k] = elements[cont];
               matrix[j][k] = val[i][j][k]; // matrix stores in a 2D array the values in (x,y) for time i
			         cont++;
    		    }
	       }
		  } 
      
      /* Call here the functions to analyse A(i,x,y) */

   }	// k==4
    
   /* Just an example for expanding kappa; must update kmax at gnlsys.cpp */
/*   else if (k==5) { // A(t,x,y,z)

      cont = 0;
    	NRmatrix<MatDoub> val(lb[0],lb[1]); // the 3D matrix is declared

      for (int i=0;i<lb[0];i++) { // i counts the time variation
         for (int j=0;j<lb[1];j++) {
			      val[i][j].resize(lb[2],lb[3]);
    			  for (int k=0;k<lb[2];k++) {
    	         for (int q=0;q<lb[3];q++) {
 				          val[i][j][k][q] = elements[cont];
    				   }
    			  }
         }
      }
    
   } // k==5
*/
	    
} // Analys
