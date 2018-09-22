#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <dirent.h>
#include <math.h>
#include "nr3.h"
#include "fourier.h"
#include "erf.h"
#include "sort.h"
#include "moment.h"
#include "gamma.h"
#include "incgammabeta.h"
#include "stattests.h"
#include "spectrum.h"
#include "gaussj.h"
#include "svd.h"
#include <vector>
#include "dfa.h"
#include "auxiliares.h"
#include "powspecscale.h"
#include "Analys.h"
#include "buildGNL.h"

#define kmax 4 // Maximum variational degree
#define nmi 1 // Maximum quantity of analytical measurements

int main(int argc, char** argv)
{

   time_t tstart, tend; // time counter
   tstart = time(0);
   
   string arg1, xmlfname;
   int nfiles = 0; // nfiles counts the number of '.gnl' files in input folder
   vector<string> files;
   VecDoub elements;
    
   if (argc > 1) {
      arg1 = argv[1];
      xmlfname = arg1 +".xml";
   } else {
      cout << "Usage is: gnlsys <input_gnl_folder_name>" << endl;
      exit(0);
   }

   if (isdir(arg1)) { // keeps only the '.gnl' files
      files = getdir(arg1);
      for (unsigned int i = 0; i < (int)files.size(); i++) {
         string s = getFileExtension(files[i]);
         if ( s == ".gnl") {
            nfiles++;
         }
         else {				
            files.erase(files.begin() + i); // ignores the 'non-gnl' files
         } //else
      } //inner for 
   } //outer for

   MatInt lb(kmax,nfiles);
   VecInt kp(nfiles), n(nfiles);
   vector<string> filename(nfiles);
   MatDoub mu; // mu stores the analytical measures for all GNLs
   mu.assign(nfiles,nmi,0);
	
   for (unsigned int i = 0; i < nfiles; i++) {
      kp[i] = 0;
      n[i] = 1;
      filename[i] = arg1 +"/" + files[i];
      ifstream myfile(filename[i].c_str());

      if (!myfile) {
         cout << "Could not open file" << endl;
	       exit(1);
      }
      
      /* Reading the file header */
      for (int j=0;j<kmax;j++) {
         myfile >> lb[j][i];
         if (lb[j][i]==0) { 
            lb[j][i]++;
         }
         else {
            kp[i]++;
         }
         n[i] = n[i] * lb[j][i];
      } //for j
        	
      elements.resize(n[i]); // resizes the elements vector to the total quantity of measurements in each GNL
        	        	
      for (int j=0;j<n[i];j++) {
         myfile >> elements[j]; // reading the values
      }
		
      myfile.close();      	
      VecInt l(kmax);
  
      for (int j=0; j<kmax; j++) {
         l[j] = lb[j][i];
      }
            
      VecDoub auxMu;
      auxMu.assign(nmi, 0);
      Analys(kp[i], l, elements, auxMu); // calling the function to analyse the data 
      for (int j=0; j<nmi; j++) {
         mu[i][j] = auxMu[j];
      }         	            	     	
   } //for i
	
   buildGNL(filename, kp, lb, mu, xmlfname); // calling the function to build the GNLs and save metadata into the XML file

   tend = time(0);
   cout << "It took " << difftime(tend, tstart) << " second(s)." << endl; // It shows the time spent in data analysis

   return 0;

}

