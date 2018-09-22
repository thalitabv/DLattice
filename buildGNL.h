using namespace std;

/* It writes the metadata and analytical measurements into the XML file */

void buildGNL(vector<string> files, VecInt_I kp, MatInt_I lb, MatDoub_I mi, string xmlfname) 
{
 
   int nfiles = kp.size();
   fstream xmlfile;
   xmlfile.open(xmlfname.c_str(), fstream::in | fstream::out | fstream::app);
   xmlfile.seekg(0, ios::end); // put the "cursor" at the end of the file
   int length = xmlfile.tellg(); // find the position of the cursor
      
   if ( length == 0 ) { 
      xmlfile <<"<?xml version=\"1.0\"?>\n"
              << "\n<GNLSystem>\n";
      for (int i=0; i<nfiles; i++) {
         xmlfile << "<GNL " // root tag
                 << "filename = \"" << files[i] << "\" "
    	         << " kappa = \"" << kp[i]+1 << "\" "
    		 << " lambda0 = \"" << lb[0][i] << "\" "
    		 << " lambda1 = \"" << lb[1][i] << "\" "
    		 << " lambda2 = \"" << lb[2][i] << "\" "
    		 << " lambda3 = \"" << lb[3][i] << "\" ";
         for (int j=0; j<mi.ncols(); j++) {	
            xmlfile << " mi" << j << " =\" " << mi[i][j] << "\" "; 
         }	
    	 xmlfile << ">\n"
    	         << "</GNL>\n";
      }

      xmlfile << "</GNLSystem>";
      xmlfile.close(); 
	
   }

}
