//It returns true if 'dir' is a directory and false otherwise
bool isdir(string dir) {

   DIR *dp;
   struct dirent *dirp;
   if ((dp = opendir(dir.c_str())) == NULL)	return false;
   else	return true;

}

//It returns 'file' extension
string getFileExtension(string file){

   string ext = "";
   for(int i=0; i<file.length(); i++){
      if(file[i] == '.'){
         for(int j = i; j<file.length(); j++){
	    ext += file[j];
	 }
	 return ext;
      }
   }
   return ext;

}

//It returns the names of the elements in dir
vector<string> getdir(string dir) {
   
   vector<string> files;
   DIR *dp;
   struct dirent *dirp;
   dp = opendir(dir.c_str());
   while ((dirp = readdir(dp)) != NULL) {
      string s = string(dirp->d_name);
      if (s != "." & s != "..")   files.push_back(s);
   }
   closedir(dp);
   return files;

}

//It shifts vector x in lag units to right and stores the result in sx
void shift(VecDoub_I &x, VecDoub_O &sx, int &lag) {

   int j=0, n=x.size();
   for(int i=0;i<(n-lag);i++) {
      sx[i] = x[i+lag];
   }
   for(int i=(n-lag);i<n;i++) {
      sx[i] = x[j++];
   }

}

//It calculates log10 of a vector x and stores the result in lx
void vlog10(VecDoub_I &x, VecDoub_O &lx) {

   int n = x.size();
   for(int i=0;i<n;i++)
   {
      lx[i] = log10(x[i]);
   }

}

//It calculates the transpose of matrix a(mxn) and stores it in at(nxm)
void transp(MatDoub_I &a, MatDoub_O &at) {

   for(int i=0;i<a.nrows();i++) {
      for(int j=0;j<a.ncols();j++) {
  	 at[j][i] = a[i][j];
      }
   }

}

//It builds a diagonal matrix dw_mxn with the elements of vector w(m)
void diag(VecDoub_I &w, MatDoub_O &dw) {

   for(int i=0;i<dw.nrows();i++) {
      for(int j=0;j<dw.ncols();j++) {
         if (i==j)  dw[i][j] = 1/w[i]; 
	 else dw[i][j] = 0; 
      }
   }

}

//It multiplies a(mxn) by b(nxp) and stores the result in c(mxp)
void mmultiply(MatDoub_I &a, MatDoub_I &b, MatDoub_O &c) {

   for(int k=0;k<a.nrows();k++) {
      for(int m=0;m<b.ncols();m++) {
         Doub soma = 0;
         for(int j=0;j<a.ncols();j++) {
            Doub tmp = a[k][j]*b[j][m];
   	    soma += tmp;
         }
         c[k][m] = soma;
      }
   }

}

//It calculates the absolute value of a complex number
Doub absZ(Doub &re, Doub &im) {
   
   Doub z = sqrt(pow(re,2)+pow(im,2));
   return z;

}
