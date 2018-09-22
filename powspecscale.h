//Calcula a escala do espectro de potencia de uma serie x
Doub powspecscale(VecDoub_I &x)
{
	int cont, cont2;
	Doub temp, Wss, lf1, lf2, lfk, dlfb, lfrange, Nlfb, p=0, n=x.size(), m=n/2, soma=0;
	VecDoub xx(n),w(n),wx(n),awx(n),s0(n),S(m),logS,f(m),logf,lfb,Slfb;

	//Ensure data has a mean of zero
	Doub ave, adev, sdev, var, skew, curt;
	moment(x, ave, adev, sdev, var, skew, curt);
	for(int i=0;i<n;i++) 
	{
		xx[i] = x[i]-ave;
		//cout << xx[i] << endl;
	}
	
	//Welch window 
	for(int i=0;i<n;i++) {
		w[i] = 1-pow(((i+1-n/2)/(n/2)),2);
		wx[i] = pow(w[i],2);
		//cout << w[i] << endl;
	}
	moment(wx, Wss, adev, sdev, var, skew, curt);
	//cout << Wss << endl;

	//calculate the spectrum 
	for(int i=0;i<n;i++)
	{
		wx[i] = w[i]*xx[i];
		//cout << wx[i] << endl;
	}
	realft(wx,1);
	/*for(int i=0;i<n;i++)
	{
		cout << wx[i] << endl;
	}*/
	temp = abs(wx[1]);
	cont = 1, cont2 = n-1;
	awx[0] = abs(wx[0]);
	for(int i=1;i<n-1;i+=2)
	{
		awx[cont] = absZ(wx[i+1],wx[i+2]);
		awx[cont2] = absZ(wx[i+1],wx[i+2]);
		cont++;
		cont2--;
	}
	awx[m] = temp;
	//cout << m << "--" << wx[m] << endl;

	for(int i=0;i<n;i++)
	{
		//cout << awx[i] << endl;
		s0[i] = (1/Wss)*(2/n)*pow(awx[i],2);
		//cout << s0[i] << endl;
	}

	//cout << "S - f" << endl;
	for(int i=0;i<m;i++)
	{
		S[i] = s0[i];
		f[i] = (i+1)/n;
		//cout << S[i] << " - " << f[i] << endl;
	}
	
	logS.resize(m);
	//cout << S.size() << endl;
	logf.resize(m);
	//cout << f.size() << endl;
	vlog10(S,logS);
	vlog10(f,logf);
	/*cout << "log(S) - log(f)" << endl;
	for(int i=0;i<m;i++)
	{
		cout << logS[i] << " - " << logf[i] << endl;
	}*/

	// estimate bin size as distance between first two points on log frequency scale 
	lf1 = log10(f[0]);
	//cout << lf1 << endl;
	lf2 = log10(f[1]);
	//cout << lf2 << endl;
	lfk = log10(f[m-1]);
	//cout << lfk << endl;
	dlfb = lf2 - lf1;
	//cout << dlfb << endl;
	
	// bin the data
	lfrange = lfk - lf1;
	//cout << lfrange << endl;
	Nlfb = ceil(lfrange/dlfb);
	//cout << Nlfb << endl;
	temp=lf1-0.5*dlfb;
	//cout << temp << endl;
	lfb.resize(Nlfb);
	Slfb.resize(Nlfb);
	for(int i=0;i<Nlfb;i++)	{
		lfb[i] = temp + (i+1)*dlfb; 
		//cout << lfb[i] << endl;
	}
	
	VecInt ind(f.size()),nind(Nlfb);
	cont=0;
	for(int i=0;i<f.size();i++) 
	{
		if ((lfb[0]-0.5*dlfb <= logf[i]) & (logf[i] <= lfb[0]+0.5*dlfb))
		{
			ind[cont] = i;
			//cout << ind[cont] << endl;
			cont++;
		}
	}
	nind[0] = cont;
	//cout << nind[0] << endl;
	for(int i=0;i<cont;i++)
	{
		soma = soma + logS[ind[i]];
	}
	Slfb[0] = soma/cont;
	//cout << Slfb[0] << endl;
	for(int i=1;i<Nlfb-1;i++)
	{
		cont = 0;
		//cout << "i = " << i << endl;
		for(int j=0;j<logf.size();j++)
		{
			if ((lfb[i]-0.5*dlfb < logf[j]) & (logf[j] <= lfb[i]+0.5*dlfb))
			{
				ind[cont] = j;
				//cout << ind[cont] << endl;
				cont++;
			}
		}
		soma = 0;
		for(int j=0;j<cont;j++)
		{
			soma = soma + logS[ind[j]];
		}
		Slfb[i] = soma/cont;
		//cout << Slfb[i] << endl;
		nind[i] = cont;
		//cout << nind[i] << endl;
	}

	cont = 0;
	for(int i=0;i<f.size();i++) 
	{
		if (lfb[Nlfb-1]-0.5*dlfb <= logf[i])
		{
			ind[cont] = i;
			//cout << ind[cont] << endl;
			cont++;
		}
	}

	nind[Nlfb-1] = cont;
	//cout << nind[Nlfb-1] << endl;
	soma = 0;
	for(int i=0;i<cont;i++)
	{
		soma = soma + logS[ind[i]];
	}
	Slfb[Nlfb-1] = soma/cont;

	//fit a line to the log-log plot
	MatDoub A(Nlfb,2,1);
	//cout << A.nrows() << "x" << A.ncols() << endl;
	for(int i=0;i<Nlfb;i++)
	{
		A[i][1] = lfb[i];
		//cout << lfb[i] << endl;
		//cout << Slfb[i] << endl;
	}
	/*for(int i=0;i<A.nrows();i++) {
		for(int j=0;j<A.ncols();j++) {
			cout << A[i][j] << " ";
		}
		cout << endl;
	}*/
	
	SVD svd(A);
	//cout << svd.u[0][0] << endl;
	MatDoub ut(svd.u.ncols(),svd.u.nrows()), dw(2,2), Ainv(2,Nlfb), tmp(svd.v.nrows(),dw.ncols()), mslfb(Nlfb,1), a(2,1);
	for(int i=0;i<Nlfb;i++) {
		mslfb[i][0] = Slfb[i];
	}
	transp(svd.u,ut);
	diag(svd.w,dw);
	/*cout << "U = " << endl;
	for(int i=0;i<svd.u.nrows();i++) {
		for(int j=0;j<svd.u.ncols();j++) {
			cout << svd.u[i][j] << " ";
		}
		cout << endl;
	}
	cout << "Ut = " << endl;
	for(int i=0;i<ut.nrows();i++) {
		for(int j=0;j<ut.ncols();j++) {
			cout << ut[i][j] << " ";
		}
		cout << endl;
	}
	cout << "V" << endl;
	for(int i=0;i<svd.v.nrows();i++) {
		for(int j=0;j<svd.v.ncols();j++) {
			cout << svd.v[i][j] << " ";
		}
		cout << endl;
	}
	cout << "*DW" << endl;
	for(int i=0;i<dw.nrows();i++) {
		for(int j=0;j<dw.ncols();j++) {
			cout << dw[i][j] << " ";
		}
		cout << endl;
	}
	cout << "=" << endl;*/
	mmultiply(svd.v,dw,tmp);
	/*for(int i=0;i<tmp.nrows();i++) {
		for(int j=0;j<tmp.ncols();j++) {
			cout << tmp[i][j] << " ";
		}
		cout << endl;
	}*/
	mmultiply(tmp,ut,Ainv);
	/*for(int i=0;i<Ainv.nrows();i++) {
		for(int j=0;j<Ainv.ncols();j++) {
			cout << Ainv[i][j] << " ";
		}
		cout << endl;
	}*/
	mmultiply(Ainv,mslfb,a);
	p = -1.*a[1][0];

	return p;
}
