////Numerical.cpp

///Some numerical routines
#include "stdafx.h"
#include "Numerical.h"

extern Vec_INT *ija_p;
extern Vec_DP *sa_p;

//#define NRANSI
//#include "nrutil.h"
//#define MAXSTP 10000
//#define TINY 1.0e-30
//
//extern int kmax,kount;
//extern double *xp,**yp,dxsav;
//
//void odeint(double vec_a, double vec_b, double vec_c, double vec_d, double vec_e, double vec_f, double ystart[], int nvar, double x1, double x2, double eps, double h1,
//	double hmin, int *nok, int *nbad,
//	void (*derivs)(double, double, double, double, double, double, double, double [], double []),
//	void (*rkqs)(double, double, double, double, double, double, double [], double [], int, double *, double, double, double [],
//	double *, double *, void (*)(double, double, double, double, double, double, double, double [], double [])))
//{
//	int nstp,i;
//	double xsav,x,hnext,hdid,h;
//	double *yscal,*y,*dydx;
//
//	yscal=vector(1,nvar);
//	y=vector(1,nvar);
//	dydx=vector(1,nvar);
//	x=x1;
//	h=SIGN(h1,x2-x1);
//	*nok = (*nbad) = kount = 0;
//	for (i=1;i<=nvar;i++) y[i]=ystart[i];
//	if (kmax > 0) xsav=x-dxsav*2.0;
//	for (nstp=1;nstp<=MAXSTP;nstp++) {
//		(*derivs)(x,vec_a, vec_b, vec_c, vec_d, vec_e, vec_f, y,dydx);
//		for (i=1;i<=nvar;i++)
//			yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
//		if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
//			xp[++kount]=x;
//			for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
//			xsav=x;
//		}
//		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
//		(*rkqs)(vec_a, vec_b, vec_c, vec_d, vec_e, vec_f, y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);
//		if (hdid == h) ++(*nok); else ++(*nbad);
//		if ((x-x2)*(x2-x1) >= 0.0) {
//			for (i=1;i<=nvar;i++) ystart[i]=y[i];
//			if (kmax) {
//				xp[++kount]=x;
//				for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
//			}
//			free_vector(dydx,1,nvar);
//			free_vector(y,1,nvar);
//			free_vector(yscal,1,nvar);
//			return;
//		}
//		if (fabs(hnext) <= hmin) 
//			nrerror("Step size too small in odeint");
//		h=hnext;
//	}
//	nrerror("Too many steps in routine odeint");
//}
//
//void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
//	double hmin, int *nok, int *nbad,
//	void (*derivs)(double, double [], double []),
//	void (*rkqs)(double [], double [], int, double *, double, double, double [],
//	double *, double *, void (*)(double, double [], double [])))
//{
//	int nstp,i;
//	double xsav,x,hnext,hdid,h;
//	double *yscal,*y,*dydx;
//
//	yscal=vector(1,nvar);
//	y=vector(1,nvar);
//	dydx=vector(1,nvar);
//	x=x1;
//	h=SIGN(h1,x2-x1);
//	*nok = (*nbad) = kount = 0;
//	for (i=1;i<=nvar;i++) y[i]=ystart[i];
//	if (kmax > 0) xsav=x-dxsav*2.0;
//	for (nstp=1;nstp<=MAXSTP;nstp++) {
//		(*derivs)(x,y,dydx);
//		for (i=1;i<=nvar;i++)
//			yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
//		if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
//			xp[++kount]=x;
//			for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
//			xsav=x;
//		}
//		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
//		(*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);
//		if (hdid == h) ++(*nok); else ++(*nbad);
//		if ((x-x2)*(x2-x1) >= 0.0) {
//			for (i=1;i<=nvar;i++) ystart[i]=y[i];
//			if (kmax) {
//				xp[++kount]=x;
//				for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
//			}
//			free_vector(dydx,1,nvar);
//			free_vector(y,1,nvar);
//			free_vector(yscal,1,nvar);
//			return;
//		}
//		if (fabs(hnext) <= hmin) 
//			nrerror("Step size too small in odeint");
//		h=hnext;
//	}
//	nrerror("Too many steps in routine odeint");
//}



////Routines to solve Runge-Kutta Integration
void rkck(Vec_I_DP &y, Vec_I_DP &dydx, const DP x,
	const DP h, Vec_O_DP &yout, Vec_O_DP &yerr,
	void derivs(const DP, Vec_I_DP &, Vec_O_DP &))
{
	static const DP a2=0.2, a3=0.3, a4=0.6, a5=1.0, a6=0.875,
		b21=0.2, b31=3.0/40.0, b32=9.0/40.0, b41=0.3, b42 = -0.9,
		b43=1.2, b51 = -11.0/54.0, b52=2.5, b53 = -70.0/27.0,
		b54=35.0/27.0, b61=1631.0/55296.0, b62=175.0/512.0,
		b63=575.0/13824.0, b64=44275.0/110592.0, b65=253.0/4096.0,
		c1=37.0/378.0, c3=250.0/621.0, c4=125.0/594.0, c6=512.0/1771.0,
		dc1=c1-2825.0/27648.0, dc3=c3-18575.0/48384.0,
		dc4=c4-13525.0/55296.0, dc5 = -277.00/14336.0, dc6=c6-0.25;
	int i;

	int n=y.size();
	Vec_DP ak2(n),ak3(n),ak4(n),ak5(n),ak6(n),ytemp(n);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+b21*h*dydx[i];
	derivs(x+a2*h,ytemp,ak2);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
	derivs(x+a3*h,ytemp,ak3);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
	derivs(x+a4*h,ytemp,ak4);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
	derivs(x+a5*h,ytemp,ak5);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
	derivs(x+a6*h,ytemp,ak6);
	for (i=0;i<n;i++)
		yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
	for (i=0;i<n;i++)
		yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
}

void rkqs(Vec_IO_DP &y, Vec_IO_DP &dydx, DP &x, const DP htry,
	const DP eps, Vec_I_DP &yscal, DP &hdid, DP &hnext,
	void derivs(const DP, Vec_I_DP &, Vec_O_DP &))
{
	const DP SAFETY=0.9, PGROW=-0.2, PSHRNK=-0.25, ERRCON=1.89e-4;
	int i;
	DP errmax,h,htemp,xnew;

	int n=y.size();
	h=htry;
	Vec_DP yerr(n),ytemp(n);
	for (;;) {
		rkck(y,dydx,x,h,ytemp,yerr,derivs);
		errmax=0.0;
		for (i=0;i<n;i++) errmax=MAX(errmax,fabs(yerr[i]/yscal[i]));
		errmax /= eps;
		if (errmax <= 1.0) break;
		htemp=SAFETY*h*pow(errmax,PSHRNK);
		h=(h >= 0.0 ? MAX(htemp,0.1*h) : MIN(htemp,0.1*h));
		xnew=x+h;
		if (xnew == x) MessageBox(NULL, "stepsize underflow in rkqs", "error", MB_OK);
	}
	if (errmax > ERRCON) hnext=SAFETY*h*pow(errmax,PGROW);
	else hnext=5.0*h;
	x += (hdid=h);
	for (i=0;i<n;i++) y[i]=ytemp[i];
}




/*--------------------------------------------------------------*/
/////functions for solving sparse linear system
void sprsin(Mat_I_DP &a, const DP thresh, Vec_O_DP &sa, Vec_O_INT &ija)
{
	int i,j,k;

	int n=a.nrows();
	int nmax=sa.size();
	for (j=0;j<n;j++) sa[j]=a[j][j];
	ija[0]=n+1;
	k=n;
	for (i=0;i<n;i++) {
		for (j=0;j<n;j++) {
			if (fabs(a[i][j]) >= thresh && i != j) {
				if (++k > nmax) MessageBox(NULL,"sprsin: sa and ija too small", "error", MB_OK);
				sa[k]=a[i][j];
				ija[k]=j;
			}
		}
		ija[i+1]=k+1;
	}
}


void linbcg(Vec_I_DP &b, Vec_IO_DP &x, const int itol, const DP tol,
	const int itmax, int &iter, DP &err)
{
	DP ak,akden,bk,bkden=1.0,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm;
	const DP EPS=1.0e-14;
	int j;

	int n=b.size();
	Vec_DP p(n),pp(n),r(n),rr(n),z(n),zz(n);
	iter=0;
	atimes(x,r,0);
	for (j=0;j<n;j++) {
		r[j]=b[j]-r[j];
		rr[j]=r[j];
	}
	//atimes(r,rr,0);
	if (itol == 1) {
		bnrm=snrm(b,itol);
		asolve(r,z,0);
	}
	else if (itol == 2) {
		asolve(b,z,0);
		bnrm=snrm(z,itol);
		asolve(r,z,0);
	}
	else if (itol == 3 || itol == 4) {
		asolve(b,z,0);
		bnrm=snrm(z,itol);
		asolve(r,z,0);
		znrm=snrm(z,itol);
	} else MessageBox(NULL, "illegal itol in linbcg", "error", MB_OK);
	//cout << fixed << setprecision(6);
	while (iter < itmax) {
		++iter;
		asolve(rr,zz,1);
		for (bknum=0.0,j=0;j<n;j++) bknum += z[j]*rr[j];
		if (iter == 1) {
			for (j=0;j<n;j++) {
				p[j]=z[j];
				pp[j]=zz[j];
			}
		} else {
			bk=bknum/bkden;
			for (j=0;j<n;j++) {
				p[j]=bk*p[j]+z[j];
				pp[j]=bk*pp[j]+zz[j];
			}
		}
		bkden=bknum;
		atimes(p,z,0);
		for (akden=0.0,j=0;j<n;j++) akden += z[j]*pp[j];
		if(akden != 0)
		    ak=bknum/akden;
		else
			ak = bknum;
		atimes(pp,zz,1);
		for (j=0;j<n;j++) {
			x[j] += ak*p[j];
			r[j] -= ak*z[j];
			rr[j] -= ak*zz[j];
		}
		asolve(r,z,0);
		if (itol == 1)
		{
			if(bnrm != 0)
				err=snrm(r,itol)/bnrm;
			else
				err=snrm(r,itol);
		}
		else if (itol == 2)
			err=snrm(z,itol)/bnrm;
		else if (itol == 3 || itol == 4) {
			zm1nrm=znrm;
			znrm=snrm(z,itol);
			if (fabs(zm1nrm-znrm) > EPS*znrm) {
				dxnrm=fabs(ak)*snrm(p,itol);
				err=znrm/fabs(zm1nrm-znrm)*dxnrm;
			} else {
				err=znrm/bnrm;
				continue;
			}
			xnrm=snrm(x,itol);
			if (err <= 0.5*xnrm) err /= xnrm;
			else {
				err=znrm/bnrm;
				continue;
			}
		}
		cout << "iter=" << setw(4) << iter+1 << setw(12) << err << endl;
		if (err <= tol) break;
	}
}

DP snrm(Vec_I_DP &sx, const int itol)
{
	int i,isamax;
	DP ans;

	int n=sx.size();
	if (itol <= 3) {
		ans = 0.0;
		for (i=0;i<n;i++) ans += sx[i]*sx[i];
		return sqrt(ans);
	} else {
		isamax=0;
		for (i=0;i<n;i++) {
			if (fabs(sx[i]) > fabs(sx[isamax])) isamax=i;
		}
		return fabs(sx[isamax]);
	}
}


void atimes(Vec_I_DP &x, Vec_O_DP &r, const int itrnsp)
{
	if (itrnsp) sprstx(*sa_p,*ija_p,x,r);
	else sprsax(*sa_p,*ija_p,x,r);
}


void asolve(Vec_I_DP &b, Vec_O_DP &x, const int itrnsp)
{
	int i;

	int n=b.size();
	for(i=0;i<n;i++) x[i]=((*sa_p)[i] != 0.0 ? b[i]/(*sa_p)[i] : b[i]);
}


void sprsax(Vec_I_DP &sa, Vec_I_INT &ija, Vec_I_DP &x, Vec_O_DP &b)
{
	int i,k;

	int n=x.size();
	if (ija[0] != n+1)
		MessageBox(NULL, "sprsax: mismatched vector and matrix", "error", MB_OK);
	for (i=0;i<n;i++) {
		b[i]=sa[i]*x[i];
		for (k=ija[i];k<ija[i+1];k++) {
			b[i] += sa[k]*x[ija[k]];
		}
	}
}


void sprstx(Vec_I_DP &sa, Vec_I_INT &ija, Vec_I_DP &x, Vec_O_DP &b)
{
	int i,j,k;

	int n=x.size();
	if (ija[0] != (n+1))
		MessageBox(NULL, "mismatched vector and matrix in sprstx", "error", MB_OK);
	for (i=0;i<n;i++) b[i]=sa[i]*x[i];
	for (i=0;i<n;i++) {
		for (k=ija[i];k<ija[i+1];k++) {
			j=ija[k];
			b[j] += sa[k]*x[i];
		}
	}
}


