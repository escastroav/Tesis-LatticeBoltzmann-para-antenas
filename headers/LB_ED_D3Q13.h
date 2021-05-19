/*
Reproduccion LBM para ED

Actualizacion:
 - Funciona para 3D, se requieren las librerías Vector.h , Distribution.h.
 - La funcion de equilibrio debe estar escrita en terminos de epsilon y mu relativos.
 - Las cantidades macroscopicas son instancias de la clase definida en Vector.h .
 - Para Efecto skin se recomienda sigma0 = 0.1, omega = M_PI*0.02. Imprimir en t=612, t=630, t=655.
 */
#include <iostream>
#include <fstream>
#include <cmath>
//#include "Distribution.h"
//#include "Vector.h"
using namespace std;

int Lx=100;
int Ly=100;
int Lz=100;

const int rS=2;
const int jS=2;
const int pS=3;
const int iS=4;

int sizef = Lx*Ly*Lz*rS*jS*pS*iS;
int sizef0 = Lx*Ly*Lz*rS;


const double C=1.0/sqrt(2.0);
const double eps0=1.0;
const double mu0=2.0;

const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;
	       
class LatticeBoltzmann{
protected:
  double V0[3]={0.0,0.0,0.0};
  double e0[3]={0.0,0.0,0.0};
  double b0[3]={0.0,0.0,0.0};
  
  double V[3][4][3];
  double vx[12];	double vy[12];	double vz[12];
  double ex[24];	double ey[24];	double ez[24];
  double bx[24];	double by[24];	double bz[24];
 
  double * f = nullptr; double * fnew = nullptr;
  double * f0 = nullptr; double * f0new = nullptr;

  //Distribution f= Distribution(Lx,Ly,Lz,false);		Distribution fnew = Distribution(Lx,Ly,Lz,false);
  //Distribution f0 = Distribution(Lx,Ly,Lz,true);	Distribution f0new = Distribution(Lx,Ly,Lz,true);
  
public:
  LatticeBoltzmann(void);
  ~LatticeBoltzmann(void);
  void ResizeDomain(int Lx0, int Ly0, int Lz0);
  double rho(int ix, int iy, int iz, bool useNew);
  //campos
  double Dx(int ix, int iy, int iz, bool useNew);
  double Dy(int ix, int iy, int iz, bool useNew);
  double Dz(int ix, int iy, int iz, bool useNew);
  double Bx(int ix, int iy, int iz, bool useNew);
  double By(int ix, int iy, int iz, bool useNew);
  double Bz(int ix, int iy, int iz, bool useNew);
  double Ex(double&Dx0, double&epsr);
  double Ey(double&Dy0, double&epsr);
  double Ez(double&Dz0, double&epsr);
  double Hx(double&Bx0, double&mur);
  double Hy(double&By0, double&mur);
  double Hz(double&Bz0, double&mur);
  double Jx(double&Ex0, double&sigma);
  double Jy(double&Ey0, double&sigma);
  double Jz(double&Ez0, double&sigma);
  //campos auxiliares
  double Jxp(double&Ex0, double&denominator);
  double Jyp(double&Ey0, double&denominator);
  double Jzp(double&Ez0, double&denominator);
  double Exp(double&Ex0, double&Jxp0, double&factor);
  double Eyp(double&Ey0, double&Jyp0, double&factor);
  double Ezp(double&Ez0, double&Jzp0, double&factor);
  //constantes dielectricas relativas
  double epsr(int ix, int iy, int iz){return 1.0;}; 
  double mur(int ix, int iy, int iz){return 1.0;};
  double sigma(int ix, int iy, int iz){return 0.0;};
  //funciones de equilibrio
  double feq(double&epsr,double&mur,double&vJp,double&eEp,double&bB,int r);
  double feq0(double&rho0);
  void Colisione(int &t);
  void Adveccione(void);
  void Inicie(void);
  //void ImponerCampos(double&rho0,double*&D0,double*&B0,double*&H0,double*&E0,double*&J0,double*&Jp0,double*&Ep0,int t);
  void Imprimase(const char* fileName,bool useNew);
};
LatticeBoltzmann::LatticeBoltzmann(void){
  int imp=0,rec=0,rec1=0,rec2=0;
  double pi4=M_PI*0.25, sq2=sqrt(2.0);
  f = new double[sizef];	fnew = new double[sizef];
  f0 = new double[sizef0];	f0new = new double[sizef0];
   
  V[0][0][0]=1;	V[1][0][0]=1;	V[2][0][0]=0;    //Plano XY p=0
  V[0][1][0]=-1;V[1][1][0]=1;	V[2][1][0]=0;	 //Plano XY p=0
  V[0][2][0]=-1;V[1][2][0]=-1;	V[2][2][0]=0;	 //Plano XY p=0
  V[0][3][0]=1;	V[1][3][0]=-1;	V[2][3][0]=0;	 //Plano XY p=0

  V[0][0][1]=1; V[1][0][1]=0; 	V[2][0][1]=1;	//Plano XZ p=1
  V[0][1][1]=-1;V[1][1][1]=0; 	V[2][1][1]=1;	//Plano XZ p=1
  V[0][2][1]=-1;V[1][2][1]=0; 	V[2][2][1]=-1;	//Plano XZ p=1
  V[0][3][1]=1; V[1][3][1]=0; 	V[2][3][1]=-1;	//Plano XZ p=1

  V[0][0][2]=0;	V[1][0][2]=1;	V[2][0][2]=1;	//Plano YZ p=2
  V[0][1][2]=0;	V[1][1][2]=-1;	V[2][1][2]=1;	//Plano YZ p=2
  V[0][2][2]=0;	V[1][2][2]=-1;	V[2][2][2]=-1;	//Plano YZ p=2
  V[0][3][2]=0;	V[1][3][2]=1;	V[2][3][2]=-1;	//Plano YZ p=2
  
  //Los vectores Velocidad V[p][i]=V^p_i como vectores
  //v0.cargue(V0[0],V0[1],V0[2]);
  for(int p=0;p<3;p++)
    for(int i=0;i<4;i++){
      //v[p*4+i].cargue(V[0][i][p],V[1][i][p],V[2][i][p]);
      vx[p*4+i] = V[0][i][p];
      vy[p*4+i] = V[1][i][p];
      vz[p*4+i] = V[2][i][p];
  }
  
  for(int p=0;p<3;p++)
    for(int i=0;i<4;i++)
      {
	rec=i%4; rec1=(i+1)%4; rec2=(i+2)%4;

	ex[rec1*2+p*8]=vx[p*4+rec2]*0.5;
	ex[1+rec1*2+p*8]=vx[p*4+rec]*0.5;
	ey[rec1*2+p*8]=vy[p*4+rec2]*0.5;
	ey[1+rec1*2+p*8]=vy[p*4+rec]*0.5;
	ez[rec1*2+p*8]=vz[p*4+rec2]*0.5;
	ez[1+rec1*2+p*8]=vz[p*4+rec]*0.5;
      }
  for(int i=0;i<4;i++)
    {/*
      b[i*2].cargue(0,0,1);		b[1+i*2].cargue(0,0,-1);

      b[i*2+8].cargue(0,-1,0);       	b[1+i*2+8].cargue(0,1,0);
      
      b[i*2+16].cargue(1,0,0);		b[1+i*2+16].cargue(-1,0,0);*/
      bx[i*2]=0;	by[i*2]=0;	bz[i*2]=1;
      bx[i*2+1]=0;	by[i*2+1]=0;	bz[i*2+1]=-1;
      bx[i*2+8]=0;	by[i*2+8]=-1;	bz[i*2+8]=0;
      bx[1+i*2+8]=0;	by[1+i*2+8]=1;	bz[1+i*2+8]=0;
      bx[i*2+16]=1;	by[i*2+16]=0;	bz[i*2+16]=0;
      bx[1+i*2+16]=-1;	by[1+i*2+16]=0;	bz[1+i*2+16]=0;
    }
 
}
void LatticeBoltzmann::ResizeDomain(int Lx0, int Ly0, int Lz0)
{
  Lx=Lx0;	Ly=Ly0;		Lz=Lz0;
  f = new double[sizef];	fnew = new double[sizef];
  f0 = new double[sizef0];	f0new = new double[sizef0];  
}
LatticeBoltzmann::~LatticeBoltzmann()
{
  delete f ;	delete fnew;
  delete f0;	delete f0new;  
}
double LatticeBoltzmann::rho(int ix, int iy, int iz, bool useNew)
{
  double sum=0;int k=0;int k0=0;
  for(int p=0;p<3;p++)
    for(int i=0;i<4;i++)      
      for(int j=0;j<2;j++){
	k=j + i*jS + p*jS*iS + iz*iS*jS*rS*pS + iy*iS*jS*rS*pS*Lz + ix*iS*jS*rS*pS*Lz*Ly;
 	sum += useNew? (*(fnew+k)) : (*(f+k));}
  k0=iz*rS + iy*rS*Lz + ix*rS*Lz*Ly;
  sum += (*f0+k0);
  return sum;
}
double LatticeBoltzmann::Dx(int ix, int iy, int iz, bool useNew)
{
  double sum=0;int k=0;
  for(int p=0;p<3;p++)
    for(int i=0;i<4;i++)      
      for(int j=0;j<2;j++){
	k=j + i*jS + p*jS*iS + iz*iS*jS*rS*pS + iy*iS*jS*rS*pS*Lz + ix*iS*jS*rS*pS*Lz*Ly;
	sum += useNew? (*(fnew+k))*ex[j+i*2+p*8] : (*(f+k))*ex[j+i*2+p*8];}
  return sum;
}
double LatticeBoltzmann::Dy(int ix, int iy, int iz, bool useNew)
{
  double sum=0;int k=0;
  for(int p=0;p<3;p++)
    for(int i=0;i<4;i++)      
      for(int j=0;j<2;j++){
	k=j + i*jS + p*jS*iS + iz*iS*jS*rS*pS + iy*iS*jS*rS*pS*Lz + ix*iS*jS*rS*pS*Lz*Ly;
	sum += useNew? (*(fnew+k))*ey[j+i*2+p*8] : (*(f+k))*ey[j+i*2+p*8] ;}
  return sum;
}
double LatticeBoltzmann::Dz(int ix, int iy, int iz, bool useNew)
{
  double sum=0;int k=0;
  for(int p=0;p<3;p++)
    for(int i=0;i<4;i++)      
      for(int j=0;j<2;j++){
	k=j + i*jS + p*jS*iS + iz*iS*jS*rS*pS + iy*iS*jS*rS*pS*Lz + ix*iS*jS*rS*pS*Lz*Ly;
	sum += useNew? (*(fnew+k))*ez[j+i*2+p*8] : (*(f+k))*ez[j+i*2+p*8] ;}		
  return sum;
}
double LatticeBoltzmann::Bx(int ix, int iy, int iz, bool useNew)
{
  double sum=0;int k=0;
  for(int p=0;p<3;p++)
    for(int i=0;i<4;i++)      
      for(int j=0;j<2;j++){
	k=j + i*jS + p*jS*iS + iS*jS*pS + iz*iS*jS*rS*pS + iy*iS*jS*rS*pS*Lz + ix*iS*jS*rS*pS*Lz*Ly;
	sum += useNew? (*(fnew+k))*bx[j+i*2+p*8] : (*(f+k))*bx[j+i*2+p*8] ;}		
  return sum;
}
double LatticeBoltzmann::By(int ix, int iy, int iz, bool useNew)
{
  double sum=0;int k=0;
  for(int p=0;p<3;p++)
    for(int i=0;i<4;i++)      
      for(int j=0;j<2;j++){
	k=j + i*jS + p*jS*iS + iS*jS*pS + iz*iS*jS*rS*pS + iy*iS*jS*rS*pS*Lz + ix*iS*jS*rS*pS*Lz*Ly;
	sum += useNew? (*(fnew+k))*by[j+i*2+p*8] : (*(f+k))*by[j+i*2+p*8] ;}    
  return sum;
}
double LatticeBoltzmann::Bz(int ix, int iy, int iz, bool useNew)
{
  double sum=0;int k=0;
  for(int p=0;p<3;p++)
    for(int i=0;i<4;i++)      
      for(int j=0;j<2;j++){
	k=j + i*jS + p*jS*iS + iS*jS*pS + iz*iS*jS*rS*pS + iy*iS*jS*rS*pS*Lz + ix*iS*jS*rS*pS*Lz*Ly;
	sum += useNew? (*(fnew+k))*bz[j+i*2+p*8] : (*(f+k))*bz[j+i*2+p*8] ;}	
  return sum;
}
double LatticeBoltzmann::Ex(double&Dx0, double&epsr)
{
  double sum=0;
  
  sum = Dx0*(1.0/epsr);
  return sum;
}
double LatticeBoltzmann::Ey(double&Dy0, double&epsr)
{
  double sum=0;
  
  sum = Dy0*(1.0/epsr);
  return sum;
}
double LatticeBoltzmann::Ez(double&Dz0, double&epsr)
{
  double sum=0;
  
  sum = Dz0*(1.0/epsr);
  return sum;
}
double LatticeBoltzmann::Hx(double&Bx0, double&mur)
{
  double sum=0;

  sum = Bx0*(1.0/mur);
  return sum;
}
double LatticeBoltzmann::Hy(double&By0, double&mur)
{
  double sum=0;

  sum = By0*(1.0/mur);
  return sum;
}
double LatticeBoltzmann::Hz(double&Bz0, double&mur)
{
  double sum=0;

  sum = Bz0*(1.0/mur);
  return sum;
}
double LatticeBoltzmann::Jx(double&Ex0, double&sigma)
{
  double sum=0;

  sum = Ex0*sigma;
  return sum; 
}
double LatticeBoltzmann::Jy(double&Ey0, double&sigma)
{
  double sum=0;

  sum = Ey0*sigma;
  return sum; 
}
double LatticeBoltzmann::Jz(double&Ez0, double&sigma)
{
  double sum=0;

  sum = Ez0*sigma;
  return sum; 
}
double LatticeBoltzmann::Jxp(double&Ex0, double&denominator)
{
  double sum=0;

  sum = Ex0*denominator;
  return sum; 
}
double LatticeBoltzmann::Jyp(double&Ey0, double&denominator)
{
  double sum=0;

  sum = Ey0*denominator;
  return sum; 
}
double LatticeBoltzmann::Jzp(double&Ez0, double&denominator)
{
  double sum=0;

  sum = Ez0*denominator;
  return sum; 
}
double LatticeBoltzmann::Exp(double&Ex0,double&Jxp0,double&factor)
{
  double sum=0;

  sum = Ex0-Jxp0*factor;
  return sum;
}
double LatticeBoltzmann::Eyp(double&Ey0,double&Jyp0,double&factor)
{
  double sum=0;

  sum = Ey0-Jyp0*factor;
  return sum;
}
double LatticeBoltzmann::Ezp(double&Ez0,double&Jzp0,double&factor)
{
  double sum=0;

  sum = Ez0-Jzp0*factor;
  return sum;
}
double LatticeBoltzmann::feq(double&epsr, double&mur, double&vJp, double&eEp, double&bB,int r)
{
  double f=0;
  if(r==0)
    {f = (0.0625*vJp) + (epsr*0.25*eEp) + (0.125*bB/(mur));}
  else if(r==1)
    {f = (0.0625*vJp) + (0.25*eEp) + (0.125*bB);}
  
  return f;  
}
double LatticeBoltzmann::feq0(double&rho0)
{
  return rho0;
}
void LatticeBoltzmann::Colisione(int &t)
{
  int k=0,k0=0,ix=0,iy=0,iz=0,r=0,j=0,i=0,p=0;
  double Epsr=0,Mur=0,Sigma=0,denominator=0,factor=0;
  double rho0=0;
  double Dx0,Bx0,Ex0,Hx0,Jx0,Jxp0,Exp0;
  double Dy0,By0,Ey0,Hy0,Jy0,Jyp0,Eyp0;
  double Dz0,Bz0,Ez0,Hz0,Jz0,Jzp0,Ezp0;
  double vJp, eEp, bB;
  
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++)
	{
	  Epsr=epsr(ix,iy,iz);		Mur=mur(ix,iy,iz);	     Sigma=sigma(ix,iy,iz);
	  denominator = Sigma/(1.0 + mu0*Sigma/(4.0*Epsr));
	  factor = mu0/(4.0*Epsr);
	  rho0=rho(ix,iy,iz,false);
	  Dx0=Dx(ix,iy,iz,false);Dy0=Dy(ix,iy,iz,false);Dz0=Dz(ix,iy,iz,false);
	  Bx0=Bx(ix,iy,iz,false);By0=By(ix,iy,iz,false);Bz0=Bz(ix,iy,iz,false);
	  Hx0=Hx(Bx0,Mur);Hy0=Hy(By0,Mur);Hz0=Hz(Bz0,Mur);
	  Ex0=Ex(Dx0,Epsr);Ey0=Ey(Dy0,Epsr);Ez0=Ez(Dz0,Epsr);
	  Jx0=Jx(Ex0,Sigma);Jy0=Jy(Ey0,Sigma);Jz0=Jz(Ez0,Sigma);
	  Jxp0=Jxp(Ex0,denominator);Jyp0=Jyp(Ey0,denominator);Jzp0=Jzp(Ez0,denominator);
	  Exp0=Exp(Ex0,Jxp0,factor);Eyp0=Eyp(Ey0,Jyp0,factor);Ezp0=Ezp(Ez0,Jzp0,factor);
	  
	  
  
	  for(r=0;r<2;r++)
	    for(p=0;p<3;p++)
	      for(i=0;i<4;i++)
		for(j=0;j<2;j++)
		  {
		    vJp = (vx[p*4+i]*Jxp0+vy[p*4+i]*Jyp0+vz[p*4+i]*Jzp0);
		    eEp = (ex[j+i*2+p*8]*Exp0+ey[j+i*2+p*8]*Eyp0+ez[j+i*2+p*8]*Ezp0);
		    bB = (bx[j+i*2+p*8]*Bx0+by[j+i*2+p*8]*By0+bz[j+i*2+p*8]*Bz0);
		    k = j + i*jS + p*jS*iS + r*iS*jS*pS + iz*iS*jS*rS*pS + iy*iS*jS*rS*pS*Lz + ix*iS*jS*rS*pS*Lz*Ly;
		    (*(fnew+k))=UmUtau*(*(f+k))+Utau*feq(Epsr,Mur,vJp,eEp,bB,r);
		    k0 = r + iz*rS + iy*rS*Lz + ix*rS*Lz*Ly;
		    (*(f0new+k0))=UmUtau*(*(f0+k0))+Utau*feq0(rho0);
		  }
	}
}
void LatticeBoltzmann::Adveccione(void)
{
  int k=0,k0=0,q=0,q0=0,jx=0,jy=0,jz=0;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int iz=0;iz<Lz;iz++)
	for(int r=0;r<2;r++)
	  for(int p=0;p<3;p++)
	    for(int i=0;i<4;i++)
	      for(int j=0;j<2;j++)
		{
		  if(ix==0 || ix==Lx-1){jx=ix;}
		  else{jx=ix+(int)V[0][i][p];};
		  if(iy==0 || iy==Ly-1){jy=iy;}
		  else{jy=iy+(int)V[1][i][p];};
		  if(iz==0 || iz==Lz-1){jz=iz;}
		  else{jz=iz+(int)V[2][i][p];};
		  k = j + i*jS + p*jS*iS + r*iS*jS*pS + iz*iS*jS*rS*pS + iy*iS*jS*rS*pS*Lz + ix*iS*jS*rS*pS*Lz*Ly;	       
		  q = j + i*jS + p*jS*iS + r*iS*jS*pS + jz*iS*jS*rS*pS + jy*iS*jS*rS*pS*Lz + jx*iS*jS*rS*pS*Lz*Ly;
		  k0 = r + iz*rS + iy*rS*Lz + ix*rS*Lz*Ly;
		  q0 = r + jz*rS + jy*rS*Lz + jx*rS*Lz*Ly;
		  (*(f+q)) = (*(fnew+k));
		  (*(f0+q0)) = (*(f0new+k0));
		}
}
void LatticeBoltzmann::Inicie(void)
{
  int k=0,k0=0;
  double Epsr,Mur;
  double rho0;
  double vJp, eEp, bB;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int iz=0;iz<Lz;iz++)
	  {
	    Epsr=epsr(ix,iy,iz);
	    Mur=mur(ix,iy,iz);
	    rho0=0;
	    vJp=eEp=bB=0;
	    for(int r=0;r<2;r++)
	      for(int p=0;p<3;p++)
		for(int i=0;i<4;i++)
		  for(int j=0;j<2;j++)
		    {
		      k = j + i*jS + p*jS*iS + r*iS*jS*pS + iz*iS*jS*rS*pS + iy*iS*jS*rS*pS*Lz + ix*iS*jS*rS*pS*Lz*Ly;
		      (*(f+k)) = feq(Epsr,Mur,vJp,eEp,bB,r);
		      k0 = r + iz*rS + iy*rS*Lz + ix*rS*Lz*Ly;
		      (*(f0+k0)) = feq0(rho0);
		    }
	  }
}
void LatticeBoltzmann::Imprimase(const char* fileName,bool useNew)
{
  ofstream outputFile(fileName);
  double Dx0, Bx0;					     
  int ix0 = Lx/2,iy0=Ly/2,iz0 = Lz/2;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
    for(int iz=0;iz<Lz;iz++)
    {
      Dx0=Dx(ix,iy,iz,useNew);	
      Bx0=Bx(ix,iy,iz,useNew);
      outputFile
	<< ix << "\t"
	<< iy << "\t"
	<< iz << "\t"
	<< Bx0 << "\t"
	<< Dx0 << "\n";	      
    }
   outputFile.close();
}
