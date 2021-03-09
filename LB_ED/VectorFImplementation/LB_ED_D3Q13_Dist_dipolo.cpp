/*
Reproduccion LBM para ED

Notas:
 - Solo funciona para simulaciones en 1D (2D no se ha probado. 3D requiere arreglos de memoria dinamica).
 - La funcion de equilibrio debe estar escrita en terminos de epsilon y mu relativos.
 - Las cantidades macroscopicas vectoriales son arreglos dinamicos, para evitar hacer 3 funciones (componentes).
 - Para Efecto skin se recomienda sigma0 = 0.1, omega = M_PI*0.02. Imprimir en t=612, t=630, t=655.
 */
//#include <iostream>
#include <fstream>
#include <cmath>
#include "Distribution.h"
//using namespace std;

const int Lx=1;
const int Ly=50;
const int Lz=50;

//Constates e0 y mu0 (INTOCABLE).
const double C=1.0/sqrt(2.0);
const double eps0=1.0;
const double mu0=2.0;

const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;
	       
class LatticeBoltzmann{
private:

  double V0[3]={0.0,0.0,0.0};
  double e0[3]={0.0,0.0,0.0};
  double b0[3]={0.0,0.0,0.0};
  
  double V[3][4][3];// V[0][i][p]=V^p_ix,  V[1][i][p]=V^p_iy, V[2][i][p]=V^p_iz
  double e[3][2][4][3]; //e[0][j][i][p]=e^p_ijx,e[1][j][i][p]=e^p_ijy,e[2][j][i][p]=e^p_ijz,
  double b[3][2][4][3]; //b[0][j][i][p]=b^p_ijx,b[1][j][i][p]=b^p_ijy,b[2][j][i][p]=b^p_ijz,

  Distribution f = Distribution(Lx,Ly,Lz,false);	Distribution fnew = Distribution(Lx,Ly,Lz,false);
  Distribution f0 = Distribution(Lx,Ly,Lz,true);	Distribution f0new = Distribution(Lx,Ly,Lz,true);
  
  //double f[Lx][Ly][Lz][2][2][4][3], fnew[Lx][Ly][Lz][2][2][4][3]; // f[ix][iy][iz][r][j][i][p]
  //double f0[Lx][Ly][Lz][2], f0new[Lx][Ly][Lz][2]; // f[ix][iy][iz][r]
public:
  LatticeBoltzmann(void);
  double rho(int ix, int iy, int iz, bool useNew);
  //campos
  double *D(int ix, int iy, int iz, bool useNew);
  double *B(int ix, int iy, int iz, bool useNew);
  double *E(double*D0,int ix, int iy, int iz, bool useNew);
  double *H(double*B0,int ix, int iy, int iz, bool useNew);
  double *J(double*E0,int ix, int iy, int iz, bool useNew);
  //campos auxiliares
  double *Jp(double*E0,int ix, int iy, int iz, bool useNew);
  double *Ep(double*E0,double*Jp0,int ix, int iy, int iz, bool useNew);
  //constantes dielectricas relativas
  double epsr(int ix, int iy, int iz){return 1.0;}; 
  double mur(int ix, int iy, int iz){return 1.0;};
  double sigma(int ix, int iy, int iz){return 0.0;};
  //funciones de equilibrio
  double feq(int ix,int iy,int iz,double *Jp0,double *Ep0,double *B0,int r,int j,int i,int p);
  double feq0(int ix, int iy, int iz, double&rho0,int r);
  void Colisione(double&rho0,double*D0,double*B0,double*H0,double*E0,double*J0,double*Jp0,double*Ep0);
  void Adveccione(void);
  void Inicie(double&rho0,double*D0,double*B0,double*H0,double*E0,double*J0,double*Jp0,double*Ep0);
  void ImponerCampos(double&rho0,double*D0,double*B0,double*H0,double*E0,double*J0,double*Jp0,double*Ep0,int t);
  void Imprimase(const char * FileName,double*D0,double*E0,bool useNew);
};
LatticeBoltzmann::LatticeBoltzmann(void){
  
  int imp=0,rec=0,rec1=0,rec2=0;
  double pi4=M_PI*0.25, sq2=sqrt(2.0);
  
  for(int i=0;i<4;i++)
    {
      imp=2*i+1;
      V[0][i][0]=cos(imp*pi4)*sq2;	V[1][i][0]=sin(imp*pi4)*sq2;	V[2][i][0]=0;	
      V[0][i][1]=cos(imp*pi4)*sq2;	V[1][i][1]=0;			V[2][i][1]=sin(imp*pi4)*sq2;	
      V[0][i][2]=0;			V[1][i][2]=cos(imp*pi4)*sq2;	V[2][i][2]=sin(imp*pi4)*sq2;
          };
  for(int p=0;p<3;p++)
    for(int i=0;i<4;i++)
      {
	rec=i%4; rec1=(i+1)%4; rec2=(i+2)%4;

	e[0][0][rec1][p]=V[0][rec2][p]*0.5;	e[1][0][rec1][p]=V[1][rec2][p]*0.5;	e[2][0][rec1][p]=V[2][rec2][p]*0.5;
	e[0][1][rec1][p]=V[0][rec][p]*0.5;	e[1][1][rec1][p]=V[1][rec][p]*0.5;	e[2][1][rec1][p]=V[2][rec][p]*0.5;
      }
  for(int i=0;i<4;i++)
    {
      b[0][0][i][0]=0;			b[1][0][i][0]=0;		b[2][0][i][0]=1;
      b[0][1][i][0]=0;			b[1][1][i][0]=0;		b[2][1][i][0]=-1;
      
      b[0][0][i][1]=0;			b[1][0][i][1]=-1;		b[2][0][i][1]=0;
      b[0][1][i][1]=0;			b[1][1][i][1]=1;		b[2][1][i][1]=0;
      
      b[0][0][i][2]=1;			b[1][0][i][2]=0;		b[2][0][i][2]=0;
      b[0][1][i][2]=-1;			b[1][1][i][2]=0;		b[2][1][i][2]=0;
    }
 
}

double LatticeBoltzmann::rho(int ix, int iy, int iz, bool useNew)
{
  double sum=0;
  
  for(int i=0;i<4;i++)
    for(int p=0;p<3;p++)
      for(int j=0;j<2;j++)	
	sum += useNew? fnew.function(0,j,p,i,ix,iy,iz) : f.function(0,j,p,i,ix,iy,iz);
  sum += f0(0,ix,iy,iz);
  return sum;
}
double *LatticeBoltzmann::D(int ix, int iy, int iz, bool useNew)
{
  double *sum=new double[3];
  sum[0]=0;sum[1]=0;sum[2]=0;
  for(int x=0;x<3;x++)
    for(int i=0;i<4;i++)
      for(int p=0;p<3;p++)
      	for(int j=0;j<2;j++)	
	  sum[x] += useNew? fnew.function(0,j,p,i,ix,iy,iz)*e[x][j][i][p] : f.function(0,j,p,i,ix,iy,iz)*e[x][j][i][p];
  return sum;
  
}
double *LatticeBoltzmann::B(int ix, int iy, int iz, bool useNew)
{
  double *sum=new double[3];
  sum[0]=0;sum[1]=0;sum[2]=0;
  for(int x=0;x<3;x++)
    for(int i=0;i<4;i++)
      for(int p=0;p<3;p++)      
	for(int j=0;j<2;j++)	
	  sum[x] += useNew? fnew.function(1,j,p,i,ix,iy,iz)*b[x][j][i][p] : f.function(1,j,p,i,ix,iy,iz)*b[x][j][i][p];
  return sum;
  
}
double *LatticeBoltzmann::E(double*D0,int ix, int iy, int iz, bool useNew)
{
  double *sum=new double[3];
  for(int x=0;x<3;x++)
    sum[x] = D0[x]/epsr(ix,iy,iz);
  return sum;
  
}
double *LatticeBoltzmann::H(double*B0,int ix, int iy, int iz, bool useNew)
{
  double *sum=new double[3];
  for(int x=0;x<3;x++)
    sum[x] = B0[x]/mur(ix,iy,iz);
  return sum;
  
}
double *LatticeBoltzmann::J(double*E0,int ix, int iy, int iz, bool useNew)
{
  double *sum=new double[3];
  for(int x=0;x<3;x++)
    sum[x] = E0[x]*sigma(ix,iy,iz);
  return sum; 
  
}
double *LatticeBoltzmann::Jp(double*E0,int ix, int iy, int iz, bool useNew)
{
  double denominator = sigma(ix,iy,iz)/(1.0 + mu0*sigma(ix,iy,iz)/(4.0*epsr(ix,iy,iz)));
  double *sum = new double[3];
  for(int x=0;x<3;x++)
    sum[x] = E0[x]*denominator;
  return sum; 
}
double *LatticeBoltzmann::Ep(double*E0,double*Jp0,int ix, int iy, int iz, bool useNew)
{
  double factor = mu0/(4.0*epsr(ix,iy,iz));
  double *sum = new double[3];
  for(int x=0;x<3;x++)
    sum[x] = E0[x]-Jp0[x]*factor;
  return sum;
}
double LatticeBoltzmann::feq(int ix,int iy, int iz,double *Jp0,double *Ep0,double *B0,int r,int j,int i,int p)
{
  double f=0, vJp=0, eEp=0, bB=0;
  for(int x=0;x<3;x++)
    {
      vJp += V[x][i][p]*Jp0[x];
      eEp += e[x][j][i][p]*Ep0[x];
      bB  += b[x][j][i][p]*B0[x];
    }
  
  //OJO. la funcion de equilibrio esta escrita en terminos de epsilon y mu relativos!!
  if(r==0)
    {f = (0.0625*vJp) + (epsr(ix,iy,iz)*0.25*eEp) + (0.125*bB/(mur(ix,iy,iz)));}
  else if(r==1)
    {f = (0.0625*vJp) + (0.25*eEp) + (0.125*bB);}
  
  return f;  
}
double LatticeBoltzmann::feq0(int ix, int iy, int iz, double&rho0,int r)
{
  return rho0;
}
void LatticeBoltzmann::Colisione(double&rho0,double*D0,double*B0,double*H0,double*E0,double*J0,double*Jp0,double*Ep0)
{
  int ix=0,iy=0,iz=0
    ,r=0,j=0,i=0,p=0;
  
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++)
	{
	  
	  rho0=rho(ix,iy,iz,false);	 D0=D(ix,iy,iz,false);	     B0=B(ix,iy,iz,false);
	  H0=H(B0,ix,iy,iz,false);	 E0=E(D0,ix,iy,iz,false);
	  J0=J(E0,ix,iy,iz,false);	 Jp0=Jp(E0,ix,iy,iz,false);
	  Ep0=Ep(E0,Jp0,ix,iy,iz,false);
	  
	  for(p=0;p<3;p++)
	    for(i=0;i<4;i++)
	      for(j=0;j<2;j++)
		for(r=0;r<2;r++)
		{ 
		  //fnew[ix][iy][iz][r][j][i][p]=UmUtau*f[ix][iy][iz][r][j][i][p]+Utau*feq(ix,iy,iz,Jp0,Ep0,B0,r,j,i,p);
		  //f0new[ix][iy][iz][r]=UmUtau*f0[ix][iy][iz][r]+Utau*feq0(ix,iy,iz,rho0,r);
		  fnew(r,j,p,i,ix,iy,iz)=UmUtau*f.function(r,j,p,i,ix,iy,iz)+Utau*feq(ix,iy,iz,Jp0,Ep0,B0,r,j,i,p);
		  f0new(r,ix,iy,iz)=UmUtau*f0.function(r,ix,iy,iz)+Utau*feq0(ix,iy,iz,rho0,r);
		}
	}
}
void LatticeBoltzmann::Adveccione(void)
{
  int jx=0,jy=0,jz=0;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int iz=0;iz<Lz;iz++)
	for(int p=0;p<3;p++)
	  for(int i=0;i<4;i++)
	    for(int j=0;j<2;j++)
	      for(int r=0;r<2;r++)
		{
		  jx=(ix+(int)V[0][i][p]+Lx)%Lx;	jy=(iy+(int)V[1][i][p]+Ly)%Ly;	jz=(iz+(int)V[2][i][p]+Lz)%Lz;
		  //f[jx][jy][jz][r][j][i][p] = fnew[ix][iy][iz][r][j][i][p];
		  //f0[jx][jy][jz][r] = f0new[ix][iy][iz][r];
		  f(r,j,p,i,jx,jy,jz) = fnew.function(r,j,p,i,ix,iy,iz);
		  f0(r,jx,jy,jz) = f0new.function(r,ix,iy,iz);
		}
}
void LatticeBoltzmann::Inicie(double&rho0,double*D0,double*B0,double*H0,double*E0,double*J0,double*Jp0,double*Ep0)
{
  //double Eo=0.001,Bo=Eo/C,alp=0.01;
  //int iz0=40;
  
  D0[1]=0;	D0[2]=0;	B0[0]=0;	B0[2]=0;
  
	    
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int iz=0;iz<Lz;iz++)
	  {
	    D0[0] = 0.0;//Eo*epsr(ix,iy,iz)*exp(-alp*(iz-iz0)*(iz-iz0));	
	    B0[1] = 0.0;//Bo*exp(-alp*(iz-iz0)*(iz-iz0));
	    rho0=rho(ix,iy,iz,false);
	    E0=E(D0,ix,iy,iz,false);	H0=H(B0,ix,iy,iz,false);
	    J0=J(E0,ix,iy,iz,false);	
	    Jp0=Jp(E0,ix,iy,iz,false);	Ep0=Ep(E0,Jp0,ix,iy,iz,false);
	    
	    for(int p=0;p<3;p++)
	      for(int i=0;i<4;i++)
		for(int j=0;j<2;j++)
		  for(int r=0;r<2;r++)
		    {
		      f(r,j,p,i,ix,iy,iz) = feq(ix,iy,iz,Jp0,Ep0,B0,r,j,i,p);
		      f0(r,ix,iy,iz) = feq0(ix,iy,iz,rho0,r);
		      
		    }
	  }
}
void LatticeBoltzmann::ImponerCampos(double&rho0,double*D0,double*B0,double*H0,double*E0,double*J0,double*Jp0,double*Ep0,int t)
{
  double Jo=0.0001,T=12.5,alp=0.75,
    Jop=0,sine=0,denominator=0,factor=0;
  int ix0=Lx/2,iy0=Ly/2,iz0=Lz/2;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int iz=0;iz<Lz;iz++)
	{
	  
	  rho0=rho(ix,iy,iz,false);
	  D0=D(ix,iy,iz,false);		B0=B(ix,iy,iz,false);
	  E0=E(D0,ix,iy,iz,false);	H0=H(B0,ix,iy,iz,false);
	  //dipole parameters
	  Jop = Jo*exp(-alp*((iy-iy0)*(iy-iy0)+(iz-iz0)*(iz-iz0)));
	  sine = Jop*sin(2*M_PI*t/T);
	  denominator = 1.0/(1.0 + mu0/(4.0*epsr(ix,iy,iz)));
	  factor = mu0/(4.0*epsr(ix,iy,iz));
	  //impose fields
	  J0[0] = sine;
	  Jp0[0] = J0[0]*denominator;
	  Ep0[0] = E0[0] - Jp0[0]*factor;
	  //eq function
	  for(int p=0;p<3;p++)
	    for(int i=0;i<4;i++)
	      for(int j=0;j<2;j++)
		for(int r=0;r<2;r++)
		  {
		    fnew(r,j,p,i,ix,iy,iz) = feq(ix,iy,iz,Jp0,Ep0,B0,r,j,i,p);
		    f0new(r,ix,iy,iz) = feq0(ix,iy,iz,rho0,r);
		  }
	}
}
void LatticeBoltzmann::Imprimase(const char * FileName,double*D0,double*E0,bool useNew)
{
  ofstream outputFile(FileName);
  double*B0=nullptr;
  //double Eo = 0.001;
  
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int iz=0;iz<Lz;iz++)
	{
	  D0=D(ix,iy,iz,useNew);	E0=E(D0,ix,iy,iz,useNew);
	  B0=B(ix,iy,iz,useNew);
	  outputFile
	    << iy << "\t"
	    << iz << "\t"	 
	    << B0[1] << "\n";
	}
  outputFile.close();
}
int main(int argc, char * argv[])
{
  if(argc != 3)
    {
      cout << "Wrong parameters!! ./a.out time fileName";
      return -1;
    }
  
  LatticeBoltzmann LB = LatticeBoltzmann();
  //Es conveniente inicializar los arreglos dinamicos una sola vez, en lugar de hacerlo en cada paso de tiempo.
  double rho0=0;
  double* D0=new double[3];	D0[0]=0;D0[1]=0;D0[2]=0;
  double* B0=new double[3];	B0[0]=0;B0[1]=0;B0[2]=0;
  double* H0=new double[3];	H0[0]=0;H0[1]=0;H0[2]=0;
  double* E0=new double[3];	E0[0]=0;E0[1]=0;E0[2]=0;
  double* Ep0=new double[3];	Ep0[0]=0;Ep0[1]=0;Ep0[2]=0;
  double* J0=new double[3];	J0[0]=0;J0[1]=0;J0[2]=0;
  double* Jp0=new double[3];	Jp0[0]=0;Jp0[1]=0;Jp0[2]=0;
  
  int t=0,tmax=atoi(argv[1]);

  
  LB.Inicie(rho0,D0,B0,H0,E0,J0,Jp0,Ep0);
  for(t=0;t<tmax;t++)
    {
      LB.Colisione(rho0,D0,B0,H0,E0,J0,Jp0,Ep0);
      LB.ImponerCampos(rho0,D0,B0,H0,E0,J0,Jp0,Ep0,t);
      LB.Adveccione();
    }
  LB.Imprimase(argv[2],D0,E0,true);
  
  delete[] D0;	delete[] B0;	delete[] H0;	delete[] E0;
  delete[] Ep0;	delete[] J0;	delete[] Jp0;
  
  return 0;
}
