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
#include "Distribution.h"
#include "Vector.h"
//using namespace std;

const int Lx=1;
const int Ly=1;
const int Lz=200;

//Constates e0 y mu0 (INTOCABLE).
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
  
  double V[3][4][3]; vector3D v[3][4],v0;// V[0][i][p]=V^p_ix,  V[1][i][p]=V^p_iy, V[2][i][p]=V^p_iz
  vector3D e[2][4][3]; //e[0][j][i][p]=e^p_ijx,e[1][j][i][p]=e^p_ijy,e[2][j][i][p]=e^p_ijz,
  vector3D b[2][4][3]; //b[0][j][i][p]=b^p_ijx,b[1][j][i][p]=b^p_ijy,b[2][j][i][p]=b^p_ijz,

  Distribution f = Distribution(Lx,Ly,Lz,false);	Distribution fnew = Distribution(Lx,Ly,Lz,false);
  Distribution f0 = Distribution(Lx,Ly,Lz,true);	Distribution f0new = Distribution(Lx,Ly,Lz,true);
  
  //double f[Lx][Ly][Lz][2][2][4][3], fnew[Lx][Ly][Lz][2][2][4][3]; // f[ix][iy][iz][r][j][i][p]
  //double f0[Lx][Ly][Lz][2], f0new[Lx][Ly][Lz][2]; // f[ix][iy][iz][r]
public:
  LatticeBoltzmann(void);
  double rho(int ix, int iy, int iz, bool useNew);
  //campos
  vector3D D(int ix, int iy, int iz, bool useNew);
  vector3D B(int ix, int iy, int iz, bool useNew);
  vector3D E(vector3D&D0, double&epsr);
  vector3D H(vector3D&B0, double&mur);
  vector3D J(vector3D&E0, double&sigma);
  //campos auxiliares
  vector3D Jp(vector3D&E0, double&denominator);
  vector3D Ep(vector3D&E0, vector3D&Jp0, double&factor);
  //constantes dielectricas relativas
  double epsr(int ix, int iy, int iz){return 1.0;}; 
  double mur(int ix, int iy, int iz){return 1.0;};
  double sigma(int ix, int iy, int iz){return 0.0;};
  //funciones de equilibrio
  double feq(double&epsr,double&mur,vector3D&Jp0,vector3D&Ep0,vector3D&B0,int r,int j,int i,int p);
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
  v0.cargue(V0[0],V0[1],V0[2]);
  for(int p=0;p<3;p++)
    for(int i=0;i<4;i++){
      v[p][i].cargue(V[0][i][p],V[1][i][p],V[2][i][p]);
  }
  
  for(int p=0;p<3;p++)
    for(int i=0;i<4;i++)
      {
	rec=i%4; rec1=(i+1)%4; rec2=(i+2)%4;

	e[0][rec1][p]=v[p][rec2]*0.5;
	e[1][rec1][p]=v[p][rec]*0.5;
      }
  for(int i=0;i<4;i++)
    {
      b[0][i][0].cargue(0,0,1);		b[1][i][0].cargue(0,0,-1);

      b[0][i][1].cargue(0,-1,0);       	b[1][i][1].cargue(0,1,0);
      
      b[0][i][2].cargue(1,0,0);		b[1][i][2].cargue(-1,0,0);
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
vector3D LatticeBoltzmann::D(int ix, int iy, int iz, bool useNew)
{
  vector3D sum; sum.cargue(0,0,0);
  
    for(int i=0;i<4;i++)
      for(int p=0;p<3;p++)
      	for(int j=0;j<2;j++)	
	  sum += useNew? fnew.function(0,j,p,i,ix,iy,iz)*e[j][i][p] : f.function(0,j,p,i,ix,iy,iz)*e[j][i][p];
  return sum;
  
}
vector3D LatticeBoltzmann::B(int ix, int iy, int iz, bool useNew)
{
  vector3D sum; sum.cargue(0,0,0);
  
    for(int i=0;i<4;i++)
      for(int p=0;p<3;p++)      
	for(int j=0;j<2;j++)	
	  sum += useNew? fnew.function(1,j,p,i,ix,iy,iz)*b[j][i][p] : f.function(1,j,p,i,ix,iy,iz)*b[j][i][p];
  return sum;
  
}
vector3D LatticeBoltzmann::E(vector3D&D0, double&epsr)
{
  vector3D sum; sum.cargue(0,0,0);
  
  sum = D0*(1.0/epsr);
  return sum;
  
}
vector3D LatticeBoltzmann::H(vector3D&B0, double&mur)
{
  vector3D sum; sum.cargue(0,0,0);

  sum = B0*(1.0/mur);
  return sum;
  
}
vector3D LatticeBoltzmann::J(vector3D&E0, double&sigma)
{
  vector3D sum; sum.cargue(0,0,0);

  sum = E0*sigma;
  return sum; 
  
}
vector3D LatticeBoltzmann::Jp(vector3D&E0, double&denominator)
{
  //double denominator = sigma(ix,iy,iz)/(1.0 + mu0*sigma(ix,iy,iz)/(4.0*epsr(ix,iy,iz)));
  vector3D sum; sum.cargue(0,0,0);

  sum = E0*denominator;
  return sum; 
}
vector3D LatticeBoltzmann::Ep(vector3D&E0,vector3D&Jp0,double&factor)
{
  //double factor = mu0/(4.0*epsr(ix,iy,iz));
  vector3D sum; sum.cargue(0,0,0);

  sum = E0-Jp0*factor;
  return sum;
}
double LatticeBoltzmann::feq(double&epsr,double&mur,vector3D&Jp0,vector3D&Ep0,vector3D&B0,int r,int j,int i,int p)
{
  double f=0, vJp=0, eEp=0, bB=0;
  vJp = (v[p][i]*Jp0);	eEp = (e[j][i][p]*Ep0);	bB = (b[j][i][p]*B0);
  /*for(int x=0;x<3;x++)
    {
      vJp += V[x][i][p]*Jp0[x];
      eEp += e[x][j][i][p]*Ep0[x];
      bB  += b[x][j][i][p]*B0[x];
      }*/
  
  //OJO. la funcion de equilibrio esta escrita en terminos de epsilon y mu relativos!!
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
  int ix=0,iy=0,iz=0,r=0,j=0,i=0,p=0;
  double Epsr=0,Mur=0,Sigma=0,denominator=0,factor=0;
  double rho0=0;
  vector3D D0,B0,E0,H0,J0,Jp0,Ep0;
  
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++)
	{
	  Epsr=epsr(ix,iy,iz);		Mur=mur(ix,iy,iz);	     Sigma=sigma(ix,iy,iz);
	  denominator = Sigma/(1.0 + mu0*Sigma/(4.0*Epsr));
	  factor = mu0/(4.0*Epsr);
	  rho0=rho(ix,iy,iz,false);	 D0=D(ix,iy,iz,false);	     B0=B(ix,iy,iz,false);
	  H0=H(B0,Mur);			 E0=E(D0,Epsr);
	  J0=J(E0,Sigma);		 Jp0=Jp(E0,denominator);
	  Ep0=Ep(E0,Jp0,factor);
	  
	  for(p=0;p<3;p++)
	    for(i=0;i<4;i++)
	      for(j=0;j<2;j++)
		for(r=0;r<2;r++)
		{ 
		  fnew(r,j,p,i,ix,iy,iz)=UmUtau*f.function(r,j,p,i,ix,iy,iz)+Utau*feq(Epsr,Mur,Jp0,Ep0,B0,r,j,i,p);
		  f0new(r,ix,iy,iz)=UmUtau*f0.function(r,ix,iy,iz)+Utau*feq0(rho0);
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
		  if(ix==0 || ix==Lx-1){jx=ix;}
		  else{jx=ix+(int)V[0][i][p];};
		  if(iy==0 || iy==Ly-1){jy=iy;}
		  else{jy=iy+(int)V[1][i][p];};
		  if(iz==0 || iz==Lz-1){jz=iz;}
		  else{jz=iz+(int)V[2][i][p];};
		  	       
		  f(r,j,p,i,jx,jy,jz) = fnew.function(r,j,p,i,ix,iy,iz);
		  f0(r,ix,iy,iz) = f0new.function(r,ix,iy,iz);
		}
}
void LatticeBoltzmann::Inicie(void)
{
  double Epsr,Mur;
  double rho0; vector3D D0,B0,E0,H0,J0,Jp0,Ep0;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int iz=0;iz<Lz;iz++)
	  {
	    Epsr=epsr(ix,iy,iz);
	    Mur=mur(ix,iy,iz);
	    rho0=0;
	    D0.cargue(0,0,0);
	    B0.cargue(0,0,0);
	    E0.cargue(0,0,0);
	    H0.cargue(0,0,0);
	    J0.cargue(0,0,0);
	    Jp0.cargue(0,0,0);
	    Ep0.cargue(0,0,0);

	    for(int p=0;p<3;p++)
	      for(int i=0;i<4;i++)
		for(int j=0;j<2;j++)
		  for(int r=0;r<2;r++)
		    {
		      f(r,j,p,i,ix,iy,iz) = feq(Epsr,Mur,Jp0,Ep0,B0,r,j,i,p);
		      f0(r,ix,iy,iz) = feq0(rho0);
		      
		    }
	  }
}
void LatticeBoltzmann::Imprimase(const char* fileName,bool useNew)
{
  ofstream outputFile(fileName);
  vector3D D0, B0;					     
  int ix0 = Lx/2,iy0=Ly/2,iz0 = Lz/2;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
    for(int iz=0;iz<Lz;iz++)
    {
      D0=D(ix,iy,iz,useNew);	
      B0=B(ix,iy,iz,useNew);
      outputFile
	<< ix << "\t"
	<< iy << "\t"
	<< iz << "\t"
	<< B0.x() << "\t"
	<< B0.y() << "\t"
	<< B0.z() << "\t"
	<< D0.x() << "\t"
	<< D0.y() << "\t"
	<< D0.z() << "\n";	      
    }
   outputFile.close();
}
