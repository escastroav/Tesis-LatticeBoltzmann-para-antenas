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
#include <GL/glew.h>
#include <GL/glut.h>
#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>

using namespace std;

#define Lx 128
#define Ly 128
#define Lz 128

#define Nx 32
#define Ny 32
#define Nz 32

const int Mx = (Lx+Nx-1)/Nx;
const int My = (Ly+Ny-1)/Ny;
const int Mz = (Lz+Nz-1)/Nz;

#define rS 2
#define jS 2
#define pS 3
#define iS 4

int sizef = Lx*Ly*Lz*rS*jS*pS*iS;
int sizef0 = Lx*Ly*Lz*rS;


const float C=1.0/sqrt(2.0);
const float eps0=1.0;
const float mu0=2.0;

const float tau=0.5;
const float Utau=1.0/tau;
const float UmUtau=1-Utau;
//constantes device
__constant__ float d_taus[3];//0->tau; 1->Utau; 2->UmUtau
__constant__ float d_c_e_m[3];//0->C;1->eps0;2->mu0

__constant__ float d_vx[12];
__constant__ float d_vy[12];
__constant__ float d_vz[12];

__constant__ float d_ex[24];
__constant__ float d_ey[24];
__constant__ float d_ez[24];

__constant__ float d_bx[24];
__constant__ float d_by[24];
__constant__ float d_bz[24];
//kernels device
__device__ float d_sigma(int ix, int iy, int iz){
  return 0.0;
}
__device__ float d_epsr(int ix, int iy, int iz){
  return 1.0;
}
__device__ float d_mur(int ix, int iy, int iz){
  return 1.0;
}
__device__ float d_rho(float * f, float * f0){
  float sum=0; int k=0; int k0=0;
  for(int p=0;p<3;p++)
    for(int i=0;i<4;i++)      
      for(int j=0;j<2;j++){
	k=j + i*jS + p*jS*iS;
 	sum += (*(f+k));}
  sum += (*f0);
  return sum;
}
__device__ float d_Dx(float * f)
{
  float sum=0;int k=0;
  for(int p=0;p<3;p++)
    for(int i=0;i<4;i++)      
      for(int j=0;j<2;j++){
	k=j + i*jS + p*jS*iS;
	sum += (*(f+k))*d_ex[j+i*2+p*8];}
  return sum;
}
__device__ float d_Dy(float * f)
{
  float sum=0;int k=0;
  for(int p=0;p<3;p++)
    for(int i=0;i<4;i++)      
      for(int j=0;j<2;j++){
	k=j + i*jS + p*jS*iS;
	sum += (*(f+k))*d_ey[j+i*2+p*8];}
  return sum;
}
__device__ float d_Dz(float * f)
{
  float sum=0;int k=0;
  for(int p=0;p<3;p++)
    for(int i=0;i<4;i++)      
      for(int j=0;j<2;j++){
	k=j + i*jS + p*jS*iS;
	sum += (*(f+k))*d_ez[j+i*2+p*8];}
  return sum;
}
__device__ float d_Bx(float * f)
{
  float sum=0;int k=0;
  for(int p=0;p<3;p++)
    for(int i=0;i<4;i++)      
      for(int j=0;j<2;j++){
	k=j + i*jS + p*jS*iS + iS*jS*pS;
	sum += (*(f+k))*d_bx[j+i*2+p*8];}
  return sum;
}
__device__ float d_By(float * f)
{
  float sum=0;int k=0;
  for(int p=0;p<3;p++)
    for(int i=0;i<4;i++)      
      for(int j=0;j<2;j++){
	k=j + i*jS + p*jS*iS + iS*jS*pS;
	sum += (*(f+k))*d_by[j+i*2+p*8];}
  return sum;
}
__device__ float d_Bz(float * f)
{
  float sum=0;int k=0;
  for(int p=0;p<3;p++)
    for(int i=0;i<4;i++)      
      for(int j=0;j<2;j++){
	k=j + i*jS + p*jS*iS + iS*jS*pS;
	sum += (*(f+k))*d_bz[j+i*2+p*8];}
  return sum;
}
__device__ float d_feq(float vJp, float eEp, float bB,int r){
  float f=0;float epsr = 1.0, mur=1.0;
  if(r==0)
    {f = (0.0625*vJp) + (epsr*0.25*eEp) + (0.125*bB/(mur));}
  else if(r==1)
    {f = (0.0625*vJp) + (0.25*eEp) + (0.125*bB);}
  
  return f;
}
//kernels principales
__global__ void d_Colisione(float * d_f, float * d_fnew, float * d_f0, float * d_f0new, int t)
{
  int dist_f = 48;//iS*jS*pS+2*iS*jS+3*jS+1;
  int dist_f0 = 2; 
  int ix = blockIdx.x*blockDim.x+threadIdx.x;
  int iy = blockIdx.y*blockDim.y+threadIdx.y;
  int iz = blockIdx.z*blockDim.z+threadIdx.z;
  float *f, *fnew, *f0, *f0new;
  f = d_f + dist_f*(ix*Ly*Lz+iy*Lz+iz); f0 = d_f0 + (ix*Ly*Lz+iy)*dist_f0;
  fnew = d_fnew + dist_f*(ix*Ly*Lz+iy*Lz+iz); f0new = d_f0new + (ix*Ly*Lz+iy*Lz+iz)*dist_f0;
  float Sigma, Epsr, Mur;
  Sigma = d_sigma(ix,iy,iz); Epsr=d_epsr(ix,iy,iz); Mur=d_mur(ix,iy,iz);
  float denominator = Sigma/(1.0 + d_c_e_m[2]*Sigma/(4.0*Epsr));
  float factor = d_c_e_m[2]/(4.0*Epsr);
  float rho = d_rho(f,f0);
  float Dx, Dy, Dz, Bx, By, Bz, Eop;
  Dx = d_Dx(f);	Bx = d_Bx(f);
  Dy = d_Dy(f);	By = d_By(f);
  Dz = d_Dz(f);	Bz = d_Bz(f);
  float Ex, Ey, Ez, Hx, Hy, Hz;
  Ex = Dx/Epsr;	Hx = Bx/Mur;
  Ey = Dy/Epsr;	Hy = By/Mur;
  Ez = Dz/Epsr;	Hz = Bz/Mur;
  float  Jxp, Jyp, Jzp, Exp, Eyp, Ezp;
  Eop = 0.0001*exp(-0.25*((ix-64)*(ix-64)+(iy-64)*(iy-64)+(iz-64)*(iz-64)));
  Jzp = Eop*sin(2*M_PI*t/25);
  Jxp = Ex*denominator;
  Jyp = Ez*denominator;
  //Jzp = Ez*denominator;
  Exp = Ex - factor*Jxp;
  Eyp = Ey - factor*Jyp;
  Ezp = Ez - factor*Jzp;
  float vJp, eEp, bB;
  int k, k0;
  for(int r=0;r<2;r++)
    for(int p=0;p<3;p++)
      for(int i=0;i<4;i++)
	for(int j=0;j<2;j++)
	  {
	    vJp = (d_vx[p*4+i]*Jxp+d_vy[p*4+i]*Jyp+d_vz[p*4+i]*Jzp);
	    eEp = (d_ex[j+i*2+p*8]*Exp+d_ey[j+i*2+p*8]*Eyp+d_ez[j+i*2+p*8]*Ezp);
	    bB = (d_bx[j+i*2+p*8]*Bx+d_by[j+i*2+p*8]*By+d_bz[j+i*2+p*8]*Bz);
	    k = j + i*jS + p*jS*iS + r*iS*jS*pS;
	    (*(fnew+k))=d_taus[2]*(*(f+k))+d_taus[1]*d_feq(vJp,eEp,bB,r);
	    k0 = r;
	    (*(f0new+k0))=d_taus[2]*(*(f0+k0))+d_taus[1]*rho;
	  }
}
__global__ void d_Adveccione(float * d_f, float * d_fnew, float * d_f0, float * d_f0new){
  int dist_f = 48;//iS*jS*pS+2*iS*jS+3*jS+1;
  int dist_f0 = 2; 
  int ix = blockIdx.x*blockDim.x+threadIdx.x;
  int iy = blockIdx.y*blockDim.y+threadIdx.y;
  int iz = blockIdx.z*blockDim.z+threadIdx.z;
  float *f, *fnew, *f0, *f0new;
  f = d_f + dist_f*(ix*Ly*Lz+iy*Lz+iz); f0 = d_f0 + (ix*Ly*Lz+iy*Lz+iz)*dist_f0;
  int jx, jy, jz, k;
  for(int r=0;r<2;r++)
    for(int p=0;p<3;p++)
      for(int i=0;i<4;i++)
	for(int j=0;j<2;j++){
	  if(ix==0 || ix==Lx-1){jx=ix;}
	  else{jx=ix+(int)d_vx[i+p*4];};
	  if(iy==0 || iy==Ly-1){jy=iy;}
	  else{jy=iy+(int)d_vy[i+p*4];};
	  if(iz==0 || iz==Lz-1){jz=iz;}
	  else{jz=iz+(int)d_vz[i+p*4];};
	  k = j + i*jS + p*jS*iS + r*iS*jS*pS;	  
	  fnew = d_fnew + dist_f*(jx*Ly*Lz+jy*Lz+jz); f0new = d_f0new + (jx*Ly*Lz+jy*Lz+jz)*dist_f0;
	  (*(f+k)) = (*(fnew+k));
	  (*(f0+r)) = (*(f0new+r));	       	
	}
}
class LatticeBoltzmann{
protected:
  float c_e_m[3]={0.0,0.0,0.0};
  float taus[3]={0.0,0.0,0.0};
  float V0[3]={0.0,0.0,0.0};
  float e0[3]={0.0,0.0,0.0};
  float b0[3]={0.0,0.0,0.0};
  
  float V[3][4][3];
  float vx[12]; float vy[12];	float vz[12];
  float ex[24];	float ey[24];	float ez[24];
  float bx[24];	float by[24];	float bz[24];
 
  float * f = nullptr; float * fnew = nullptr;
  float * f0 = nullptr; float * f0new = nullptr;

  float * d_f = nullptr; float * d_fnew = nullptr;
  float * d_f0 = nullptr; float * d_f0new = nullptr;

  //Distribution f= Distribution(Lx,Ly,Lz,false);		Distribution fnew = Distribution(Lx,Ly,Lz,false);
  //Distribution f0 = Distribution(Lx,Ly,Lz,true);	Distribution f0new = Distribution(Lx,Ly,Lz,true);
  float Eo=0.0001,Eop,lambda=16,T=lambda/C,omega=2*M_PI/T,k=omega/C;
  int ix0=Lx/2,iy0=Ly/2,iz0=Lz/2;
  
public:
  LatticeBoltzmann(void);
  ~LatticeBoltzmann(void);
  void ResizeDomain(int Lx0, int Ly0, int Lz0);
  float rho(int ix, int iy, int iz, bool useNew);
  //campos
  float Dx(int ix, int iy, int iz, bool useNew);
  float Dy(int ix, int iy, int iz, bool useNew);
  float Dz(int ix, int iy, int iz, bool useNew);
  float Bx(int ix, int iy, int iz, bool useNew);
  float By(int ix, int iy, int iz, bool useNew);
  float Bz(int ix, int iy, int iz, bool useNew);
  float Ex(float&Dx0, float&epsr);
  float Ey(float&Dy0, float&epsr);
  float Ez(float&Dz0, float&epsr);
  float Hx(float&Bx0, float&mur);
  float Hy(float&By0, float&mur);
  float Hz(float&Bz0, float&mur);
  float Jx(float&Ex0, float&sigma);
  float Jy(float&Ey0, float&sigma);
  float Jz(float&Ez0, float&sigma);
  //campos auxiliares
  float Jxp(float&Ex0, float&denominator);
  float Jyp(float&Ey0, float&denominator);
  float Jzp(float&Ez0, float&denominator);
  float Exp(float&Ex0, float&Jxp0, float&factor);
  float Eyp(float&Ey0, float&Jyp0, float&factor);
  float Ezp(float&Ez0, float&Jzp0, float&factor);
  //constantes dielectricas relativas
  float epsr(int ix, int iy, int iz){return 1.0;}; 
  float mur(int ix, int iy, int iz){return 1.0;};
  float sigma(int ix, int iy, int iz){return 0.0;};
  //funciones de equilibrio
  float feq(float&epsr,float&mur,float&vJp,float&eEp,float&bB,int r);
  float feq0(float&rho0);
  void Colisione(int &t);
  void h_Colisione(int&t);
  void Adveccione(void);
  void h_Adveccione(void);
  void Inicie(void);
  //void ImponerCampos(float&rho0,float*&D0,float*&B0,float*&H0,float*&E0,float*&J0,float*&Jp0,float*&Ep0,int t);
  void Imprimase(const char* fileName,bool useNew);
};
LatticeBoltzmann::LatticeBoltzmann(void){
  int imp=0,rec=0,rec1=0,rec2=0;
  float pi4=M_PI*0.25, sq2=sqrt(2.0);
  f = new float[sizef];	fnew = new float[sizef];
  f0 = new float[sizef0]; f0new = new float[sizef0];
  c_e_m[0] = C; c_e_m[1] = eps0; c_e_m[2] = mu0;
  taus[0] = tau; taus[1] = Utau; taus[2] = UmUtau;
   
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
  //vectores al device
  cudaMemcpyToSymbol(d_c_e_m,c_e_m,3*sizeof(float),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_taus,taus,3*sizeof(float),0,cudaMemcpyHostToDevice);
  
  cudaMemcpyToSymbol(d_vx,vx,12*sizeof(float),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_vy,vy,12*sizeof(float),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_vz,vz,12*sizeof(float),0,cudaMemcpyHostToDevice);

  cudaMemcpyToSymbol(d_ex,ex,24*sizeof(float),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_ey,ey,24*sizeof(float),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_ez,ez,24*sizeof(float),0,cudaMemcpyHostToDevice);
  
  cudaMemcpyToSymbol(d_bx,bx,24*sizeof(float),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_by,by,24*sizeof(float),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_bz,bz,24*sizeof(float),0,cudaMemcpyHostToDevice);
  //funciones al device
  cudaMalloc((void**) &d_f,sizef*sizeof(float));
  cudaMalloc((void**) &d_fnew,sizef*sizeof(float));
  cudaMalloc((void**) &d_f0,sizef0*sizeof(float));
  cudaMalloc((void**) &d_f0new,sizef0*sizeof(float));
  //variables dipolo
  Eo=0.0001;
  T=25;
  lambda=C*T;
  omega=2*M_PI/T;
  k=omega/C;
  
  ix0=Lx/2;
  iy0=Ly/2;
  iz0=Lz/2;
}
void LatticeBoltzmann::ResizeDomain(int Lx0, int Ly0, int Lz0)
{
  f = new float[sizef];	fnew = new float[sizef];
  f0 = new float[sizef0];	f0new = new float[sizef0];  
}
LatticeBoltzmann::~LatticeBoltzmann()
{
  delete f ;	delete fnew;
  delete f0;	delete f0new;
  cudaFree(d_f); cudaFree(d_fnew);
  cudaFree(d_f0); cudaFree(d_f0new);  
}
float LatticeBoltzmann::rho(int ix, int iy, int iz, bool useNew)
{
  float sum=0;int k=0;int k0=0;
  for(int p=0;p<3;p++)
    for(int i=0;i<4;i++)      
      for(int j=0;j<2;j++){
	k=j + i*jS + p*jS*iS + iz*iS*jS*rS*pS + iy*iS*jS*rS*pS*Lz + ix*iS*jS*rS*pS*Lz*Ly;
 	sum += useNew? (*(fnew+k)) : (*(f+k));}
  k0=iz*rS + iy*rS*Lz + ix*rS*Lz*Ly;
  sum += (*f0+k0);
  return sum;
}
float LatticeBoltzmann::Dx(int ix, int iy, int iz, bool useNew)
{
  float sum=0;int k=0;
  for(int p=0;p<3;p++)
    for(int i=0;i<4;i++)      
      for(int j=0;j<2;j++){
	k=j + i*jS + p*jS*iS + iz*iS*jS*rS*pS + iy*iS*jS*rS*pS*Lz + ix*iS*jS*rS*pS*Lz*Ly;
	sum += useNew? (*(fnew+k))*ex[j+i*2+p*8] : (*(f+k))*ex[j+i*2+p*8];}
  return sum;
}
float LatticeBoltzmann::Dy(int ix, int iy, int iz, bool useNew)
{
  float sum=0;int k=0;
  for(int p=0;p<3;p++)
    for(int i=0;i<4;i++)      
      for(int j=0;j<2;j++){
	k=j + i*jS + p*jS*iS + iz*iS*jS*rS*pS + iy*iS*jS*rS*pS*Lz + ix*iS*jS*rS*pS*Lz*Ly;
	sum += useNew? (*(fnew+k))*ey[j+i*2+p*8] : (*(f+k))*ey[j+i*2+p*8] ;}
  return sum;
}
float LatticeBoltzmann::Dz(int ix, int iy, int iz, bool useNew)
{
  float sum=0;int k=0;
  for(int p=0;p<3;p++)
    for(int i=0;i<4;i++)      
      for(int j=0;j<2;j++){
	k=j + i*jS + p*jS*iS + iz*iS*jS*rS*pS + iy*iS*jS*rS*pS*Lz + ix*iS*jS*rS*pS*Lz*Ly;
	sum += useNew? (*(fnew+k))*ez[j+i*2+p*8] : (*(f+k))*ez[j+i*2+p*8] ;}		
  return sum;
}
float LatticeBoltzmann::Bx(int ix, int iy, int iz, bool useNew)
{
  float sum=0;int k=0;
  for(int p=0;p<3;p++)
    for(int i=0;i<4;i++)      
      for(int j=0;j<2;j++){
	k=j + i*jS + p*jS*iS + iS*jS*pS + iz*iS*jS*rS*pS + iy*iS*jS*rS*pS*Lz + ix*iS*jS*rS*pS*Lz*Ly;
	sum += useNew? (*(fnew+k))*bx[j+i*2+p*8] : (*(f+k))*bx[j+i*2+p*8] ;}		
  return sum;
}
float LatticeBoltzmann::By(int ix, int iy, int iz, bool useNew)
{
  float sum=0;int k=0;
  for(int p=0;p<3;p++)
    for(int i=0;i<4;i++)      
      for(int j=0;j<2;j++){
	k=j + i*jS + p*jS*iS + iS*jS*pS + iz*iS*jS*rS*pS + iy*iS*jS*rS*pS*Lz + ix*iS*jS*rS*pS*Lz*Ly;
	sum += useNew? (*(fnew+k))*by[j+i*2+p*8] : (*(f+k))*by[j+i*2+p*8] ;}    
  return sum;
}
float LatticeBoltzmann::Bz(int ix, int iy, int iz, bool useNew)
{
  float sum=0;int k=0;
  for(int p=0;p<3;p++)
    for(int i=0;i<4;i++)      
      for(int j=0;j<2;j++){
	k=j + i*jS + p*jS*iS + iS*jS*pS + iz*iS*jS*rS*pS + iy*iS*jS*rS*pS*Lz + ix*iS*jS*rS*pS*Lz*Ly;
	sum += useNew? (*(fnew+k))*bz[j+i*2+p*8] : (*(f+k))*bz[j+i*2+p*8] ;}	
  return sum;
}
float LatticeBoltzmann::Ex(float&Dx0, float&epsr)
{
  float sum=0;
  
  sum = Dx0*(1.0/epsr);
  return sum;
}
float LatticeBoltzmann::Ey(float&Dy0, float&epsr)
{
  float sum=0;
  
  sum = Dy0*(1.0/epsr);
  return sum;
}
float LatticeBoltzmann::Ez(float&Dz0, float&epsr)
{
  float sum=0;
  
  sum = Dz0*(1.0/epsr);
  return sum;
}
float LatticeBoltzmann::Hx(float&Bx0, float&mur)
{
  float sum=0;

  sum = Bx0*(1.0/mur);
  return sum;
}
float LatticeBoltzmann::Hy(float&By0, float&mur)
{
  float sum=0;

  sum = By0*(1.0/mur);
  return sum;
}
float LatticeBoltzmann::Hz(float&Bz0, float&mur)
{
  float sum=0;

  sum = Bz0*(1.0/mur);
  return sum;
}
float LatticeBoltzmann::Jx(float&Ex0, float&sigma)
{
  float sum=0;

  sum = Ex0*sigma;
  return sum; 
}
float LatticeBoltzmann::Jy(float&Ey0, float&sigma)
{
  float sum=0;

  sum = Ey0*sigma;
  return sum; 
}
float LatticeBoltzmann::Jz(float&Ez0, float&sigma)
{
  float sum=0;

  sum = Ez0*sigma;
  return sum; 
}
float LatticeBoltzmann::Jxp(float&Ex0, float&denominator)
{
  float sum=0;

  sum = Ex0*denominator;
  return sum; 
}
float LatticeBoltzmann::Jyp(float&Ey0, float&denominator)
{
  float sum=0;

  sum = Ey0*denominator;
  return sum; 
}
float LatticeBoltzmann::Jzp(float&Ez0, float&denominator)
{
  float sum=0;

  sum = Ez0*denominator;
  return sum; 
}
float LatticeBoltzmann::Exp(float&Ex0,float&Jxp0,float&factor)
{
  float sum=0;

  sum = Ex0-Jxp0*factor;
  return sum;
}
float LatticeBoltzmann::Eyp(float&Ey0,float&Jyp0,float&factor)
{
  float sum=0;

  sum = Ey0-Jyp0*factor;
  return sum;
}
float LatticeBoltzmann::Ezp(float&Ez0,float&Jzp0,float&factor)
{
  float sum=0;

  sum = Ez0-Jzp0*factor;
  return sum;
}
float LatticeBoltzmann::feq(float&epsr, float&mur, float&vJp, float&eEp, float&bB,int r)
{
  float f=0;
  if(r==0)
    {f = (0.0625*vJp) + (epsr*0.25*eEp) + (0.125*bB/(mur));}
  else if(r==1)
    {f = (0.0625*vJp) + (0.25*eEp) + (0.125*bB);}
  
  return f;  
}
float LatticeBoltzmann::feq0(float&rho0)
{
  return rho0;
}
void LatticeBoltzmann::Colisione(int &t)
{
  int k=0,k0=0,ix=0,iy=0,iz=0,r=0,j=0,i=0,p=0;
  float Epsr=0,Mur=0,Sigma=0,denominator=0,factor=0;
  float rho0=0;
  float Dx0,Bx0,Ex0,Hx0,Jx0,Jxp0,Exp0;
  float Dy0,By0,Ey0,Hy0,Jy0,Jyp0,Eyp0;
  float Dz0,Bz0,Ez0,Hz0,Jz0,Jzp0,Ezp0;
  float vJp, eEp, bB;
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
	  Eop = Eo*exp(-0.25*((ix-ix0)*(ix-ix0)+(iy-iy0)*(iy-iy0)+(iz-iz0)*(iz-iz0)));
	  
	  Jzp0 = Eop*sin(omega*t);
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
void LatticeBoltzmann::h_Colisione(int&t){
  dim3 ThreadsPerBlock(8,8,8);
  dim3 BlocksPerGrid(16,16,16);
  d_Colisione<<<BlocksPerGrid,ThreadsPerBlock>>>(d_f, d_fnew, d_f0, d_f0new,t);
}
void LatticeBoltzmann::h_Adveccione(void){
  dim3 ThreadsPerBlock(8,8,8);
  dim3 BlocksPerGrid(16,16,16);
  d_Adveccione<<<BlocksPerGrid,ThreadsPerBlock>>>(d_f, d_fnew, d_f0, d_f0new);
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
  //float Eo=0.001,Bo=Eo/C,alp=0.01;
  //int iz0=40;
  int k=0,k0=0;
  float Epsr,Mur;
  float rho0;
  float vJp, eEp, bB;
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
  cudaMemcpy(d_f,f,sizef*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(d_f0,f0,sizef0*sizeof(float),cudaMemcpyHostToDevice);
  
}
void LatticeBoltzmann::Imprimase(const char* fileName,bool useNew)
{
  cudaMemcpy(f,d_f,sizef*sizeof(float),cudaMemcpyDeviceToHost);
  cudaMemcpy(f0,d_f0,sizef0*sizeof(float),cudaMemcpyDeviceToHost);
  cudaMemcpy(fnew,d_fnew,sizef*sizeof(float),cudaMemcpyDeviceToHost);
  cudaMemcpy(f0new,d_f0new,sizef0*sizeof(float),cudaMemcpyDeviceToHost);
  ofstream outputFile(fileName);
  float By0;
  float Dz0;
  float sigma=0;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int iz=0;iz<Lz;iz++)
	{
	  By0=By(ix,iy,iz,useNew);
	  Dz0=Dz(ix,iy,iz,useNew);
	  outputFile
	    << ix << "\t"
	    << iy << "\t"
	    << iz << "\t"
	    << Dz0 << "\t"
	    << By0 << "\n";	      
	}
   outputFile.close();
}

int main(int argc, char * argv[])
{
  if(argc != 3)
    {
      cout << "Wrong parameters!! ./a.out time";
      return -1;
    }
  const char* fileName = argv[2];
  LatticeBoltzmann Dip = LatticeBoltzmann();
  //Dip.ResizeDomain(Lx0,Ly0,Lz0);
  int t=0,tmax=atoi(argv[1]);
  
  
  Dip.Inicie();
  for(t=0;t<tmax;t++)
    {
      Dip.h_Colisione(t);
      Dip.h_Adveccione();
    }
  Dip.Imprimase(fileName,true);
  return 0;
}
