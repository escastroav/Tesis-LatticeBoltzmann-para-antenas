#include "headers/LB_ED_D3Q13.h"

int Lx0=100;
int Ly0=100;
int Lz0=100;

class Dipolo : public LatticeBoltzmann
{
private:
  double Eo=0.0001,Eop,lambda=16,T=lambda/C,omega=2*M_PI/T,k=omega/C;
  int ix0=Lx0/2,iy0=Ly0/2,iz0=Lz0/2;

public:
  Dipolo(void);
  void ColisioneDipolo(int&t);
  void FreeBoundAdvection(void);
  void InicieDipolo(void);
  void ImprimirDipolo(const char* fileName,bool useNew);
};
Dipolo::Dipolo(void)
{
  Eo=0.0001;
  T=25;
  lambda=C*T;
  omega=2*M_PI/T;
  k=omega/C;
  
  ix0=Lx0/2;
  iy0=Ly0/2;
  iz0=Lz0/2;
}
void Dipolo::ColisioneDipolo(int &t)
{
  int ix=0,iy=0,iz=0,r=0,j=0,i=0,p=0;
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
		    fnew(r,j,p,i,ix,iy,iz)=UmUtau*f.function(r,j,p,i,ix,iy,iz)+Utau*feq(Epsr,Mur,vJp,eEp,bB,r);
		    f0new(r,ix,iy,iz)=UmUtau*f0.function(r,ix,iy,iz)+Utau*feq0(rho0);
		  }
	}
}
void Dipolo::FreeBoundAdvection(void)
{
  int jx=0,jy=0,jz=0;
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
		  	       
		  f(r,j,p,i,jx,jy,jz) = fnew.function(r,j,p,i,ix,iy,iz);
		  f0(r,ix,iy,iz) = f0new.function(r,ix,iy,iz);
		}
}
void Dipolo::InicieDipolo(void)
{
  //double Eo=0.001,Bo=Eo/C,alp=0.01;
  //int iz0=40;
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
		      f(r,j,p,i,ix,iy,iz) = feq(Epsr,Mur,vJp,eEp,bB,r);
		      f0(r,ix,iy,iz) = feq0(rho0);		      
		    }
	  }
}
void Dipolo::ImprimirDipolo(const char* fileName,bool useNew)
{
  ofstream outputFile(fileName);
  double By0;
  double sigma=0;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int iz=0;iz<Lz;iz++)
	{
	  By0=By(ix,iy,iz,useNew);
	  outputFile
	    << ix << "\t"
	    << iy << "\t"
	    << iz << "\t"
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
  Dipolo Dip = Dipolo();
  Dip.ResizeDomain(Lx0,Ly0,Lz0);
  int t=0,tmax=atoi(argv[1]);
  
  
  Dip.InicieDipolo();
  for(t=0;t<tmax;t++)
    {
      Dip.ColisioneDipolo(t);
      Dip.FreeBoundAdvection();
    }
  Dip.ImprimirDipolo(fileName,true);
  return 0;
}
