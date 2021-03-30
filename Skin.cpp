#include "headers/LB_ED_D3Q13.h"

int Lx0=1;
int Ly0=1;
int Lz0=1000;

class SkinWave : public LatticeBoltzmann
{
private:
  double Eo=0.001,Bo=Eo/C,omega=M_PI*0.02;
  double sigma0=0.1; 
public:
  SkinWave(void);
  double SigmaSkin(int iz);
  void ColisioneSkin(void);
  void ForcedWave(int&t);
  void PeriodicAdvection(void);
  void InicieSkin(void);
  void ImprimirSkin(const char* fileName,bool useNew);
};
SkinWave::SkinWave(void)
{
  Eo=0.001;
  Bo=Eo/C;
  omega=M_PI*0.01;
  
}
double SkinWave::SigmaSkin(int iz)
{
  return sigma0*0.5*(1+tanh((iz-Lz/4)));
}
void SkinWave::ColisioneSkin(void)
{
  int ix=0,iy=0,iz=0,r=0,j=0,i=0,p=0;
  double Epsr=0,Mur=0,Sigma=0,denominator=0,factor=0;
  double rho0=0;
  vector3D D0,B0,E0,H0,J0,Jp0,Ep0;
  
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++)
	{
	  Epsr=epsr(ix,iy,iz);		Mur=mur(ix,iy,iz);	     Sigma=SigmaSkin(iz);
	  denominator = Sigma/(1.0 + mu0*Sigma/(4.0*Epsr));
	  factor = mu0/(4.0*Epsr);
	  rho0=rho(ix,iy,iz,false);	 D0=D(ix,iy,iz,false);	     B0=B(ix,iy,iz,false);
	  H0=H(B0,Mur);	 E0=E(D0,Epsr);
	  J0=J(E0,Sigma);	 Jp0=Jp(E0,denominator);
	  Ep0=Ep(E0,Jp0,factor);
	  
	  for(p=0;p<3;p++)
	    for(i=0;i<4;i++)
	      for(j=0;j<2;j++)
		for(r=0;r<2;r++)
		{ 
		  //fnew[ix][iy][iz][r][j][i][p]=UmUtau*f[ix][iy][iz][r][j][i][p]+Utau*feq(ix,iy,iz,Jp0,Ep0,B0,r,j,i,p);
		  //f0new[ix][iy][iz][r]=UmUtau*f0[ix][iy][iz][r]+Utau*feq0(ix,iy,iz,rho0,r);
		  fnew(r,j,p,i,ix,iy,iz)=UmUtau*f.function(r,j,p,i,ix,iy,iz)+Utau*feq(Epsr,Mur,Jp0,Ep0,B0,r,j,i,p);
		  f0new(r,ix,iy,iz)=UmUtau*f0.function(r,ix,iy,iz)+Utau*feq0(rho0);
		}
	}
}
void SkinWave::ForcedWave(int&t)
{
  int ix=0,iy=0,iz=0;
  vector3D D0,B0,J0,E0,H0,Jp0,Ep0;
  double Epsr=epsr(ix,iy,iz),Mur=mur(ix,iy,iz),Sigma=SigmaSkin(iz),rho0=0;
  double denominator=0,factor=0;
  denominator = Sigma/(1.0 + mu0*Sigma/(4.0*Epsr));
  factor = mu0/(4.0*Epsr);
  D0[1]=0;	D0[2]=0;	B0[0]=0;	B0[2]=0;
  D0[0] = Eo*sin(omega*t);	
  B0[1] = Bo*sin(omega*t);
  rho0=rho(ix,iy,iz,false);
  E0=E(D0,Epsr);	H0=H(B0,Mur);
  J0=J(E0,Sigma);	
  Jp0=Jp(E0,denominator);	Ep0=Ep(E0,Jp0,factor);
  
  for(int p=0;p<3;p++)
    for(int i=0;i<4;i++)
      for(int j=0;j<2;j++)
	for(int r=0;r<2;r++)
	  {
	    fnew(r,j,p,i,ix,iy,iz) = feq(Epsr,Mur,Jp0,Ep0,B0,r,j,i,p);
	    f0new(r,ix,iy,iz) = feq0(rho0);	    
	  }
  
  
}
void SkinWave::PeriodicAdvection(void)
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
		  f(r,j,p,i,jx,jy,jz) = fnew.function(r,j,p,i,ix,iy,iz);
		  f0(r,ix,iy,iz) = f0new.function(r,ix,iy,iz);
		}
}
void SkinWave::InicieSkin(void)
{
  double Epsr,Mur,Sigma;
  double rho0; vector3D D0,B0,E0,H0,J0,Jp0,Ep0;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int iz=0;iz<Lz;iz++)
	  {
	    Epsr=epsr(ix,iy,iz);
	    Mur=mur(ix,iy,iz);
	    Sigma=SigmaSkin(iz);
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
void SkinWave::ImprimirSkin(const char* fileName,bool useNew)
{
  ofstream outputFile(fileName);
  vector3D D0;					     
  double sigma=0;
  for(int iz=0;iz<Lz;iz++)
    {
      D0=D(0,0,iz,useNew);	     
      sigma = SigmaSkin(iz);
      
      outputFile
	<< iz << "\t"	
	<< D0.x() << "\t"	
	<< sigma << "\n";	      
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
  SkinWave skin = SkinWave();
  skin.ResizeDomain(Lx0,Ly0,Lz0);
  int t=0,tmax=atoi(argv[1]);
  
  
  skin.InicieSkin();
  for(t=0;t<tmax;t++)
    {
      skin.ColisioneSkin();
      skin.ForcedWave(t);
      skin.PeriodicAdvection();
    }
  skin.ImprimirSkin(fileName,true);
  return 0;
}
