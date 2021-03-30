#include "headers/LB_ED_D3Q13.h"

int Lx0=1;
int Ly0=1;
int Lz0=200;

class Interfase : public LatticeBoltzmann
{
private:
  double Eo=0.001,Bo=Eo/C,alpha=0.01;
public:
  Interfase(void);
  double EpsrInt(int iz);
  void ColisioneInterfase(void);
  void PeriodicAdvection(void);
  void InicieInterfase(void);
  void ImprimirInt(const char* fileName,bool useNew);
};
Interfase::Interfase(void)
{
  Eo=0.001;
  Bo=Eo/C;
  alpha=0.01;
}
double Interfase::EpsrInt(int iz)
{
  return (1.75+0.75*tanh((double)(iz-Lz/2)));
}
void Interfase::ColisioneInterfase(void)
{
  int ix=0,iy=0,iz=0,r=0,j=0,i=0,p=0;
  double Epsr=0,Mur=0,Sigma=0,denominator=0,factor=0;
  double rho0=0;
  vector3D D0,B0,E0,H0,J0,Jp0,Ep0;
  
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++)
	{
	  Epsr=EpsrInt(iz);		Mur=mur(ix,iy,iz);	     Sigma=sigma(ix,iy,iz);
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
void Interfase::PeriodicAdvection(void)
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
void Interfase::InicieInterfase(void)
{
  double Epsr,Mur,Sigma,denominator,factor;
  double rho0; vector3D D0,B0,E0,H0,J0,Jp0,Ep0;
  int iz0 = 40;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int iz=0;iz<Lz;iz++)
	{
	  Epsr=EpsrInt(iz);
	  Mur=mur(ix,iy,iz);
	  Sigma=sigma(ix,iy,iz);
	  denominator = Sigma/(1.0 + mu0*Sigma/(4.0*Epsr));
	  factor = mu0/(4.0*Epsr);
	  
	  rho0=rho(ix,iy,iz,false);
	  D0[0] =Eo*Epsr*exp(-alpha*(iz-iz0)*(iz-iz0));	
	  B0[1] =Bo*exp(-alpha*(iz-iz0)*(iz-iz0));
	  H0=H(B0,Mur);	 E0=E(D0,Epsr);
	  J0=J(E0,Sigma);	 Jp0=Jp(E0,denominator);
	  Ep0=Ep(E0,Jp0,factor);
	  
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
void Interfase::ImprimirInt(const char* fileName,bool useNew)
{
  ofstream outputFile(fileName);
  vector3D D0,E0;					     
  double Epsr=0;
  for(int iz=0;iz<Lz;iz++)
    {
      D0=D(0,0,iz,useNew);		     
      Epsr = EpsrInt(iz);
      E0=E(D0,Epsr);
      outputFile
	<< iz << "\t"	
	<< E0.x()/Eo << "\t"	
	<< Epsr << "\n";	      
    }
   outputFile.close();
}
int main(int argc, char * argv[])
{
  if(argc != 4)
    {
      cout << "Wrong parameters!! ./a.out time filet0 filet140";
      return -1;
    }
  const char* fileName = argv[3];
  Interfase inter = Interfase();
  inter.ResizeDomain(Lx0,Ly0,Lz0);
  int t=0,tmax=atoi(argv[1]);
  
  
  inter.InicieInterfase();
  inter.ImprimirInt(argv[2],false);
  for(t=0;t<tmax;t++)
    {
      inter.ColisioneInterfase();
      inter.PeriodicAdvection();
    }
  inter.ImprimirInt(fileName,true);
  return 0;
}
