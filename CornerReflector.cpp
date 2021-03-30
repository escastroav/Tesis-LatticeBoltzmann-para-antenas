#include "headers/LB_ED_D3Q13.h"

class Dipolo : public LatticeBoltzmann
{
private:
  double Eo=0.0001,Eop,lambda=16,T=lambda/C,omega=2*M_PI/T,k=omega/C;
  double sigma0=10,thickness=0.8;
  int offset=40;
  int ix0=Lx/2,iy0=Ly/2,iz0=Lz/2;

  double Emean=0;
public:
  Dipolo(void);
  double SigmaReflectors(int ix, int iy, int iz);
  void ColisioneDipolo(int&t);
  void FreeBoundAdvection(void);
  void InicieDipolo(void);
  void ImprimirDipolo(const char* fileName,bool useNew);
  double EtimeMean(vector3D&E0,int t);
};
Dipolo::Dipolo(void)
{
  Eo=0.0001;
  lambda=16;
  T=lambda/C;
  omega=2*M_PI/T;
  k=omega/C;
  
  thickness=0.75;
  offset=3;
  sigma0=lambda/(offset*offset*M_PI*mu0*C);

  ix0=Lx/2;
  iy0=Ly/2;
  iz0=Lz/2;
}
double Dipolo::SigmaReflectors(int ix, int iy, int iz)
{
  //en Y=0
  double Yplate = 0.5*(1-tanh(thickness*(-ix+offset)));
  //en X=0
  double Xplate = 0.5*(1-tanh(thickness*(-iy+offset)));
  //bothplates
  double plates = sigma0*(1-Yplate*Xplate);

  return plates;
}
double Dipolo::EtimeMean(vector3D&E0,int t)
{
  return norma(E0)*t/T;
}
void Dipolo::ColisioneDipolo(int &t)
{
  int ix=0,iy=0,iz=0,r=0,j=0,i=0,p=0;
  double Epsr=0,Mur=0,Sigma=0,denominator=0,factor=0;
  double rho0=0;
  vector3D D0,B0,E0,H0,J0,Jp0,Ep0;
  
  double Eo=0.0001,Eop,lambda=16,T=lambda/C,omega=2*M_PI/T,k=omega/C;
  int ix0=(int)(lambda*0.5/sqrt(2))+offset,iy0=(int)(lambda*0.5/sqrt(2))+offset,iz0=Lz/2,l4=(int)(lambda*0.25);
  
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++)
	{
	  Epsr=epsr(ix,iy,iz);		Mur=mur(ix,iy,iz);	     Sigma=SigmaReflectors(ix,iy,iz);
	  denominator = Sigma/(1.0 + mu0*Sigma/(4.0*Epsr));
	  factor = mu0/(4.0*Epsr);
	  rho0=rho(ix,iy,iz,false);	 D0=D(ix,iy,iz,false);	     B0=B(ix,iy,iz,false);
	  H0=H(B0,Mur);	 E0=E(D0,Epsr);
	  J0=J(E0,Sigma);	 //Jp0=Jp(E0,ix,iy,iz,false);

	  Eop = Eo*exp(-0.25*((ix-ix0)*(ix-ix0)+(iy-iy0)*(iy-iy0)));
	  
	  if(iz>=iz0-l4 && iz<=iz0+l4)
	    Eop *= cos(k*(iz-iz0));
	  else
	    Eop *= exp(-0.75*(iz-iz0)*(iz-iz0));
	  
	  Jp0[2] = Eop*sin(omega*t) + Jp(E0,denominator).z();
	  Ep0=Ep(E0,Jp0,factor);
	  
	  Emean+=EtimeMean(E0,t);
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
void Dipolo::FreeBoundAdvection(void)
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
void Dipolo::InicieDipolo(void)
{
  //double Eo=0.001,Bo=Eo/C,alp=0.01;
  //int iz0=40;
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
void Dipolo::ImprimirDipolo(const char* fileName,bool useNew)
{
  ofstream outputFile(fileName);
  vector3D D0, B0;					     
  int ix0 = Lx/2,iy0=Ly/2,iz0 = Lz/2;
  double sigma=0;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
    for(int iz=0;iz<Lz;iz++)
    {
      D0=D(ix,iy,iz,useNew);	
      B0=B(ix,iy,iz,useNew);
      sigma=SigmaReflectors(ix,iy,iz);
      outputFile
	<< ix << "\t"
	<< iy << "\t"
	<< iz << "\t"
	<< D0.z() << "\t"
	<< sigma << "\t"
	<< Emean << "\n";	      
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
