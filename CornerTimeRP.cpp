#include "headers/LB_ED_D3Q13.h"
#include <string>
int L=400;
int Lx0=L;
int Ly0=L;
int Lz0=1;
double size=1;

class Corner : public LatticeBoltzmann
{
private:
  double Eo=0.0001,Eop,lambda=16,T=lambda/C,omega=2*M_PI/T,k=omega/C;
  double sigma0=10,thickness=0.8,order=0.5,delta=1;
  int offset=1;
  int ix0=Lx0/2,iy0=Ly0/2,iz0=Lz0/2,n=100;

  double Eref=0,sheet=0;
  vector<double> Emax;
  vector<double> Emin;
public:
  Corner(double order0);
  double GetT(){return T;};
  double GetOrder(){return order;};
  double SigmaReflectors(int ix, int iy, int iz);
  void ColisioneCorner(int&t);
  void FreeBoundAdvection(void);
  void InicieCorner(void);
  void PerfectConductor(void);
  void RadiationPattern(int&t,double r);
  void ImprimirPattern(const char* fileName,double r,bool useNew);
  void ImprimirCorner(const char* fileName,bool useNew,int&t);
  vector3D Interpolation(int ix1, int ix2,int iy1,int iy2,double x,double y,bool getB);
  vector3D GetSField(int ix,int iy, int iz);
  vector3D S(vector3D&E0, vector3D&H0);
};
Corner::Corner(double order0)
{
  Eo=0.0001;
  lambda=16;
  T=lambda/C;
  omega=2*M_PI/T;
  k=omega/C;
  
  thickness=1;
  offset=1;
  sigma0=10;//lambda/(delta*delta*M_PI*mu0*C);
  order=order0*0.5;

  ix0=(int)(lambda*order/sqrt(2))+offset;
  iy0=(int)(lambda*order/sqrt(2))+offset;
  iz0=Lz0/2;

  Eref=2.091104e-5;
  n=100;
  Emax.resize(n);
  Emin.resize(n);
}
double Corner::SigmaReflectors(int ix, int iy, int iz)
{
  //en Y=0
  double Yplate = 0.5*(1-tanh(thickness*(-ix+offset)));
  //en X=0
  double Xplate = 0.5*(1-tanh(thickness*(-iy+offset)));
  //bothplates
  double plates = sigma0*(1-Yplate*Xplate);

  return plates;
}
vector3D Corner::S(vector3D&E0,vector3D&H0)
{
  vector3D s = E0^H0;
  return s;
}
vector3D Corner::GetSField(int ix,int iy, int iz)
{
  vector3D D0, B0, E0, H0, S0;
  double Epsr=1.0,Mur=1.0;
  D0=D(ix,iy,iz0,false);	E0=E(D0,Epsr);
  B0=B(ix,iy,iz0,false);	H0=H(B0,Mur);
  S0=S(E0,H0);
  return S0;
}
void Corner::ColisioneCorner(int &t)
{
  int ix=0,iy=0,iz=0,r=0,j=0,i=0,p=0;
  double Epsr=0,Mur=0,Sigma=0,denominator=0,factor=0;
  double rho0=0;
  vector3D D0,B0,E0,H0,J0,Jp0,Ep0,S0;
  
  
  int l4=(int)(lambda*0.25);
  
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++)
	{
	  Epsr=epsr(ix,iy,iz);		Mur=mur(ix,iy,iz);	     Sigma=SigmaReflectors(ix,iy,iz);
	  denominator = Sigma/(1.0 + mu0*Sigma/(4.0*Epsr));
	  factor = mu0/(4.0*Epsr);
	  rho0=rho(ix,iy,iz,false);	 D0=D(ix,iy,iz,false);	     B0=B(ix,iy,iz,false);
	  H0=H(B0,Mur);	 E0=E(D0,Epsr);		S0=S(E0,H0);
	  J0=J(E0,Sigma);	 //Jp0=Jp(E0,ix,iy,iz,false);

	  Eop = Eo*exp(-0.25*((ix-ix0)*(ix-ix0)+(iy-iy0)*(iy-iy0)));
	  
	  if(iz>=iz0-l4 && iz<=iz0+l4)
	    Eop *= cos(k*(iz-iz0));
	  else
	    Eop *= 0;
	  
	   if(Sigma > 0.24*sigma0)
	     {//conductor
	      Jp0=Jp(E0,denominator);
	    }
	   else//dipolo
	     {
	       Jp0[0] = 0; Jp0[1] = 0;
	       Jp0[2] = Eop*sin(omega*t);// + Jp(E0,denominator).z();
	     }
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
void Corner::FreeBoundAdvection(void)
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
		  //jx=(ix+(int)V[0][i][p]+Lx)%Lx;jy=(iy+(int)V[1][i][p]+Ly)%Ly;jz=(iz+(int)V[2][i][p]+Lz)%Lz;
		  f(r,j,p,i,jx,jy,jz) = fnew.function(r,j,p,i,ix,iy,iz);
		  f0(r,ix,iy,iz) = f0new.function(r,ix,iy,iz);
		}
}
void Corner::InicieCorner(void)
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
void Corner::PerfectConductor(void)
{
  double Epsr,Mur,Sigma,denominator,factor;
  double rho0; vector3D D0,B0,E0,H0,J0,Jp0,Ep0;
  int iyp=offset;
  for(iyp=0;iyp<=offset;iyp++){
  for(int ix=offset;ix<Lx;ix++)
    {
      Epsr=epsr(ix,iyp,0);
      Mur=mur(ix,iyp,0);
      Sigma=sigma(ix,iyp,0);
      denominator = Sigma/(1.0 + mu0*Sigma/(4.0*Epsr));
      factor = mu0/(4.0*Epsr);
      
      D0[0]=0;	D0[2]=0;	D0[1]=D(ix,iyp,0,false).y();
      B0[1]=0;	B0[0]=B(ix,iyp,0,false).x();
      B0[2]=B(ix,iyp,0,false).z();
      
      rho0=rho(ix,iyp,0,false);
      E0=E(D0,Epsr);
      H0=H(B0,Mur);
      J0=J(E0,Sigma);
      Jp0=Jp(J0,denominator);
      Ep0=Ep(E0,Jp0,factor);
      
      for(int p=0;p<3;p++)
	for(int i=0;i<4;i++)
	  for(int j=0;j<2;j++)
	    for(int r=0;r<2;r++)
	      {
		fnew(r,j,p,i,ix,iyp,0) = feq(Epsr,Mur,Jp0,Ep0,B0,r,j,i,p);
		f0new(r,ix,iyp,0) = feq0(rho0);		
	      }
    }
  for(int iy=offset;iy<Ly;iy++)
    {
      Epsr=epsr(iyp,iy,0);
      Mur=mur(iyp,iy,0);
      Sigma=sigma(iyp,iy,0);
      denominator = Sigma/(1.0 + mu0*Sigma/(4.0*Epsr));
      factor = mu0/(4.0*Epsr);
      
      D0[1]=0;	D0[2]=0;	D0[0]=D(iyp,iy,0,false).x();
      B0[0]=0;	B0[1]=B(iyp,iy,0,false).y();
      B0[2]=B(iyp,iy,0,false).z();
      
      rho0=rho(iyp,iy,0,false);
      E0=E(D0,Epsr);
      H0=H(B0,Mur);
      J0=J(E0,Sigma);
      Jp0=Jp(J0,denominator);
      Ep0=Ep(E0,Jp0,factor);
      
      for(int p=0;p<3;p++)
	for(int i=0;i<4;i++)
	  for(int j=0;j<2;j++)
	    for(int r=0;r<2;r++)
	      {
		fnew(r,j,p,i,iyp,iy,0) = feq(Epsr,Mur,Jp0,Ep0,B0,r,j,i,p);
		f0new(r,iyp,iy,0) = feq0(rho0);		
	      }
    }}
}
vector3D Corner::Interpolation(int ix1, int ix2,int iy1,int iy2,double x,double y,bool getB)
{
  double u = (x-ix1)/(ix2-ix1);
  double v = (y-iy1)/(iy2-iy1);
  vector3D F0;

  double
    f11=(1-u)*(1-v),
    f12=(1-u)*v,
    f21=(1-v)*u,
    f22=u*v;
  vector3D F11,F12,F21,F22;
  if(getB){
    F11=B(ix1,iy1,iz0,false),
    F12=B(ix1,iy2,iz0,false),
    F21=B(ix2,iy1,iz0,false),
    F22=B(ix2,iy2,iz0,false);
  }else{
    F11=D(ix1,iy1,iz0,false),
    F12=D(ix1,iy2,iz0,false),
    F21=D(ix2,iy1,iz0,false),
    F22=D(ix2,iy2,iz0,false);
  }
  F0=f11*F11+f12*F12+f21*F21+f22*F22;
  return F0;
}
void Corner::RadiationPattern(int&t,double r)
{
  vector3D E0, B0, S0;
  double x=0,y=0;
  int ix1=0, iy1=0, ix2=0, iy2=0,ia=0;
  double angle=0,deltaAng=M_PI*0.5/n,Exy=0;
  double r2=(r)*lambda;
  Emax.resize(n);
  for(ia=0; ia<n;ia++)
    {
      angle = ia*deltaAng;
      x = r2*cos(angle)+offset;	ix1=floor(x);	ix2=ix1+1;
      y = r2*sin(angle)+offset;	iy1=floor(y);	iy2=iy1+1;
      E0=Interpolation(ix1,ix2,iy1,iy2,x,y,false);
      B0=Interpolation(ix1,ix2,iy1,iy2,x,y,true);
      S0=S(E0,B0);
      Exy=norma(S0);
      if(t=0){
	Emax.at(ia)=Exy;
	Emin.at(ia)=Exy;
      }
      if(Exy>Emax.at(ia))
	Emax.at(ia)=Exy;
      if(Exy<Emin.at(ia))
	Emin.at(ia)=Exy;	      
    }
}
void Corner::ImprimirCorner(const char* fileName,bool useNew,int&t)
{
  ofstream outputFile(fileName);
  vector3D D0, B0, S0, Jp0;
  double Sigma=0,denominator=0;
  //double Epsr=1.0,Mur=1.0,E2=0;
 
  for(double ix=0;ix<Lx0;ix+=1.0)
    for(double iy=0;iy<Ly0;iy+=1.0)
      {       
	Sigma = SigmaReflectors(ix,iy,iz0);
	denominator = Sigma/(1.0 + mu0*Sigma/4.0);
	D0=D(ix,iy,iz0,false);//	E0=E(D0,Epsr);
	B0=B(ix,iy,iz0,false);//	H0=H(B0,Mur);
	S0=S(D0,B0);
	if(Sigma > 0.24*sigma0)
	  {//conductor
	    Jp0=Jp(D0,denominator);
	    sheet=1;
	    }
	else//dipolo
	  {
	    Jp0[0] = 0; Jp0[1] = 0;
	    Jp0[2] = Eop*sin(omega*t);// + Jp(E0,denominator).z();
	    sheet=0;
	  }	
	outputFile
	<< ix << "\t"
	<< iy << "\t"
	<< Sigma << "\t"
	<< sheet << "\t"
	<< norma(Jp0) << "\t"
	<< norma(S0) << "\n";
      }
  outputFile.close();
}
void Corner::ImprimirPattern(const char* fileName,double r,bool useNew)
{
  ofstream outputFile(fileName);
  int ia=0;
  double x=0,y=0,angle=0,deltaAng=M_PI*0.5/n,Eteorico=0,r2=(r)*lambda;
  for(ia=0;ia<n;ia++)
    {
      angle=ia*deltaAng;
      x = r2*cos(angle)+offset;
      y = r2*sin(angle)+offset;
      Eteorico = 2*(cos(order*2*M_PI*cos(angle+M_PI*.25))-cos(order*2*M_PI*sin(angle+M_PI*0.25)));
      outputFile
	<< x << "\t"
	<< y << "\t"
	<< angle << "\t"	      
	<< Emax.at(ia) << "\t"	      
	<< Emin.at(ia) << "\t"	      
	<< Eteorico << "\n";	      
    }
  outputFile.close();
}
int main(int argc, char * argv[])
{
  if(argc != 6)
    {
      cout << "Wrong parameters!! ./a.out time Sfile patternfile animate(0 or 1) 2order";
      return -1;
      }
  const char* fileName = argv[2];
  const char* patternFile = argv[3];
  const char* frameFile = ".dat";
  double orderN = atoi(argv[5]);
  Corner corner = Corner(orderN);
  corner.ResizeDomain(Lx0,Ly0,Lz0);
  double T0 = corner.GetT(), order0 = corner.GetOrder(),NN = atoi(argv[1]);
  int t=0,tmax=NN*T0,T3=(NN-1)*T0,tp=0,frameSkip = 2;
  double r0 = NN/2;
  bool anim = argv[4];
    
  corner.InicieCorner();
  for(t=0;t<tmax;t++)
    {
      corner.ColisioneCorner(t);
      //corner.PerfectConductor();
      corner.FreeBoundAdvection();
      /*if(anim && t%frameSkip==0)
	{
	  frameFile = to_string(t).c_str();
	  corner.ImprimirCorner(frameFile,true);
	  }*/
      if(t>=T3)
	{
	  tp=t-T3;
	  corner.RadiationPattern(tp,r0);
	}
    }
  corner.ImprimirPattern(patternFile,r0,true);
  corner.ImprimirCorner(fileName,true,tmax);
  return 0;
}
