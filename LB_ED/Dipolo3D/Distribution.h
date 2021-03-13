#include <iostream>
#include <vector>

const int rS = 2;
const int jS = 2;
const int pS = 3;
const int iS = 4;

/*
lattice:
int Lx
int Ly
int Lz
*/

using namespace std;
//Class definition for distribution
class Distribution
{
 public:
  int size=0;
  bool zero = false;
  double Lx=0,Ly=0,Lz=0;
  vector<double> array;
 
  //constructor
  Distribution(int Lx0, int Ly0, int Lz0, bool zero0);
  ~Distribution(void);
  //acces function element
  double function(int r, int j, int p, int i, int ix, int iy, int iz);
  double function(int r, int ix, int iy, int iz);
  //set function element
  double & operator()(int r, int j, int p, int i, int ix, int iy, int iz);
 
  double & operator()(int r, int ix, int iy, int iz);
  

  
 private:
};
//Implementating functions
Distribution::Distribution(int Lx0, int Ly0, int Lz0, bool zero0)
{
  zero = zero0;
  Lx=Lx0; Ly=Ly0; Lz=Lz0;
  size = zero? rS*Lx*Ly*Lz : rS*jS*pS*iS*Lx*Ly*Lz;
  array.resize(size);
}

Distribution::~Distribution(void)
{
  array.clear();
}

double Distribution::function(int r, int j, int p, int i, int ix, int iy, int iz)
{  
  double function = array.at(ix + iy*Lx + iz*Lx*Ly + i*Lx*Ly*Lz + p*Lx*Ly*Lz*iS + j*Lx*Ly*Lz*iS*pS + r*Lx*Ly*Lz*iS*pS*jS);

  return function;
}

double Distribution::function(int r, int ix, int iy, int iz)
{  
  double function = array.at(ix + iy*Lx + iz*Lx*Ly + r*Ly*Lx*Lz);

  return function;
}

double & Distribution::operator()(int r, int j, int p, int i, int ix, int iy, int iz)
{
  return array.at(ix + iy*Lx + iz*Lx*Ly + i*Lx*Ly*Lz + p*Lx*Ly*Lz*iS + j*Lx*Ly*Lz*iS*pS + r*Lx*Ly*Lz*iS*pS*jS);
}

double & Distribution::operator()(int r, int ix, int iy, int iz)
{
  return array.at(ix + iy*Lx + iz*Lx*Ly + r*Ly*Lx*Lz);   
}
