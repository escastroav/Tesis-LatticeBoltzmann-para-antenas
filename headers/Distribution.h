/*
Distribution.h:
Encargada de definir las funciones de distribucion en un arreglo tipo <vector> unidimensional.
El indice va de menor indice a mayor indice (r,j,p,i,Lx,Ly,Lz) y se accede al índice mediante:
ix + iy*Lx + iz*Lx*Ly + i*Lx*Ly*Lz + p*Lx*Ly*Lz*4 + j*Lx*Ly*Lz*4*3 + r*Lx*Ly*Lz*4*3*2
*/
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
  int k = j + i*jS + p*jS*iS + r*iS*jS*pS + iz*iS*jS*rS*pS + iy*iS*jS*rS*pS*Lz + ix*iS*jS*rS*pS*Lz*Ly;
  double function = array.at(k);

  return function;
}

double Distribution::function(int r, int ix, int iy, int iz)
{
  int k = r + iz*rS + iy*rS*Lz + ix*rS*Lz*Ly;
  double function = array.at(k);

  return function;
}

double & Distribution::operator()(int r, int j, int p, int i, int ix, int iy, int iz)
{
  int k = j + i*jS + p*jS*iS + r*iS*jS*pS + iz*iS*jS*rS*pS + iy*iS*jS*rS*pS*Lz + ix*iS*jS*rS*pS*Lz*Ly;
  return array.at(k);
}

double & Distribution::operator()(int r, int ix, int iy, int iz)
{
  int k = r + iz*rS + iy*rS*Lz + ix*rS*Lz*Ly;
  return array.at(k);   
}
