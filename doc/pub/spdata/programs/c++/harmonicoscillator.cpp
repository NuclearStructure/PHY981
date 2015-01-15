/*
  Solves the one-particle Schrodinger equation
  for a potential specified in function
  potential(). This example is for the harmonic oscillator in 3d
*/
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>

using namespace  std;
using namespace  arma;
// input and output files as global variables
ofstream ofile;  
ifstream ifile;

double potential(double);
void output(double, double, int, vec& );


// Begin of main program   

int main(int argc, char* argv[])
{
  char *outfilename, *infilename;
  int       i, j, Dim, lOrbital;
  double    RMin, RMax, Step, DiagConst, NondiagConst, OrbitalFactor; 
  // With spherical coordinates RMin = 0 always
  RMin = 0.0;
  if (argc <= 2) {
    cout << "Bad Usage: " << argv[0] << 
      " read also input and output files on same line" << endl;
    exit(1);
  }
  if (argc > 2) {
    infilename=argv[1];
    outfilename=argv[2];
    ifile.open(infilename); 
    ofile.open(outfilename); 
  }
  //  Read in data from file
  ifile >>  RMax;  ifile >> lOrbital;  ifile >> Dim;  
  mat Hamiltonian = zeros<mat>(Dim,Dim);
  // Integration step length
  Step    = RMax/ Dim;
  DiagConst = 2.0 / (Step*Step);
  NondiagConst =  -1.0 / (Step*Step);
  OrbitalFactor = lOrbital * (lOrbital + 1.0);
  
  // local memory for r and the potential w[r] 
  vec r(Dim); vec w(Dim);
  for(i = 0; i < Dim; i++) {
    r(i) = RMin + (i+1) * Step;
    w(i) = potential(r(i)) + OrbitalFactor/(r(i) * r(i));
  }


  // Setting up tridiagonal matrix and brute diagonalization using Armadillo
  Hamiltonian(0,0) = DiagConst + w(0);
  Hamiltonian(0,1) = NondiagConst;
  for(i = 1; i < Dim-1; i++) {
    Hamiltonian(i,i-1)    = NondiagConst;
    Hamiltonian(i,i)    = DiagConst + w(i);
    Hamiltonian(i,i+1)    = NondiagConst;
  }
  Hamiltonian(Dim-1,Dim-2) = NondiagConst;
  Hamiltonian(Dim-1,Dim-1) = DiagConst + w(Dim-1);
  // diagonalize and obtain eigenvalues
  vec Eigval(Dim);
  eig_sym(Eigval, Hamiltonian);
  output(RMin , RMax, Dim, Eigval);
  ifile.close();  // close input file
  ofile.close();  // close output file
  return 0;
}  //  end of main function


/*
  The function potential()
  calculates and return the value of the 
  potential for a given argument x.
  The potential here is for the hydrogen atom
*/        

double potential(double x)
{
  return x*x;

} // End: function potential()  


void output(double RMin , double RMax, int Dim, vec& d)
{
  int i;
  ofile << "RESULTS:" << endl;
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile <<"Rmin = " << setw(15) << setprecision(8) << RMin << endl;  
  ofile <<"Rmax = " << setw(15) << setprecision(8) << RMax << endl;  
  ofile <<"Number of steps = " << setw(15) << Dim << endl;  
  ofile << "Five lowest eigenvalues:" << endl;
  for(i = 0; i < 5; i++) {
    ofile << setw(15) << setprecision(8) << d[i] << endl;
  }
}  // end of function output         

