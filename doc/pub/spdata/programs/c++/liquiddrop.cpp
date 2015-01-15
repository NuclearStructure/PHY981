// Liquid drop model code
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

// Begin of main program   

int main(int argc, char* argv[])
{
  char *outfilename, *infilename;
  int i, n, z, j, k, ndata, maxa, maxz, maxn;
  double bedata, error, a1, a2, a3, a4, a;

  a1 = 15.8; a2 = 18.3; a3 = 0.714; a4 = 23.2;
  ndata = 2375; maxz= 110; maxn = 159;  maxa = 269; 
  mat liquiddrop(maxz+1, maxn+1), be(maxz+1, maxn+1); 
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

  //  Read in data
  for( i = 1; i <= ndata; i++) {
    ifile >>  z;  ifile >> a;  ifile >> bedata;  ifile >> error; 
    ifile >>j; ifile >> n;
    be(z,n) = bedata/a;
  }

  // Compute separation energies and shell gaps 
  for( i = 1; i <=  maxz; i++){
    for( j = 1; j<=  maxn; j++){
      a = i+j;
      // Liquid drop model
      liquiddrop(i,j) = a1-a2*pow(a,0.66666)/a -a3*i*i*pow(a,-4.0/3.0)-a4*(j-i)*(j-i)*pow(a,-2.0);
    }
  }
  // Print the results
  for( j = 1; j <=  maxz; j++){
    for( k = 1; k <=  maxn; k++){
      if (be(j,k) != 0.0) {
	double error = abs(be(j,k)-liquiddrop(j,k));
	ofile << setiosflags(ios::showpoint | ios::uppercase);
	ofile << setw(15) << setprecision(8) << j;
	ofile << setw(15) << setprecision(8) << k;
	ofile << setw(15) << setprecision(8) << k+j;
	ofile << setw(15) << setprecision(8) << be(j,k);
	ofile << setw(15) << setprecision(8) << liquiddrop(j,k);
	ofile << setw(15) << setprecision(8) << error << endl;
      }
    }
  }
  ifile.close();  // close input file
  ofile.close();  // close output file
  return 0;
}  //  end of main function









