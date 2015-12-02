#include <fstream>
#include <cmath>

using namespace std;
void k(double* const y, double* k1, double* k2, double* k3, const double eta, const double dx);
void f(double* const y, const double eta, double*);

int main(){
 const double dx=0.001; const double psi0=1e-5; const double eta=1.0/2.0; const double x0=0; const double xf=100; const double n=(xf-x0)/dx;
 double y[2],k1[2],k2[2],k3[2]; y[0]=psi0; y[1]=sqrt(eta)*psi0; k1[0]=0;k1[1]=0;
 ofstream out("Data.txt");
 
 for(int i=0;i<n;i++){
  out << i*dx << " " << y[0] << " " << y[1] << " " << sqrt(2*eta)/cosh(sqrt(eta)*(i*dx-17.2)) << endl;
  k(y,k1,k2,k3,eta,dx);
  y[0]+= dx/6.0*(k1[0]+4*k2[0]+k3[0]);
  y[1]+= dx/6.0*(k1[1]+4*k2[1]+k3[1]);
 }
 
 out.close();
 return 0;
}

void k(double* const y, double* k1, double* k2, double* k3, const double eta, const double dx){
  f(y,eta,k1);
  double kk[2]; kk[0]=y[0]+dx*1/2.0*k1[0]; kk[1]=y[1]+dx*1/2.0*k1[1];
  f(kk,eta,k2);
  kk[0]=y[0]-dx*k1[0]+2.0*dx*k2[0]; kk[1]=y[1]-dx*k1[1]+2.0*dx*k2[1];
  f(kk,eta,k3);
}

void f(double* const y, const double eta, double* k){
  k[0]=y[1];
  k[1]=(eta-pow(y[0],2))*y[0];
}