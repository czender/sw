// Program to demonstrate g++ complex division 
// g++ complex.cc -o complex
// xlC complex.cc -o complex
#include <iostream>
#include <complex>
int main(const int argc,char **argv){
  std::complex<int> icx(1,1);
  std::complex<float> fcx(1,1);
  std::complex<double> dcx(1,1);
  std::cout << "Hello, World" << std::endl;
  std::complex<float>dof_flt=static_cast<std::complex<float> >(dcx)/fcx; // Not supported by egcs
  std::complex<double>dof_dbl=dcx/static_cast<std::complex<double> >(fcx); // Not supported by egcs
  //  std::cout << "icx = " << icx << std::endl;
  std::cout << "std::sqrt(fcx) = " << std::sqrt(fcx) << std::endl;
  std::cout << "dof_flt = " << dof_flt << std::endl;
  std::cout << "dof_dbl = " << dof_dbl << std::endl;
}








