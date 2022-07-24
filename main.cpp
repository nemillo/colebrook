/*
Colebrook-White Equation Solver
Newton-Cotes Method

Rafael OLONA POBLET N° Matricula:21299  Colegiado N° 25156
*/ 
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <cstdlib>

int main(int argc, char *argv[])
{
  
  double Q=1.07; //Flow in m3/s
  double D=1.0857; //Internal Diameter in m
  double epsilon = 0.002; // Absolute Roughness in m
  double nu = 0.0000013; // Kinematic Viscosity in m2/s
     
  double S = M_PI*D*D/4; //Internal Surface
  double vel = Q/S; //Flow Velocity
  double Re = vel*D/nu; //Reynolds Number
     
  std::cout << "Reynolds Number = "<< Re << "\n";
  std::cout << "Relative Roughness = " << epsilon/D << "\n" ;
  std::cout << "Flow Velocity = " << vel << " m/s" << "\n";
  std::cout << "Flow Surface = " << S << " m^2" << "\n" << "\n" ;
  
      if (Re < 2300.0) {
          std::cout << "f Darcy Factor (Laminar Regime)" << 64/Re << "\n";
      }
      else {
          double er = 1;
          double der = 1;          
          double C1 = 2.51/Re;
          double C2 = epsilon/(D*3.715);
          double x = -2*log10(C2); // x = 1/sqrt(f) ; x0 = -2*log10(C2) 
          int n = 0;
          while(abs(er) > .000001) {
              n = n +1;
              er = x + 2*log10(C1*x+C2); // Error -> x - (-2*log10(C1*x+C2))
              der = 1 + 2*C1/(C1*x+C2); // Error derivative ->  1 - (-2*C2/(C1*x+C2))
              x = x -er/der; //xn+1 = xn - error/error derivative (Newton-Cotes method)
              std::cout << "Iteration nr " << n << "  Error: " << er << "  Error Derivative: " << der << "  f Darcy Factor: " << 1/(x*x) << "\n" << "\n";
          }
          std::cout << "Total Iterations = " << n << "\n";
          std::cout << "f Darcy Factor (Turbulent Regime) = " << 1/(x*x) << "\n";
      }        
                          
    
  system("PAUSE");	
  return 0;
}
