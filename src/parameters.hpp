#define _USE_MATH_DEFINES
#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP


#include <cmath>
#include <string>
#include "ND_Array/ND_Array.hpp"

using namespace std;







// Spatial domain
constexpr double ax = 0.0;
constexpr double bx = 1.0;

constexpr double ay = 0.0;
constexpr double by = 1.0;



// Time domain
constexpr double t0 = 0.0;


// Equation
const string BC = "outflow";
const string Equation_Name = "Euler";
const string int_flux_Name = "LLF";
//const string u0_Name = "Riemann_2D_1";  constexpr double t_fin = 0.25;
const string u0_Name = "Riemann_2D_2";  constexpr double t_fin = 0.3;
//const string u0_Name = "Riemann_2D_3";  constexpr double t_fin = 0.25;
//const string u0_Name = "Riemann_2D_4";    constexpr double t_fin = 1.1;



constexpr int System_dim = 4;


// Only for Euler equation
constexpr double gamma_ = 7.0/5.0;





// Scheme parameters
constexpr int order = 3;
constexpr int Ny = 75;
constexpr int Nx = 75;
constexpr double dx = (bx-ax)/Nx;
constexpr double dy = (by-ay)/Ny;
extern double dt_old;
//double CFL = 2.0/(2*order-1);
constexpr double CFL = 0.5;
constexpr bool viscosity_bool = false;
constexpr std::string Artificial_Viscosity_Name = "Dialation";    constexpr double c_vis = 1.0; constexpr double c_max = 0.5;
//std::string Artificial_Viscosity_Name = "Entropy";    double c_vis = 1.0; double c_max = 0.5;





typedef ND::Array<double, Ny+2,Nx+2, order,order, System_dim> grid_vec_qq;
typedef ND::Array<double, Ny+2,Nx+2, order, System_dim> grid_vec_q;
typedef ND::Array<double, Ny+2,Nx+2, System_dim> grid_vec_;

typedef ND::Array<double, order,order, System_dim> vec_qq;
typedef ND::Array<double, order, System_dim> vec_q;
typedef ND::Array<double, System_dim> vec_;

typedef ND::Array<double, Ny+2,Nx+2, order,order> grid_point_qq;
typedef ND::Array<double, Ny+2,Nx+2, order> grid_point_q;
typedef ND::Array<double, Ny+2,Nx+2> grid_point_;

typedef ND::Array<double, order,order> point_qq;
typedef ND::Array<double, order> point_q;
//typedef ND::Array<double, > point_;

typedef ND::Array<double, System_dim,System_dim> point_Jac;






// Saving data
constexpr bool Save_data = true;
constexpr int sampling = 20;         // ignored if Save_data==false
constexpr bool Save_data_first = true;
constexpr bool Save_data_last = true;



constexpr std::string exp_name = "";
extern double t;




extern ND::Array<double, order> x_ref;
extern ND::Array<double, order> w_ref;
extern ND::Array<double, order,order> D;


extern vec_ (*f1)(vec_ const&);
extern vec_ (*f2)(vec_ const&);
extern point_Jac (*df1)(vec_ const&);
extern point_Jac (*df2)(vec_ const&);






extern double (*E)(grid_vec_qq const&);

// fluxes
extern void (*EC_flux1)(vec_ const&, vec_ const&, vec_ &);
extern void (*EC_flux2)(vec_ const&, vec_ const&, vec_ &);
extern vec_ (*int_flux1)(vec_ const& ,vec_ const&);
extern vec_ (*int_flux2)(vec_ const& ,vec_ const&);
extern vec_qq (*entropy_var)(vec_qq const&);

extern void (*artificial_viscosity)(grid_vec_qq const& u, grid_point_qq & nu);
//extern void (*artificial_viscosity)(grid_vec_qq const& u, grid_point_qq & nu);


// Entropy based viscosity
extern grid_vec_qq u_old;



constexpr double power_ln(double aL,double aR)
{
    //res = (aR-aL)/(log(aR)-log(aL));

    double t = aL/aR;
    double f = (t-1)/(t+1);
    double u = f*f;

    double F;

    if(u < 1e-2) F = 1.+u/3.+u*u/5.+u*u*u/7.;
    else F = log(t)/2/f;

    return (aL+aR)/2/F;
}


#endif