#include "parameters.hpp"
#include "init.hpp"
#include "Gauss_Lobatto.hpp"
#include "Equation.hpp"
#include "entropy_variables.hpp"
#include "Artificial_viscosity.hpp"

using namespace std;

Array<order> x_ref = lgl_nodes();
Array<order> w_ref = lgl_weights(x_ref);
Array<order,order> D = init_D(x_ref);

double t;
double dt_old = t_fin;


vec_ (*f1)(vec_ const&) = f1_Euler;
vec_ (*f2)(vec_ const&) = f2_Euler;
point_Jac (*df1)(vec_ const&) = df1_Euler;
point_Jac (*df2)(vec_ const&) = df2_Euler;
double (*E)(grid_vec_qq const&) = entropy_Euler;
void (*EC_flux1)(vec_ const&, vec_ const&, vec_ &) = EC_flux_Euler1;
void (*EC_flux2)(vec_ const&, vec_ const&, vec_ &) = EC_flux_Euler2;
vec_ (*int_flux1)(vec_ const& ,vec_ const&) = LLF1_Euler;
vec_ (*int_flux2)(vec_ const& ,vec_ const&) = LLF2_Euler;
vec_qq (*entropy_var)(vec_qq const&) = entropy_var_Euler;


void (*artificial_viscosity)(grid_vec_qq const& u, grid_point_qq & nu) = dilation_based_viscosity_Euler;
//void (*artificial_viscosity)(grid_vec_qq const& u, grid_point_qq & nu) = entropy_based_viscosity_Euler;


// Entropy based viscosity
grid_vec_qq u_old;
