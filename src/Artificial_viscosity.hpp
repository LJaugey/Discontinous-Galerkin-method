#ifndef ARTIFICIAL_VISCOSITY_HPP
#define ARTIFICIAL_VISCOSITY_HPP

#include "parameters.hpp"

void dilation_based_viscosity_Euler(grid_vec_qq const& u, Array<Ny+2,Nx+2, order,order> & nu);
void entropy_based_viscosity_Euler(grid_vec_qq const& u, Array<Ny+2,Nx+2, order,order> & nu);

#endif