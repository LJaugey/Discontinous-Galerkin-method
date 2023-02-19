#ifndef ARTIFICIAL_VISCOSITY_HPP
#define ARTIFICIAL_VISCOSITY_HPP

#include "parameters.hpp"

void dilation_based_viscosity_Euler(grid_vec_qq const& u, grid_point_qq & nu);
void entropy_based_viscosity_Euler(grid_vec_qq const& u, grid_point_qq & nu);

#endif