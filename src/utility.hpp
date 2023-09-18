#ifndef UTILITY_HPP
#define UTILITY_HPP

#include "parameters.hpp"

ND::Array<double, Ny,Nx, System_dim> quadrature_average_2D(grid_vec_qq const& u);

ND::Array<double, Ny,Nx> quadrature_average_2D(grid_point_qq const& nu);

#endif