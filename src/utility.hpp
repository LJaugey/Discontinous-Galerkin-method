#ifndef UTILITY_HPP
#define UTILITY_HPP

#include "parameters.hpp"

Array<Ny,Nx, System_dim> quadrature_average_2D(grid_vec_qq const& u);

Array<Ny,Nx> quadrature_average_2D(grid_point_qq const& nu);

#endif