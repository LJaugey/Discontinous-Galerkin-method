#ifndef GAUSS_LOBATTO_HPP
#define GAUSS_LOBATTO_HPP

#include <string>

#include "parameters.hpp"

grid_vec_qq Gauss_Lobatto_2D(ND::Array<double, Nx+1> x, ND::Array<double, Ny+1> y, ND::Array<double, System_dim> (*f)(double,std::string, double,std::string));

ND::Array<double, order> lgl_nodes();
ND::Array<double, order> lgl_weights(ND::Array<double, order> x);

#endif
