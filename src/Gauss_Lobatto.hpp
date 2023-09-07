#ifndef GAUSS_LOBATTO_HPP
#define GAUSS_LOBATTO_HPP

#include <string>

#include "parameters.hpp"

grid_vec_qq Gauss_Lobatto_2D(ND::Array<Nx+1> x, ND::Array<Ny+1> y, ND::Array<System_dim> (*f)(double,std::string, double,std::string));

ND::Array<order> lgl_nodes();
ND::Array<order> lgl_weights(ND::Array<order> x);

#endif
