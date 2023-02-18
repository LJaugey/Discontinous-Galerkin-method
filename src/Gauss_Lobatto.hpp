#ifndef GAUSS_LOBATTO_HPP
#define GAUSS_LOBATTO_HPP

#include <string>

#include "parameters.hpp"

grid_vec_qq Gauss_Lobatto_2D(Array<Nx+1> x, Array<Ny+1> y, Array<System_dim> (*f)(double,std::string, double,std::string));

Array<order> lgl_nodes();
Array<order> lgl_weights(Array<order> x);

#endif
