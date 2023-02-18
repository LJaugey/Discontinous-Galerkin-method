#ifndef ENTROPY_VARIABLES_HPP
#define ENTROPY_VARIABLES_HPP

#include "parameters.hpp"


vec_qq entropy_var_Custom(vec_qq const& u);
vec_qq entropy_var_Transport(vec_qq const& u);
vec_qq entropy_var_Burger(vec_qq const& u);
vec_qq entropy_var_Euler(vec_qq const& u);

#endif