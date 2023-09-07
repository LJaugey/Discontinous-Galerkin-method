#ifndef INITIAL_CONDITION_HPP
#define INITIAL_CONDITION_HPP

#include "parameters.hpp"
#include <string>

ND::Array<System_dim> u0_custom(double x, std::string limit);
ND::Array<System_dim> u0_sin(double x, std::string limit);
ND::Array<System_dim> u0_pos_sin(double x, std::string limit);
ND::Array<System_dim> u0_Euler_smooth(double x, std::string limit1, double y, std::string limit2);
ND::Array<System_dim> u0_Sod(double x, std::string limit1, double y, std::string limit2);
ND::Array<System_dim> Riemann_2D_1(double x, std::string limit1, double y, std::string limit2);
ND::Array<System_dim> Riemann_2D_2(double x, std::string limit1, double y, std::string limit2);
ND::Array<System_dim> Heat(double x, std::string limit1, double y, std::string limit2);

ND::Array<System_dim> (*init_u0())(double, std::string, double, std::string);

#endif
