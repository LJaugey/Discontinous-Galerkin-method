
#ifndef EQUATION_HPP
#define EQUATION_HPP

#include "parameters.hpp"


// ======================== Custom equations ======================== //
double entropy_Custom(grid_vec_qq const& u);
vec_ f1_Custom(vec_ const& u);
vec_ f2_Custom(vec_ const& u);
point_Jac df1_Custom(vec_ const& u);
point_Jac df2_Custom(vec_ const& u);
vec_ EC_flux_Custom1(vec_ const& uL,vec_ const& uR);
vec_ EC_flux_Custom2(vec_ const& uL,vec_ const& uR);

vec_ LLF1_Custom(vec_ const& uL,vec_ const& uR);
vec_ LLF2_Custom(vec_ const& uL,vec_ const& uR);
vec_ EC_Diffusion_Custom(vec_ const& uL,vec_ const& uR);

// ======================== Transport equations ======================== //
double entropy_Transport(grid_vec_qq const& u);
vec_ f1_Transport(vec_ const& u);
vec_ f2_Transport(vec_ const& u);
point_Jac df1_Transport(vec_ const& u);
point_Jac df2_Transport(vec_ const& u);
vec_ EC_flux_Transport1(vec_ const& uL,vec_ const& uR);
vec_ EC_flux_Transport2(vec_ const& uL,vec_ const& uR);

// ======================== Burger's equations ======================== //
double entropy_Burger(grid_vec_qq const& u);
vec_ f1_Burger(vec_ const& u);
vec_ f2_Burger(vec_ const& u);
point_Jac df1_Burger(vec_ const& u);
point_Jac df2_Burger(vec_ const& u);
vec_ EC_flux_Burger1(vec_ const& uL,vec_ const& uR);
vec_ EC_flux_Burger2(vec_ const& uL,vec_ const& uR);


// ======================== Euler equations ======================== //
double entropy_Euler(grid_vec_qq const& u);
vec_ f1_Euler(vec_ const& u);
vec_ f2_Euler(vec_ const& u);
point_Jac df1_Euler(vec_ const& u);
point_Jac df2_Euler(vec_ const& u);
void EC_flux_Euler1(vec_ const& uL,vec_ const& uR, vec_ & res);
void EC_flux_Euler2(vec_ const& uL,vec_ const& uR, vec_ & res);

vec_ LLF1_Euler(vec_ const& uL,vec_ const& uR);
vec_ LLF2_Euler(vec_ const& uL,vec_ const& uR);
vec_ EC_Diffusion_Euler(vec_ const& uL,vec_ const& uR);



#endif