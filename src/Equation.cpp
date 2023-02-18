#include "Equation.hpp"
#include "utility.hpp"

#include <iostream>

using namespace std;


// ======================== Custom equations ======================== //
double entropy_Custom(grid_vec_qq const& u)
{
    double res = 0.0;

    for(size_t i = 1; i<u.size()-1; i++)
    {
        for(size_t j = 1; j<u.size(1)-1; j++)
        {
            for(size_t l = 0; l<u.size(2); l++)
            {
                for(size_t m = 0; m<u.size(3); m++)
                {
                    res += 0.0;
                }
            }
        }
    }

    res *= 0.25*dy*dx/(Nx*Ny);

    return res;
}

// f
vec_ f1_Custom(vec_ const& u)
{
    return vec_();
}
vec_ f2_Custom(vec_ const& u)
{
    return vec_();
}


// df
point_Jac df1_Custom(vec_ const& u)
{
    return point_Jac();
}
point_Jac df2_Custom(vec_ const& u)
{
    return point_Jac();
}


// EC flux
vec_ EC_flux_Custom1(vec_ const& uL,vec_ const& uR)
{
    return vec_();
}
vec_ EC_flux_Custom2(vec_ const& uL,vec_ const& uR)
{
    return vec_();
}











// ======================== Transport equations ======================== //
double entropy_Transport(grid_vec_qq const& u)
{
    double res = 0.0;

    for(size_t i = 1; i<u.size()-1; i++)
    {
        for(size_t j = 1; j<u.size(1)-1; j++)
        {
            for(size_t l = 0; l<u.size(2); l++)
            {
                for(size_t m = 0; m<u.size(3); m++)
                {
                    res += w_ref(l)*w_ref(m)*u(i,j,l,m,0)*u(i,j,l,m,0);
                }
            }
        }
    }

    res *= 0.5*0.25*dy*dx/(Nx*Ny);

    return res;
}

// f
vec_ f1_Transport(vec_ const& u)
{
    return u;
}
vec_ f2_Transport(vec_ const& u)
{
    return u;
}


// df
point_Jac df1_Transport(vec_ const& u)
{
    point_Jac result;
    result.fill(1.0);

    return result;
}
point_Jac df2_Transport(vec_ const& u)
{
    point_Jac result;
    result.fill(1.0);

    return result;
}


// EC flux
vec_ EC_flux_Transport(vec_ const& uL,vec_ const& uR)
{
    return 0.5*(uL + uR);
}






// ======================== Burger's equations ======================== //
double entropy_Burger(grid_vec_qq const& u)
{
    double res = 0.0;

    for(size_t i = 1; i<u.size()-1; i++)
    {
        for(size_t j = 1; j<u.size(1)-1; j++)
        {
            for(size_t l = 0; l<u.size(2); l++)
            {
                for(size_t m = 0; m<u.size(3); m++)
                {
                    res += w_ref(l)*w_ref(m)*u(i,j,l,m,0)*u(i,j,l,m,0);
                }
            }
        }
    }

    res *= 0.5*0.25*dy*dx/(Nx*Ny);

    return res;
}

// f
vec_ f1_Burger(vec_ const& u)
{
    return 0.5*u*u;
}
vec_ f2_Burger(vec_ const& u)
{
    return 0.5*u*u;
}


// df
point_Jac df1_Burger(vec_ const& u)
{
    point_Jac result(u);

    return result;
}
point_Jac df2_Burger(vec_ const& u)
{
    point_Jac result(u);

    return result;
}


// EC flux
vec_ EC_flux_Burger(vec_ const& uL,vec_ const& uR)
{
    return (uL*uL + uL*uR + uR*uR)/6;
}





// ======================== Euler equations ======================== //
double entropy_Euler(grid_vec_qq const& u)
{
    double res = 0.0;

    double P = 0.0;
    double s = 0.0;

    for(size_t i = 1; i<u.size()-1; i++)
    {
        for(size_t j = 1; j<u.size(1)-1; j++)
        {
            for(size_t l = 0; l<u.size(2); l++)
            {
                for(size_t m = 0; m<u.size(3); m++)
                {
                    P = (gamma_-1)*(u(i,j,l,m,2)-0.5*(u(i,j,l,m,1)*u(i,j,l,m,1)/u(i,j,l,m,0)));
                    s = log(P)-gamma_*log(u(i,j,l,m,0));

                    res += w_ref(l)*w_ref(m)*(u(i,j,l,m,0)*s);
                }
            }
        }
    }

    res *= -0.25*dy*dx/((gamma_-1)*Nx*Ny);

    return res;
}

// f
vec_ f1_Euler(vec_ const& u)
{
    vec_ result;

    double p = (gamma_-1)*(u(3)-0.5*(u(1)*u(1) + u(2)*u(2))/u(0));

    result(0) = u(1);
    result(1) = u(1)*u(1)/u(0) + p;
    result(2) = u(1)*u(2)/u(0);
    result(3) = u(1)*(u(3)+p)/u(0);

    return result;
}

vec_ f2_Euler(vec_ const& u)
{
    vec_ result;

    double p = (gamma_-1)*(u(3)-0.5*(u(1)*u(1) + u(2)*u(2))/u(0));

    result(0) = u(2);
    result(1) = u(2)*u(1)/u(0);
    result(2) = u(2)*u(2)/u(0) + p;
    result(3) = u(2)*(u(3)+p)/u(0);

    return result;
}


// df
point_Jac df1_Euler(vec_ const& u)
{
    point_Jac result;

    double v1 = u(1)/u(0);
    double v2 = u(2)/u(0);

    double V2 = (v1*v1+v2*v2);
    double H = gamma_*u(3)/u(0) - 0.5*(gamma_-1)*V2;

    result(0,0) = 0; 
    result(0,1) = 1;
    result(0,2) = 0;
    result(0,3) = 0;

    result(1,0) = 0.5*(gamma_-1)*V2-v1*v1; 
    result(1,1) = (3-gamma_)*v1;
    result(1,2) = (1-gamma_)*v2;
    result(1,3) = (gamma_-1);

    result(2,0) = -v1*v2;
    result(2,1) = v2;
    result(2,2) = v1;
    result(2,3) = 0;

    result(3,0) = v1*(0.5*(gamma_-1)*V2-H);
    result(3,1) = H-(gamma_-1)*v1*v1;
    result(3,2) = (1-gamma_)*v1*v2;
    result(3,3) = gamma_*v1;

    return result;
}
point_Jac df2_Euler(vec_ const& u)
{
    point_Jac result;
    
    double v1 = u(1)/u(0);
    double v2 = u(2)/u(0);

    double V2 = (v1*v1+v2*v2);
    double H = gamma_*u(3)/u(0) - 0.5*(gamma_-1)*V2;

    result(0,0) = 0; 
    result(0,1) = 0;
    result(0,2) = 1;
    result(0,3) = 0;

    result(1,0) = -v1*v2;
    result(1,1) = v2;
    result(1,2) = v1;
    result(1,3) = 0;

    result(2,0) = 0.5*(gamma_-1)*V2-v2*v2; 
    result(2,1) = (1-gamma_)*v1;
    result(2,2) = (3-gamma_)*v2;
    result(2,3) = (gamma_-1);

    result(3,0) = v2*(0.5*(gamma_-1)*V2-H);
    result(3,1) = (1-gamma_)*v1*v2;
    result(3,2) = H-(gamma_-1)*v2*v2;
    result(3,3) = gamma_*v2;

    return result;
}


// EC flux
void EC_flux_Euler1(vec_ const& uL,vec_ const& uR, vec_ & result)
{
    // left
    //double pL = (gamma_-1)*(uL(3)-0.5*((uL(1)*uL(1)+uL(2)*uL(2))/uL(0)));

    //double zL1 = uL(0);
    double zL2 = uL(1)/uL(0);
    double zL3 = uL(2)/uL(0);
    //double zL4 = uL(0)/pL;
    //double zL4 = uL(0)/((gamma_-1)*(uL(3)-0.5*((uL(1)*uL(1)+uL(2)*uL(2))/uL(0))));
    double zL4 = uL(0)*uL(0)/((gamma_-1)*(uL(0)*uL(3)-0.5*(uL(1)*uL(1)+uL(2)*uL(2))));


    // right
    //double pR = (gamma_-1)*(uR(3)-0.5*((uR(1)*uR(1)+uR(2)*uR(2))/uR(0)));

    //double zR1 = uR(0);
    double zR2 = uR(1)/uR(0);
    double zR3 = uR(2)/uR(0);
    //double zR4 = uR(0)/pR;
    //double zR4 = uR(0)/((gamma_-1)*(uR(3)-0.5*((uR(1)*uR(1)+uR(2)*uR(2))/uR(0))));
    double zR4 = uR(0)*uR(0)/((gamma_-1)*(uR(0)*uR(3)-0.5*(uR(1)*uR(1)+uR(2)*uR(2))));


    //double z1av = 0.5*(zL1 + zR1);
    double z2av = 0.5*(zL2 + zR2);
    double z3av = 0.5*(zL3 + zR3);
    //double z4av = 0.5*(zL4 + zR4);

    //result(0) = z2av * power_ln(zL1,zR1);
    result(0) = z2av * power_ln(uL(0),uR(0));
    //result(1) = z1av/z4av + z2av*result(0);
    //result(1) = 0.5*(zL1 + zR1)/z4av + z2av*result(0);
    result(1) = (uL(0) + uR(0))/(zL4 + zR4) + z2av*result(0);
    result(2) = z3av*result(0);
    result(3) = (1.0/((gamma_-1)*power_ln(zL4,zR4)) - 0.25*(zL2*zL2+zR2*zR2 + zL3*zL3+zR3*zR3))*result(0) + z2av*result(1) + z3av*result(2);
}


void EC_flux_Euler2(vec_ const& uL,vec_ const& uR, vec_ & result)
{
    // left
    double pL = (gamma_-1)*(uL(3)-0.5*((uL(1)*uL(1)+uL(2)*uL(2))/uL(0)));

    double zL1 = uL(0);
    double zL2 = uL(1)/uL(0);
    double zL3 = uL(2)/uL(0);
    double zL4 = uL(0)/pL;

    // right
    double pR = (gamma_-1)*(uR(3)-0.5*((uR(1)*uR(1)+uR(2)*uR(2))/uR(0)));

    double zR1 = uR(0);
    double zR2 = uR(1)/uR(0);
    double zR3 = uR(2)/uR(0);
    double zR4 = uR(0)/pR;


    //double z1av = 0.5*(zL1 + zR1);
    double z2av = 0.5*(zL2 + zR2);
    double z3av = 0.5*(zL3 + zR3);
    double z4av = 0.5*(zL4 + zR4);

    result(0) = z3av * power_ln(zL1,zR1);
    result(1) = z2av*result(0);
    //double F3 = z1av/z4av + z3av*result(0);
    result(2) = 0.5*(zL1 + zR1)/z4av + z3av*result(0);
    result(3) = (1.0/((gamma_ -1)*power_ln(zL4,zR4)) - 0.25*(zL2*zL2+zR2*zR2 + zL3*zL3+zR3*zR3))*result(0) + z2av*result(1) + z3av*result(2);
}


vec_ LLF1_Euler(vec_ const& uL,vec_ const& uR)
{
    double v_L = uL(1)/uL(0);
    double p_L = (gamma_-1)*(uL(3)-0.5*(uL(1)*uL(1) + uL(2)*uL(2))/uL(0));
    double a_L = sqrt(gamma_*p_L/uL(0));

    double v_R = uR(1)/uR(0);
    double p_R = (gamma_-1)*(uR(3)-0.5*(uR(1)*uR(1) + uR(2)*uR(2))/uR(0));;
    double a_R = sqrt(gamma_*p_R/uR(0));

    double coeff = max(abs(v_L)+a_L,abs(v_R)+a_R);

    return 0.5*(f1(uL) + f1(uR) - coeff*(uR-uL));
}


vec_ LLF2_Euler(vec_ const& uL,vec_ const& uR)
{
    double v_L = uL(2)/uL(0);
    double p_L = (gamma_-1)*(uL(3)-0.5*(uL(1)*uL(1) + uL(2)*uL(2))/uL(0));
    double a_L = sqrt(gamma_*p_L/uL(0));

    double v_R = uR(2)/uR(0);
    double p_R = (gamma_-1)*(uR(3)-0.5*(uR(1)*uR(1) + uR(2)*uR(2))/uR(0));;
    double a_R = sqrt(gamma_*p_R/uR(0));

    double coeff = max(abs(v_L)+a_L,abs(v_R)+a_R);

    return 0.5*(f2(uL) + f2(uR) - coeff*(uR-uL));
}
