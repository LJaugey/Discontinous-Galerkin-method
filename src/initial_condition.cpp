#include "initial_condition.hpp"

#include <string>
#include <iostream>

using namespace std;

ND::Array<System_dim> u0_custom(double x, string limit1, double y, string limit2)
{
    ND::Array<System_dim> result;

    return result;
}


ND::Array<System_dim> u0_sin(double x, string limit1, double y, string limit2)
{
    ND::Array<System_dim> result;

    result(0) = sin(2*M_PI*x);

    return result;
}


ND::Array<System_dim> u0_pos_sin(double x, string limit1, double y, string limit2)
{
    ND::Array<System_dim> result;

    result(0) = 1+sin(2*M_PI*x);

    return result;
}


ND::Array<System_dim> u0_Euler_smooth(double x, string limit1, double y, string limit2)
{
    ND::Array<System_dim> result;
    
    // X direction
    //result(0) = 1+0.2*sin(2*M_PI*x);
    //result(1) = 1+0.2*sin(2*M_PI*x);
    //result(2) = 0;
    //result(3) = (1+0.2*sin(2*M_PI*x))*(0.5 + 1.0/((gamma_-1)*(1+0.2*sin(2*M_PI*x))));
    
    // Y direction
    //result(0) = 1+0.2*sin(2*M_PI*y);
    //result(1) = 0;
    //result(2) = 1+0.2*sin(2*M_PI*y);
    //result(3) = (1+0.2*sin(2*M_PI*y))*(0.5 + 1.0/((gamma_-1)*(1+0.2*sin(2*M_PI*y))));

    // 2D test
    //double vx = 0.4;
    //double vy = 0.6;
    //result(0) = 1+0.2*sin(2*M_PI*(x+y));
    //result(1) = vx*(1+0.2*sin(2*M_PI*(x+y)));
    //result(2) = vy*(1+0.2*sin(2*M_PI*(x+y)));
    //result(3) = 0.5*result(0)*(vx*vx+vy*vy) + 1.0/(gamma_-1);

    double r2 = pow(x+5,2) + pow(y+5.0,2);
    double T = 1-25.0*(gamma_-1)*exp(1-r2)/(8.0*gamma_*M_PI*M_PI);

    double vx = 1-5.0*exp(0.5*(1-r2))*(y+5.0)/(2*M_PI);
    double vy = 1+5.0*exp(0.5*(1-r2))*(x+5.0)/(2*M_PI);

    double rho = pow(T,1.0/(gamma_-1.0));
    double p = pow(rho,gamma_);

    result(0) = rho;
    result(1) = rho*vx;
    result(2) = rho*vy;
    result(3) = 0.5*result(0)*(vx*vx+vy*vy) + p/(gamma_-1.0);

    return result;
}


ND::Array<System_dim> u0_Sod(double x, string limit1, double y, string limit2)
{
    ND::Array<System_dim> result;
    
    if(x<0.5)
    {
        result(0) = 1.0;
        result(1) = 0.0;
        result(2) = 0.0;
        result(3) = 1.0/(gamma_-1);
    }
    else if(x>0.5)
    {
        result(0) = 0.125;
        result(1) = 0.0;
        result(2) = 0.0;
        result(3) = 0.1/(gamma_-1);
    }
    else
    {
        if(limit1=="end")
        {
            result(0) = 1.0;
            result(1) = 0.0;
            result(2) = 0.0;
            result(3) = 1.0/(gamma_-1);
        }
        else if(limit1=="start")
        {
            result(0) = 0.125;
            result(1) = 0.0;
            result(2) = 0.0;
            result(3) = 0.1/(gamma_-1);
        }
        else    //use the average otherwise
        {
            result(0) = 0.5*(1.0 + 0.125);
            result(1) = 0.5*(0.0 + 0.0);
            result(2) = 0.5*(0.0 + 0.0);
            result(3) = 0.5*(1.0/(gamma_-1) + 0.1/(gamma_-1));
        }
    }
    
    return result;
}



ND::Array<System_dim> Riemann_2D_1(double x, string limit1, double y, string limit2)
{
    ND::Array<System_dim> result;

    double v1 = 0.0;
    double v2 = 0.0;
    double P = 0.0;
    
    if((x<0.5 || (x==0.5 && limit1=="end")) && (y<0.5 || (y==0.5 && limit2=="end")))
    {
        result(0) = 0.8;
        v1 = 0.0;
        v2 = 0.0;
        P = 1.0;
    }
    else if((x>0.5 || (x==0.5 && limit1=="start")) && (y<0.5 || (y==0.5 && limit2=="end")))
    {
        result(0) = 1.0;
        v1 = 0.0;
        v2 = 0.7276;
        P = 1.0;
    }
    else if((x<0.5 || (x==0.5 && limit1=="end")) && (y>0.5 || (y==0.5 && limit2=="start")))
    {
        result(0) = 1.0;
        v1 = 0.7276;
        v2 = 0.0;
        P = 1.0;
    }
    else if((x>0.5 || (x==0.5 && limit1=="start")) && (y>0.5 || (y==0.5 && limit2=="start")))
    {
        result(0) = 0.5313;
        v1 = 0.0;
        v2 = 0.0;
        P = 0.4;
    }
    else
    {
        cerr<<"Warning, some initial values are not initialized.\n Make sure to use even Nx and Ny."<<endl;
    }

    result(1) = result(0)*v1;
    result(2) = result(0)*v2;
    result(3) = 0.5*result(0)*(v1*v1+v2*v2) + P/(gamma_-1);

    return result;
}




ND::Array<System_dim> Riemann_2D_2(double x, string limit1, double y, string limit2)
{
    ND::Array<System_dim> result;

    //cout<<x<<","<<y<<endl;
    double v1 = 0.0;
    double v2 = 0.0;
    double P = 0.0;
    
    if((x<0.5 || (x==0.5 && limit1=="end")) && (y<0.5 || (y==0.5 && limit2=="end")))
    {
        result(0) = 1.0;
        v1 = -0.75;
        v2 = 0.5;
        P = 1.0;
    }
    else if((x>0.5 || (x==0.5 && limit1=="start")) && (y<0.5 || (y==0.5 && limit2=="end")))
    {
        result(0) = 3.0;
        v1 = -0.75;
        v2 = -0.5;
        P = 1.0;
    }
    else if((x<0.5 || (x==0.5 && limit1=="end")) && (y>0.5 || (y==0.5 && limit2=="start")))
    {
        result(0) = 2.0;
        v1 = 0.75;
        v2 = 0.5;
        P = 1.0;
    }
    else if((x>0.5 || (x==0.5 && limit1=="start")) && (y>0.5 || (y==0.5 && limit2=="start")))
    {
        result(0) = 1.0;
        v1 = 0.75;
        v2 = -0.5;
        P = 1.0;
    }
    else
    {
        cerr<<"Warning, some initial values are not initialized.\n Make sure to use even Nx and Ny."<<endl;
    }

    result(1) = result(0)*v1;
    result(2) = result(0)*v2;
    result(3) = 0.5*result(0)*(v1*v1+v2*v2) + P/(gamma_-1);

    return result;
}

ND::Array<System_dim> Riemann_2D_3(double x, string limit1, double y, string limit2)
{
    ND::Array<System_dim> result;

    double v1 = 0.0;
    double v2 = 0.0;
    double P = 0.0;
    
    if((x<0.5 || (x==0.5 && limit1=="end")) && (y<0.5 || (y==0.5 && limit2=="end")))
    {
        result(0) = 1.1;
        v1 = 0.8939;
        v2 = 0.8939;
        P = 1.1;
    }
    else if((x>0.5 || (x==0.5 && limit1=="start")) && (y<0.5 || (y==0.5 && limit2=="end")))
    {
        result(0) = 0.5065;
        v1 = 0.0;
        v2 = 0.8939;
        P = 0.35;
    }
    else if((x<0.5 || (x==0.5 && limit1=="end")) && (y>0.5 || (y==0.5 && limit2=="start")))
    {
        result(0) = 0.5065;
        v1 = 0.8939;
        v2 = 0.0;
        P = 0.35;
    }
    else if((x>0.5 || (x==0.5 && limit1=="start")) && (y>0.5 || (y==0.5 && limit2=="start")))
    {
        result(0) = 1.1;
        v1 = 0.0;
        v2 = 0.0;
        P = 1.1;
    }
    else
    {
        cerr<<"Warning, some initial values are not initialized.\n Make sure to use even Nx and Ny."<<endl;
    }

    result(1) = result(0)*v1;
    result(2) = result(0)*v2;
    result(3) = 0.5*result(0)*(v1*v1+v2*v2) + P/(gamma_-1);

    return result;
}

ND::Array<System_dim> Riemann_2D_4(double x, string limit1, double y, string limit2)
{
    ND::Array<System_dim> result;

    double v1 = 0.0;
    double v2 = 0.0;
    double P = 0.0;
    
    if((x<0.5 || (x==0.5 && limit1=="end")) && (y<0.5 || (y==0.5 && limit2=="end")))
    {
        result(0) = 0.1379;
        v1 = 1.206;
        v2 = 1.206;
        P = 0.029;
    }
    else if((x>0.5 || (x==0.5 && limit1=="start")) && (y<0.5 || (y==0.5 && limit2=="end")))
    {
        result(0) = 0.5323;
        v1 = 0.0;
        v2 = 1.206;
        P = 0.3;
    }
    else if((x<0.5 || (x==0.5 && limit1=="end")) && (y>0.5 || (y==0.5 && limit2=="start")))
    {
        result(0) = 0.5323;
        v1 = 1.206;
        v2 = 0.0;
        P = 0.3;
    }
    else if((x>0.5 || (x==0.5 && limit1=="start")) && (y>0.5 || (y==0.5 && limit2=="start")))
    {
        result(0) = 1.5;
        v1 = 0.0;
        v2 = 0.0;
        P = 1.5;
    }
    else
    {
        cerr<<"Warning, some initial values are not initialized.\n Make sure to use even Nx and Ny."<<endl;
    }

    result(1) = result(0)*v1;
    result(2) = result(0)*v2;
    result(3) = 0.5*result(0)*(v1*v1+v2*v2) + P/(gamma_-1);

    return result;
}



ND::Array<System_dim> Heat(double x, string limit1, double y, string limit2)
{
    ND::Array<System_dim> result;
    
    result(0) = cos(2*M_PI*x);
    
    return result;
}


ND::Array<System_dim> (*init_u0())(double, string, double, string)
{
    if(u0_Name == "Custom")
    {
        return u0_custom;
    }
    else if(u0_Name == "sin")
    {
        return u0_sin;
    }
    else if(u0_Name == "pos_sin")
    {
        return u0_pos_sin;
    }
    else if(u0_Name == "Euler_smooth")
    {
        return u0_Euler_smooth;
    }
    else if(u0_Name == "Sod")
    {
        return u0_Sod;
    }
    else if(u0_Name == "Riemann_2D_1")
    {
        return Riemann_2D_1;
    }
    else if(u0_Name == "Riemann_2D_2")
    {
        return Riemann_2D_2;
    }
    else if(u0_Name == "Riemann_2D_3")
    {
        return Riemann_2D_3;
    }
    else if(u0_Name == "Riemann_2D_4")
    {
        return Riemann_2D_4;
    }
    else if(u0_Name == "Heat")
    {
        return Heat;
    }
    else
    {
        cerr<<"Error, the specified initial condition is not known"<<endl;
    }
    return u0_custom;
}