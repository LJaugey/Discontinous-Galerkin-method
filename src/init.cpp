#include "init.hpp"


using namespace std;

Array<order,order> init_D(Array<order> x_ref)
{
    Array<order,order> D_;

    for(int i=0; i<order; i++)
    {
        for(int j=0; j<order; j++)
        {
            if(i==j)
            {
                for(int k=0; k<order; k++)
                {
                    if(k!=j)    D_(i,j) += 1.0/(x_ref(i)-x_ref(k));
                }
            }
            else
            {
                D_(i,j) = 1.0/(x_ref(j)-x_ref(i));

                for(int k=0; k<order; k++)
                {
                    if((k != i) && (k != j))    D_(i,j) *= (x_ref(i)-x_ref(k))/(x_ref(j)-x_ref(k));
                }
            }
        }
    }
    return D_;
}


/*
void init_Equation()
{
    //else if(Artificial_Viscosity_Name=="Entropy")
    //{
    //    artificial_viscosity = entropy_based_viscosity;
    //}
    //if(Equation_Name=="Transport")
    //{
    //    System_dim = 1;
    //    E = entropy_Transport;
    //    f1 = f1_Transport;
    //    f2 = f2_Transport;
    //    df1 = df1_Transport;
    //    df2 = df2_Transport;
    //    EC_flux1 = EC_flux_Transport1;
    //    EC_flux2 = EC_flux_Transport2;
    //    if(int_flux_Name == "LLF")
    //    {
    //        int_flux1 = LLF1_scal;
    //        int_flux2 = LLF2_scal;
    //    }
    //    else if(int_flux_Name == "EC_Diffusion")
    //    {
    //        int_flux1 = EC_Diffusion_scal;
    //        int_flux2 = EC_Diffusion_scal;
    //    }
    //    entropy_var = entropy_var_Transport;
    //}
    //else if(Equation_Name=="Burger")
    //{
    //    System_dim = 1;
    //    E = entropy_Burger;
    //    f1 = f1_Burger;
    //    f2 = f2_Burger;
    //    df1 = df1_Burger;
    //    df2 = df2_Burger;
    //    EC_flux1 = EC_flux_Burger1;
    //    EC_flux2 = EC_flux_Burger2;
    //    if(int_flux_Name == "LLF")
    //    {
    //        int_flux1 = LLF1_scal;
    //        int_flux2 = LLF2_scal;
    //    }
    //    else if(int_flux_Name == "EC_Diffusion")
    //    {
    //        int_flux1 = EC_Diffusion_scal;
    //        int_flux2 = EC_Diffusion_scal;
    //    }
    //    entropy_var = entropy_var_Burger;
    //}
    //else if(Equation_Name=="Euler")
    if(Equation_Name=="Euler")
    {
        System_dim = 4;
        E = entropy_Euler;
        f1 = f1_Euler;
        f2 = f2_Euler;
        df1 = df1_Euler;
        df2 = df2_Euler;
        EC_flux1 = EC_flux_Euler1;
        EC_flux2 = EC_flux_Euler2;
        if(int_flux_Name == "LLF")
        {
            int_flux1 = LLF1_Euler;
            int_flux2 = LLF2_Euler;
        }
        else if(int_flux_Name == "EC_Diffusion")
        {
            //int_flux1 = EC_Diffusion_Euler;
            //int_flux2 = EC_Diffusion_Euler;
        }
        entropy_var = entropy_var_Euler;
        if(Artificial_Viscosity_Name=="Dialation")  artificial_viscosity = dilation_based_viscosity_Euler;
        else if(Artificial_Viscosity_Name=="Entropy")   artificial_viscosity = entropy_based_viscosity_Euler;
    }
    //else if(Equation_Name=="Custom")
    //{
    //    System_dim = 1;
    //    E = entropy_Custom;
    //    f1 = f1_Custom;
    //    f2 = f2_Custom;
    //    df1 = df1_Custom;
    //    df2 = df2_Custom;
    //    EC_flux1 = EC_flux_Custom1;
    //    EC_flux2 = EC_flux_Custom2;
    //    if(int_flux_Name == "LLF")
    //    {
    //        int_flux1 = LLF1_Custom;
    //        int_flux2 = LLF2_Custom;
    //    }
    //    else if(int_flux_Name == "EC_Diffusion")
    //    {
    //        //int_flux1 = EC_Diffusion_Custom;
    //        //int_flux2 = EC_Diffusion_Custom;
    //    }
    //    entropy_var = entropy_var_Custom;
    //}
}
*/