import matplotlib.pyplot as plt
import numpy as np
import sys
import os


# matplotlib latex interpreter
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

if(len(sys.argv)>1):
    Equation_Name = sys.argv[1]

if(Equation_Name=="Euler"):
    System_dim = 4
    gamma_ = 7/5
elif(Equation_Name=="Transport" or Equation_Name=="Burger"):
    System_dim = 1



ax = 0
bx = 1

ay = 0
by = 1

t_fin = 0.2

L1_err = []
L2_err = []
Linf_err = []

dx_list = []

Nx_min = 10
Nx_max = 40
delta_Nx = 10

ana_int = lambda x,t: x - 0.2*np.cos(2*np.pi*(x-t))/(2*np.pi)

#os.system("make clean")
#os.system("make")

for i, sim in enumerate(range(Nx_min,Nx_max+1,delta_Nx)):

    print(i)
    print(sim)

    #os.system("./run Nx "+str(sim)+" "+str(i))
    #os.system("srun -n 1 -c 28 ./run Nx "+str(sim)+" "+str(i))

    Nx = sim
    Ny = 1

    dx = (bx-ax)/Nx
    dy = (by-ay)/Ny

    dx_list.append(dx)



    ### u

    data_u = open("output/u"+str(i)+".out", 'r').read()

    Nt = len(data_u.split('\n\n\n'))-1

    u_flatten = data_u.split()

    u = np.empty((Nt,Ny,Nx,System_dim))


    for it in range(Nt):
        for i in range(Ny):
            for j in range(Nx):
                for k in range(System_dim):
                
                    u[it,i,j,k] = float(u_flatten[it*(Nx*Ny)*(System_dim) + i*Nx*System_dim + j*System_dim + k])


    u = np.where(np.isnan(u), np.nanmin(u), u)

    x_ = np.linspace(ax+0.5*dx, bx-0.5*dx, Nx)
    y_ = np.linspace(ay+0.5*dy, by-0.5*dy, Ny)

    
    ana = (ana_int(x_+0.5*dx, t_fin) - ana_int(x_-0.5*dx, t_fin))/dx

    diff = u[0,0,:,0]-ana

    
    L1_err.append(np.sum(np.abs(diff))*dx)
    L2_err.append(np.sqrt(np.sum(diff**2))*dx)
    Linf_err.append(np.max(np.abs(diff)))

    #if(sim == Nx_max):
    #
    #    plt.figure()
    #    plt.plot(x_,u[0,:,:,0]-ana)


dx_list = np.array(dx_list)

print(dx_list)
print()

L1_err = np.array(L1_err)
L2_err = np.array(L2_err)
Linf_err = np.array(Linf_err)


print(L1_err)
print(L2_err)
print(Linf_err)
print()


print(np.log(L1_err[1:]/L1_err[0:-1])/(np.log(dx_list[1:]/dx_list[0:-1])))
print(np.log(L2_err[1:]/L2_err[0:-1])/(np.log(dx_list[1:]/dx_list[0:-1])))
print(np.log(Linf_err[1:]/Linf_err[0:-1])/(np.log(dx_list[1:]/dx_list[0:-1])))
print()

plt.figure()
plt.loglog(dx_list, L1_err)
plt.grid(which='both')


plt.figure()
plt.loglog(dx_list, L2_err)
plt.grid(which='both')


plt.figure()
plt.loglog(dx_list, Linf_err)
plt.grid(which='both')


plt.show()
