import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os


# matplotlib latex interpreter
Latex_interpreter = False

# Run c++ code
Compute_bool = True

# Number of threads (for the c++ code)
thread_num = 16



Nx = 50
Ny = 50


#Prob_name = "RP_1"
Prob_name = "RP_2"
#Prob_name = "RP_3"

#P_order = "P1"
P_order = "P2"
#P_order = "P4"

vis_name = "no_vis"
#vis_name = "dilation"
#vis_name = "entropy"


ax = 0
bx = 1

ay = 0
by = 1


System_dim = 4
gamma = 7/5



#========================================#



if(Latex_interpreter):
    
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')


dx = (bx-ax)/Nx
dy = (by-ay)/Ny


exp_name = Prob_name+"_"+P_order+"_"+vis_name

data_file = "output_"+exp_name+"/"

if(Compute_bool):
    os.system("make -j"+str(thread_num))
    os.system("rm -r "+data_file)
    os.system("OMP_NUM_THREADS="+str(thread_num)+" ./run")
    os.system("mv output/ "+data_file)
    os.system("mkdir output")


### t
data_t = open(data_file+"t.out", 'r').read()
t = np.array(data_t.split()).astype(float)

Nt = t.shape[0]

print("Nt = ",Nt)

### u

data_u = open(data_file+"u.out", 'r').read()

u_flatten = data_u.split()

u = np.empty((Nt,Ny,Nx,System_dim))

for it in range(Nt):
    for i in range(Ny):
        for j in range(Nx):
            for k in range(System_dim):
            
                u[it,i,j,k] = float(u_flatten[it*(Nx*Ny)*(System_dim) + i*Nx*System_dim + j*System_dim + k])


u = np.where(np.isnan(u), np.nanmin(u), u)



data_x = open(data_file+"x.out", 'r').read()
x_flatten = np.array(data_x.split()).astype(float)
x = np.reshape(x_flatten, (Nx,-1))

order = x.shape[1]

data_y = open(data_file+"y.out", 'r').read()
y_flatten = np.array(data_y.split()).astype(float)
y = np.reshape(y_flatten, (Ny,-1))



x_ = np.linspace(ax+0.5*dx, bx-0.5*dx, Nx)
y_ = np.linspace(ay+0.5*dy, by-0.5*dy, Ny)


[X,Y] = np.meshgrid(x_,y_)

rho = u[:,:,:,0]
momentum1 = u[:,:,:,1]
momentum2 = u[:,:,:,2]
Energy = u[:,:,:,3]

Pressure = (gamma-1)*(Energy-0.5*((momentum1**2 + momentum2**2)/rho))

fig1, ax1 = plt.subplots()
#fig1, ax1 = plt.subplots(figsize=[3.4,3.4])
#fig2, ax2 = plt.subplots(figsize=[4.5, 3.4])
#fig3, ax3 = plt.subplots(figsize=[4.5, 3.4])
#fig4, ax4 = plt.subplots(figsize=[4.5, 3.4])


plt.figure(fig1.number)
Cplot1 = ax1.contour(X,Y,rho[0], 100, algorithm='threaded', linewidths=0.5)
ax1.set_aspect('equal')
#plt.tight_layout(rect=(0.0,0.0,1.0,0.975))
#cbar1 = fig1.colorbar(Cplot1)
#plt.xlabel("$x$")
#plt.ylabel("$y$")
plt.title("Density at $t = $"+str(t[0]))

plt.pause(0.001)

#plt.pause(7.5)
plt.pause(1.0)

for it in range(Nt):
    plt.figure(fig1.number)
    for coll in Cplot1.collections: 
        plt.gca().collections.remove(coll)
    #cbar1.remove()
    
    plt.figure(fig1.number)
    Cplot1 = ax1.contour(X,Y,rho[it], 100, algorithm='threaded', linewidths=0.5)
    #cbar1 = fig1.colorbar(Cplot1)
    #plt.xlabel("$x$")
    #plt.ylabel("$y$")
    plt.title("Density at $t = $"+str(t[it]))

    if(Nt>20):
        plt.pause(0.001)
    else:
        plt.pause(0.2)

    # v1, v2, E
    """
    plt.figure(fig2.number)
    Cplot2 = ax2.contourf(X,Y,momentum1[it], 100)
    cbar2 = fig2.colorbar(Cplot2)
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.title("Momentum in $x$ direction at $t = $"+str(t[it]))

    plt.pause(0.001)

    plt.figure(fig2.number)
    for coll in Cplot2.collections: 
        plt.gca().collections.remove(coll)
    cbar2.remove()

    
    plt.figure(fig3.number)
    Cplot3 = ax3.contourf(X,Y,momentum2[it], 100)
    cbar3 = fig3.colorbar(Cplot3)
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.title("Momentum in $y$ direction at $t = $"+str(t[it]))

    plt.pause(0.001)

    plt.figure(fig3.number)
    for coll in Cplot3.collections: 
        plt.gca().collections.remove(coll)
    cbar3.remove()

    
    plt.figure(fig4.number)
    Cplot4 = ax4.contourf(X,Y,Pressure[it], 100)
    cbar4 = fig4.colorbar(Cplot4)
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.title("Pressure at $t = $"+str(t[it]))
    
    plt.pause(0.001)

    plt.figure(fig4.number)
    for coll in Cplot4.collections: 
        plt.gca().collections.remove(coll)
    cbar4.remove()
"""


fig1.savefig(exp_name+"_rho"+".pdf", format="pdf")


### nu

data_nu = open(data_file+"nu.out", 'r').read()
nu_flatten = np.array(data_nu.split()).astype(float)
nu = np.reshape(nu_flatten, (Ny,Nx))

print("max nu: ", np.nanmax(nu))


fig,ax = plt.subplots(figsize=[4.5, 3.4])
plt.tight_layout(rect=(0.005,0.01,1.0,0.99))

plt.contourf(X,Y,nu, 50)
plt.contourf(X,Y,nu, 50)
cmap = mpl.cm.viridis
norm = mpl.colors.Normalize(vmin=np.nanmin(nu), vmax=np.nanmax(nu))
plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap))
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title("viscosity $\\nu$")

fig.savefig(exp_name+"_nu"+".pdf", format="pdf")


plt.show()

