from cProfile import label
from matplotlib import scale
from pylab import *
import mpl_toolkits.mplot3d.axes3d as p3

Nx = 25
Ny = 25
radius = 8
Niter =  1500 
phi = zeros((Ny,Nx))
errors = [None]*Niter

n = arange(Niter)
x = linspace(-0.5, 0.5, Nx)
y = linspace(-0.5, 0.5, Ny)
Y,X = meshgrid(y,x)
ii = where((X*X + Y*Y) <= 0.35*0.35)
phi[ii] = 1

for k in n:
    oldphi = phi.copy()
    phi[1:-1,1:-1] = 0.25*(phi[1:-1,0:-2] + phi[1:-1,2:] + phi[0:-2,1:-1] + phi[2:,1:-1])
    
    phi[1:-1,0] = phi[1:-1,1]
    phi[1:-1,-1] = phi[1:-1,-2]
    phi[0,1:-1] = phi[1,1:-1]
    phi[ii] = 1

    errors[k] = (abs(phi-oldphi)).max()

best_est = lstsq(c_[n,ones(Niter)], log(errors), rcond=None)[0]
b = best_est[0]
a = best_est[1]

best_est_500 = lstsq(c_[n[500:],ones(Niter-500)], log(errors[500:]), rcond=None)[0]
b_500 = best_est_500[0]
a_500 = best_est_500[1]

Jx = zeros((Ny,Nx))
Jy = zeros((Nx,Ny))
Jx[1:-1,1:-1] = 0.5*(phi[1:-1,0:-2] - phi[1:-1,2:])
Jy[1:-1,1:-1] = 0.5*(phi[0:-2,1:-1] - phi[2:,1:-1])

figure(1)
plot(ii[0]/Nx-0.5, ii[1]/Ny-0.5, 'or', label="V = 1")
xlim(-0.5,0.5)
ylim(-0.5,0.5)
xlabel(r'$X\rightarrow$')
ylabel(r'$Y\rightarrow$')
grid(True)
legend()

figure(2)
semilogy(n,errors,label="error")
xlabel(r'$X\rightarrow$')
ylabel(r'$Errror\rightarrow$')
grid(True)
legend()

figure(3)
loglog(n,errors,label="error")
xlabel(r'$X\rightarrow$')
ylabel(r'$Error\rightarrow$')
grid(True)
legend()

figure(4)
semilogy(n[::50],errors[::50],'or', label="error (interval=50)")
semilogy(n,exp(a+(b*n)), '--b', label="fit1")
semilogy(n[500:],exp(a_500+(b_500*n[500:])),'--g', label="fit2")
xlabel(r'$X\rightarrow$')
ylabel(r'$Error\rightarrow$')
grid(True)
legend()

figure(5)
loglog(n[::50],errors[::50],'or',label="error (interval=50)")
loglog(n,exp(a+(b*n)), '--b', label="fit1")
loglog(n[500:],exp(a_500+(b_500*n[500:])),'--g', label="fit2")
xlabel(r'$X\rightarrow$')
ylabel(r'$Error\rightarrow$')
grid(True)
legend()

figure(6)
contourf(-Y,-X,phi)
plot(ii[0]/Nx-0.5, ii[1]/Ny-0.5, 'or', label="V = 1")
xlim(-0.5,0.5)
ylim(-0.5,0.5)
xlabel(r'$X\rightarrow$')
ylabel(r'$Y\rightarrow$')
grid(True)
legend()

fig6 = figure(7)
ax = p3.Axes3D(fig6)
title("The 3-D surface plot of the potential")
surf = ax.plot_surface(Y,X,phi.T, rstride=1, cstride=1, cmap=cm.jet)

figure(8)
plot(ii[0]/Nx-0.5, ii[1]/Ny-0.5, 'or', label="V = 1")
quiver(x,y,-Jx[::-1,:],-Jy[::-1,:],scale=4)
xlim(-0.5,0.5)
ylim(-0.5,0.5)
xlabel(r'$X\rightarrow$')
ylabel(r'$Y\rightarrow$')
grid(True)
legend()

show()

print(a,b)
print(a_500,b_500)