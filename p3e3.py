from matplotlib.pylab import *
import numpy as np
#------------------- grafico discretizacion-----------------------------------
L = 1.0
n = 100 
dx = L / n
x = linspace(0, L, n+1)
dt = 2.
Nt = 50000
u = zeros((Nt, n+1))
u[:,0] = 0. 
u[:,-1] = 20. 
u[0, 1:n] = 20
K = 79.5 
c = 450.  
rho = 7800. 
alpha = K*dt/(c*rho*dx**2)

for k in range(Nt-1):
	t = dt*k
	print(f"k = {k} t = {t}")
	u[k,0] = u[k,1] - 5 * dx
	for i in range(1,n):
		u[k+1, i] = u[k,i] + alpha*(u[k, i+1] - 2*u[k, i] + u[k, i-1])

	if k % 1000 == 0:
		plot(x,u[k,:])

title ("k = {}  t = {} s".format(k, k*dt))
show ()




#------------------------------replica 4.5----------------------------------

L = 1.04
n = 20  
dx = L / n
xu = 2
Nt = 100000
boizq = 0       
boder = 0       
boin = 20       
K = 0.001495
c = 1.023
rho = 2476
αf = K/(c * rho)
nodos = [0.104, 0.208, 0.416]
vdt = [1, 5, 10, 50, 100]

figure()

for x in nodos:
    for dt in vdt:
        α = K * dt / (c * rho * dx ** 2)
        Tem = Nt//dt 
        u = zeros((Tem, n + 1))
        u[:,0] = boizq
        u[:,-1] = boder
        u[0,1:n] = boin

        for k in range(Tem - 1):

            t = dt*k

            for i in range(1, n):
                u[k + 1, i] = u[k, i] + α * (u[k, i + 1] - 2 * u[k, i] + u[k, i - 1])

        vT = linspace(0, Nt, Tem)

        plot(vT, u[:, xu], label = f"Malla 20 delta t={dt} (s)")


    Fourier = []
    for time in range(Nt):
        sf = 0
        for ni in range(1, 50):
            sf += 40 * (1-(-1)**ni) / (ni*np.pi) * np.sin(ni*np.pi* x/L) * np.exp(-αf * (ni*np.pi/L)**2*time)
        Fourier.append(sf)


    vtf = linspace(0, Nt, Nt)

    plot(vtf, Fourier, label="Serie de Fourier", color="black", linestyle="--")

    title(f" replica 4.5 x = {x}")
    xticks([18000, 36000, 54000, 72000, 90000], ["5", "10", "15", "20", "25"])
    legend(loc="best")
    ylabel("Temperatura [C]")
    xlabel("Tiempo [horas]")

    xu *= 2

    show()


#--------------------------------------replica 4.7---------------------


L = 1.04
nn = [10, 20, 40, 60, 100] 
dt = 60
xu = [1,2,4,6,10,  2,4,8,12,20,  4,8,16,24,40]
Nt = 100000
boizq = 0       
boder = 0       
boin = 20       
K = 0.001495
c = 1.023
rho = 2476
αf = K/(c * rho)
it = 0
nodos = [0.104, 0.208, 0.416]
vdt = [1, 5, 10, 50, 100]

figure()

for x in nodos:
    for n in nn:
        dx = L/n
        α = K*dt/ (c * rho * dx ** 2)
        Tem = Nt // dt  
        u = zeros((Tem, n+1))
        u[:,0] = boizq
        u[:,-1] = boder
        u[0,1:n] = boin

        for k in range(Tem - 1):

            t = dt*k

            for i in range(1, n):
                u[k + 1, i] = u[k, i] + α * (u[k, i + 1] - 2 * u[k, i] + u[k, i - 1])

        vT = linspace(0, Nt, Tem)  

        plot(vT, u[ : ,xu[it]], label = f"M{n}: Malla {n}")

        it = it + 1


    Fourier = []
    for time in range(Nt):
        sf = 0
        for ni in range(1, 50):
            sf += 40 * (1-(-1)**ni) / (ni*np.pi) * np.sin(ni*np.pi* x/L) * np.exp(-αf * (ni*np.pi/L)**2*time)
        Fourier.append(sf)


    vtf = linspace(0, Nt, Nt)

    plot(vtf, Fourier, label="Serie de Fourier", color="black", linestyle="--")

    title(f"replica 4.7 x = {x}")
    xticks([18000, 36000, 54000, 72000, 90000], ["10", "20", "30", "40", "50"])
    legend()
    ylabel("Temperatura [C]")
    xlabel("Tiempo [horas]")

    show()







