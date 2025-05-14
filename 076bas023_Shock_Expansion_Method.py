# Now I am become Death, the destroyer of worlds.#

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize

#Tangent Ogive Profile Parameters
rho = 1.5 # radius of the circle (in m)
L = 1 # length of the ogive profile (in m)
R = rho - np.sqrt(rho**2 - L**2)#base radius of the tangent ogive profile
print(R)

#profile independent data calculation
x = np.linspace(0, L, 1000)
y = np.sqrt(rho**2 - (L-x)**2) + R-rho

#Freestream Flow Parameters
M1 = 6.2 #Mach number
g = 1.4 #specific heat ratio
#atmospheric properties corresponding to alt = 30km (from sea level)
T1 = 226.6 #Temperature in K (in celcius -46.64 )
P1 = 1197 #Pressure in Pa
rhoa1 = 0.01841 #density in kg/m^3
a = np.sqrt((g*P1)/rhoa1) #speed of sound
u1 = M1 * a #flow velocity 

#Function to find the inclination angle(theta/deflection angle) at any point on the profile
def theta_calculation(x):
    theta = np.arctan((L - x) / np.sqrt(rho**2 - (L - x)**2))
    return theta

#Function to find beta for a given theta
def beta_calculation(theta_val):
    theta_function = lambda beta_val, theta_val : np.arctan(2*(M1**2*np.sin(beta_val)**2-1)/(np.tan(beta_val)*(M1**2*(g+np.cos(2*beta_val))+2))) - theta_val
    beta_initial_guess = 0.3
    beta_solution = optimize.fsolve(theta_function, beta_initial_guess, args=(theta_val), xtol = 1e-6)
    return beta_solution[0]

#Function to find the expansion angle (Prandtl_Meyer_Function)
def expansion_angle(del_theta):
    v = np.sqrt((g+1)/(g-1)) * np.arctan(np.sqrt((g-1)/(g+1) * (M[0]**2 - 1))) - np.arctan(np.sqrt(M[0]**2 - 1))
    del_theta_function = lambda Mi, del_theta :  np.sqrt((g+1)/(g-1))*np.arctan(np.sqrt((g-1)/(g+1) * (Mi**2 - 1))) - np.arctan(np.sqrt(Mi**2 - 1)) - v - del_theta
    Mi_initial_guess = 1.5
    Mi_solution = optimize.fsolve(del_theta_function, Mi_initial_guess, args = (del_theta), xtol = 1e-6)
    return Mi_solution[0]

#Flow properties due to the oblique shockwave at nose and the shock parameters
theta_n = theta_calculation(x[0])  #deflection at nose
#print(theta_n, np.rad2deg(theta_n))
beta_n = beta_calculation(theta_n) # shockwave angle
#print(beta_n, np.rad2deg(beta_n))

Mn1 = M1 * np.sin(beta_n) # Mach normal component before shock
Mn2 = np.sqrt((1+ (g-1)/2 * Mn1**2)/(g*Mn1**2-(g-1)/2))# Mach normal component behind the shock
M2_n = Mn2/np.sin(beta_n - theta_n)# Flow Mach Behind the shock
P2_n = (1 + (2*g)/(g+1)*(M1**2*np.sin(beta_n)**2-1)) * P1 # Pressure behind the shock
Cp_n = (P2_n - P1)/ (0.5 * rhoa1 * u1**2) #coeff of pressure at nose
#print(Mn1, Mn2, Cp_n, M2_n, P2_n)


#Cp calculation along the profile
P = np.empty(len(x))
M = np.empty(len(x))
Cp = np.empty(len(x))
M[0] = M2_n # Mn
P[0] = P2_n # Pn
Cp[0] = Cp_n
for i in range(1, len(x)-1):
    del_theta1 = theta_n - theta_calculation(x[i])
    M[i] = expansion_angle(del_theta1)
    P[i] = (((1 + ((g-1)/2) * M[0]**2)/(1 + ((g-1)/2) * M[i]**2))**(g/(g-1))) * P[0]
    Cp[i] = (P[i] - P1)/ (0.5 * rhoa1 * u1**2)

Cp[len(x)-1] = 2*Cp[len(x)-2] - Cp[len(x)-3]
#print(M, Cp)

plt.plot(x/L, Cp, 'k', label = 'Cp on the upper profile', linewidth = 1.5 )

# Set font properties for labels
font_propl = {
    'family': 'serif', 
    'size': 9, 
    'weight': 'bold'  
}
# Set font properties for title
font_propt = {
    'family': 'serif', 
    'size': 14,  
    'weight': 'bold'
}

# Set plot title and labels
plt.title('Cp vs x/L', fontdict = font_propt)
plt.xlabel(' Distance ratio along symmetric axis, x/L', fontdict = font_propl)
plt.ylabel(' Coefficient of Pressure, $C_{p}$', fontdict = font_propl)
plt.xlim([0, 1.2])
plt.ylim([0, 1.3])
plt.grid(True)
plt.legend()
plt.show()


