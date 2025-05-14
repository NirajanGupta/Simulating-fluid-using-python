# “All that had gone before was not a thousandth of what was yet to come; the story of this star had barely begun.”

import numpy as np
import scipy.optimize as optimize
import matplotlib.pyplot as plt

g = 1.4 # specific heat ratio (gamma)
M = 5.7 # upstream mach number
theta_1 = 20 # negative angle
theta_2 = 35

# Possible range of beta and theta for M | Baseline plot | locus of all the possible shock wave formations for M
beta = np.arange(np.arcsin(1/M), np.pi/2, np.pi/10000)
theta = np.rad2deg(np.arctan(2*(M**2*np.sin(beta)**2-1)/(np.tan(beta)*(M**2*(g+np.cos(2*beta))+2))))
PR = 1 + (2*g)/(g+1)*(M**2*np.sin(beta)**2-1)
plt.plot(theta, PR, 'k--', label='baseline' )
plt.plot(-theta, PR, 'k--', 0, 1, "o", markerfacecolor='k', linewidth=1.5)

#Function to find beta for a given theta
def calculate_beta(theta_val):
    theta_function = lambda beta_val, theta_val:np.arctan(2*(M**2*np.sin(beta_val)**2-1)/(np.tan(beta_val)*(M**2*(g+np.cos(2*beta_val))+2)))*180/np.pi - theta_val
    beta_initial_guess = 0
    beta_solution = optimize.fsolve(theta_function, beta_initial_guess, args=(theta_val))
    return beta_solution[0]

#Calculation of flow and shock wave properties for theta_1
beta_1 = calculate_beta(theta_1)
print(beta_1)
print(np.rad2deg(beta_1))

Mn1_1 = M * np.sin(beta_1)
Mn2_1 = np.sqrt((1+(g-1)/2*Mn1_1**2)/(g*Mn1_1**2-(g-1)/2))
M2_1 = Mn2_1/np.sin(beta_1-np.deg2rad(theta_1))
PR_1 = 1 + (2*g)/(g+1)*(M**2*np.sin(beta_1)**2-1)

beta1 = np.arange(np.arcsin(1/M2_1), np.pi/2, np.pi/10000)
theta1 = np.rad2deg(np.arctan(2*(M2_1**2*np.sin(beta1)**2-1)/(np.tan(beta1)*(M2_1**2*(g+np.cos(2*beta1))+2))))
PR1 = (1 + (2*g)/(g+1)*(M2_1**2*np.sin(beta1)**2-1))*PR_1
plt.plot(-(theta_1+theta1), PR1, 'b', label = '$θ_{1}$')
plt.plot(-(theta_1-theta1), PR1, 'b',-theta_1, PR_1, "o", markerfacecolor='b', linewidth=1.5)

#Calculation of flow and shock wave properties for theta_2
beta_2 = calculate_beta(theta_2)
print(np.rad2deg(beta_2))

Mn1_2 = M * np.sin(beta_2)
Mn2_2 = np.sqrt((1+(g-1)/2*Mn1_2**2)/(g*Mn1_2**2-(g-1)/2))
M2_2 = Mn2_2/np.sin(beta_2-np.deg2rad(theta_2))
PR_2 = 1 + (2*g)/(g+1)*(M**2*np.sin(beta_2)**2-1)

beta2 = np.arange(np.arcsin(1/M2_2), np.pi/2, np.pi/10000)
theta2 = np.rad2deg(np.arctan(2*(M2_2**2*np.sin(beta2)**2-1)/(np.tan(beta2)*(M2_2**2*(g+np.cos(2*beta2))+2))))
PR2 = (1 + (2*g)/(g+1)*(M2_2**2*np.sin(beta2)**2-1))*PR_2
plt.plot(theta_2+theta2, PR2, 'r',label = '$θ_{2}$')
plt.plot(theta_2-theta2, PR2, 'r', theta_2, PR_2, "o", markerfacecolor='r', linewidth=1.5)

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
plt.title('Shock Polar Plot', fontdict = font_propt)
plt.xlabel(' Deflection angle, θ (in degree)', fontdict = font_propl)
plt.ylabel(' Pressure ratio, $P_{2}$/$P_{1}}$', fontdict = font_propl)
plt.xlim([-60, 60])
plt.ylim([0, 150])
plt.grid(True)
plt.legend()
plt.show()