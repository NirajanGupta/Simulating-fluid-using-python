import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

N = 51
L = 1.0
h = L / (N - 1)
dt = 0.0001
alpha = dt/(h**2)

#Initializing the problem
T = np.zeros((N, N))
T[N-1, :] = 1.0

T_new = np.zeros((N, N))
T_new[N-1, :] = 1.0

num_error = 1
epsilon = 1e-5
iterations = 0



while num_error > epsilon:
    for i in range(1,N-1):
        for j in range(1,N-1):
            T_new[i,j] = T[i,j] + alpha * (T[i+1,j] + T[i-1,j] + T[i,j+1] + T[i,j-1] - 4*T[i,j])
    

    #resetting the nume_error and calculating the new one
    num_error = 0
    for i in range(1, N-1):
        for j in range(1, N-1):
            num_error += np.abs(T_new[i,j] - T[i,j])


    #iteraton and advancement
    iterations += 1
    T = T_new.copy()

    # Plotting numerical error
    # if iterations % 1000 == 0:
    #     plt.figure(10)
    #     plt.semilogy(iterations, num_error, 'ko')
    #     plt.pause(0.01)


# Plotting the final temperature distribution
x_dom = np.arange(N)*h
y_dom =np.arange(N)*h
[X,Y] = np.meshgrid(x_dom, y_dom)

plt.figure(11)
plt.contourf(X, Y, T)
plt.grid(False)
plt.title('Temperature Distribution')
plt.xlabel('X')
plt.ylabel('Y')
plt.colorbar(label='Temperature')
plt.show()

    

 
    
    




