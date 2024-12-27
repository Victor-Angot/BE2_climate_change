import matplotlib.pyplot as plt
import numpy as np
from math import exp

N = 10
T_N_hom = 255
T_N_inhom = 214
H = 2
verif_model = False

def eps(z):
    if verif_model:
        return 1
    else:
        return exp(-z/H)

def compute_Tn_hom(n, T_N):
    return T_N*(N+1-n)**(1/4)

def compute_Tn_inhom(n, T_N):
    if n > N:
        return None
    epsilon = [eps(z) for z in range(0, N+1)]
    epsilon_sum = 0
    for i in range(n+1, N+1):
        epsilon_sum += epsilon[i]/(2 - epsilon[i])
    
    T_n = T_N * ((2 - epsilon[N]) * (1/(2-epsilon[n]) + epsilon_sum))**(1/4)
    return T_n

def load_data_from_txt(file_path):
    z_km = []  
    Ta_K = [] 
    with open(file_path, 'r') as file:
        next(file)
        for line in file:
            z, Ta = line.split()
            z_km.append(float(z))
            Ta_K.append(float(Ta)) 
    return z_km, Ta_K

if __name__ == '__main__':

    file_path = 'data_headers.txt'  
    z_km, Ta_K = load_data_from_txt(file_path)

    z = [i for i in range(0, N+1)]
    T_hom = [compute_Tn_hom(n, T_N_hom) for n in z]
    T_inhom = [compute_Tn_inhom(n, T_N_inhom) for n in z]

    if not verif_model:
        plt.plot(z_km, Ta_K, 'ro')
    plt.plot(z, T_hom, 'go')
    if verif_model:
        plt.plot(z, T_inhom, 'bo', marker='x')
    else:
        plt.plot(z, T_inhom, 'bo')
    if verif_model:
        plt.legend(['Homogeneous model', 'Inhomogeneous model ($\epsilon$ = 1)'])
    else:
        plt.legend(['Experimental data', 'Homogeneous model ($T_N$ = 255K)', 'Inhomogeneous model ($T_N$ = 214K)'])
    plt.xlabel('Elevation (km)')
    plt.ylabel('Temperature (K)')
    plt.title('Atmosphere temperature evolution with elevation')

    plt.show()
