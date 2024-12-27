from Q3 import compute_Tn_inhom
from math import *
import matplotlib.pyplot as plt
import numpy as np

scenario_epsilon_50percentup = False
H = 2
sigma = 5.67 * 10 ** -8
N = 10 
time_scale = 25 # number of years
delta_t_f = time_scale * 365 * 24 * 3600
#5 * 10 ** 9  # final time in seconds
delta_t = 36000  # 10 hours in seconds
layer_thickness = 1  
albedo = 0.3  
C = 10 ** 7
Es = 1370 / 4 
Es_prim = Es * (1 - albedo)

T_i = 500
T_initial = (N + 1) * [T_i]

def increase_eps(eps):
    eps *= 1.5
    if eps >= 1:
        eps = 1
    return eps

def upward_energy_flux(N, layer_thickness, albedo, C, T):
    """Create the list of upward energy fluxes from layers 0 to N at time t (corresponding to T)"""
    epsilon_0 = exp(-0 * layer_thickness / H)
    if scenario_epsilon_50percentup:
       epsilon_0 = increase_eps(epsilon_0)

    upward_flux = [epsilon_0 * sigma * (T[0] ** 4)]  # upward flux for layer 0

    for n in range(1, N + 1):
        epsilon_n = exp(-(n * layer_thickness) / H)  # calculate epsilon for layer n
        if scenario_epsilon_50percentup:
            epsilon_n = increase_eps(epsilon_n) # Q4.2 assumption
        upward_flux.append(epsilon_n * sigma * (T[n] ** 4) + (1 - epsilon_n) * upward_flux[n - 1])

    return upward_flux

def downward_energy_flux(N, layer_thickness, albedo, C, T):
    """Create the list of downward energy fluxes from layers N to 0 at time t (corresponding to T)"""
    epsilon_N = exp(-(N * layer_thickness) / H)  # calculate epsilon for layer N
    if scenario_epsilon_50percentup:
        epsilon_N = increase_eps(epsilon_N) # Q4.2 assumption

    downward_flux = [epsilon_N * sigma * (T[N] ** 4)]

    for n in range(1, N):
        epsilon = exp(-((N - n) * layer_thickness) / H)  # calculate epsilon for layer n
        if scenario_epsilon_50percentup:
            epsilon = increase_eps(epsilon) # Q4.2 assumption
        downward_flux.append(epsilon * sigma * (T[N - n] ** 4) + (1 - epsilon) * downward_flux[n - 1])

    downward_flux.append(0)  # add zero flux for the ground layer
    downward_flux = np.flip(downward_flux)  # reverse the order of the list

    return downward_flux

def next_temperature(delta_t, N, layer_thickness, albedo, C, T):
    """Calculate the list of temperatures delta_t after time t corresponding to the list of temperatures T"""
    upward_flux_list = upward_energy_flux(N, layer_thickness, albedo, C, T)  # calculate upward fluxes for layers 0 to N
    downward_flux_list = downward_energy_flux(N, layer_thickness, albedo, C, T)  # calculate downward fluxes for layers 0 to N

    next_temp_list = [T[0] + ((Es_prim + downward_flux_list[1] - upward_flux_list[0]) / C) * delta_t]

    for n in range(1, N):
        next_temp_list.append(T[n] + ((downward_flux_list[n + 1] + upward_flux_list[n - 1] - upward_flux_list[n] - downward_flux_list[n]) / C) * delta_t)

    next_temp_list.append(T[N] + ((upward_flux_list[N - 1] - upward_flux_list[N] - downward_flux_list[N]) / C) * delta_t)

    return next_temp_list


def temperature_matrix(delta_t, delta_t_f, N, layer_thickness, albedo, C, T_initial):
    """Build the matrix of temperatures of layers 0 to N over time (from t_0 to t_0 + delta_t_f in steps of delta_t)"""
    temp_matrix = []  # create the temperature matrix with time as x-axis and layer as y-axis

    for n in range(0, N + 1):
        temp_matrix.append([T_initial[n]])

    num_time_steps = floor(delta_t_f / delta_t)
    current_temp_list = T_initial

    for _ in range(1, num_time_steps + 1):
        current_temp_list = next_temperature(delta_t, N, layer_thickness, albedo, C, current_temp_list)
        for n in range(N + 1):
            temp_matrix[n].append(current_temp_list[n])

    return temp_matrix

temp_matrix_1 = temperature_matrix(delta_t, delta_t_f, N, layer_thickness, albedo, C, T_initial)

plt.figure()
for i in range(N + 1):
    plt.plot(range(len(temp_matrix_1[0])), temp_matrix_1[i], label="Layer {0}".format(i))

plt.legend(title="Legend", loc='upper right')
plt.xlabel("Iterations")
plt.ylabel("Temperature (K)")
plt.title("Atmosphere layers temperature evolution over {} years ($T_i$ = {}K, $\Delta_t$ = {}h)".format(time_scale, T_i, delta_t/3600))
plt.show()

# List of layer temperatures in the steady state
stable_temperatures = [temp_matrix_1[layer][-1] for layer in range(len(temp_matrix_1))]

T_N = 220
Temp_inhomo2 = [compute_Tn_inhom(i, T_N) for i in range(0, N)] + [T_N]

plt.figure()
plt.plot(range(N + 1), stable_temperatures, label="Non-stationary")
plt.plot(range(N + 1), Temp_inhomo2, label="Stationary")
plt.xlabel("Layer")
plt.ylabel("Temperature (K)")
plt.legend()
plt.title("Atmosphere layers temperature in the stable state")
plt.show()
