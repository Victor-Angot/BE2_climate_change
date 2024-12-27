import matplotlib.pyplot as plt

T_N = 220
N = 10

def compute_T(n):
    return T_N*((N+1-n)**(1/4))

n = [i for i in range(0, N+1)]
T = [compute_T(i) for i in n]

plt.plot(n, T, 'bo')
plt.xlabel('Layer number')
plt.ylabel('Temperature (K)')
plt.title('Atmosphere temperature evolution in the layers')

plt.show()

