import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq, fftshift, ifft
import scipy.linalg

n = 10000
L = 100
dx = L/n
x = np.arange(-L/2, L/2, dx, dtype="complex_")
k = (2*np.pi/L)*np.arange(-n/2, n/2)
k = fftshift(k)

def propagate(psi, dt):
    psi_hat = fft(psi)
    return ifft(psi_hat*np.exp(1j*dt*k**2))

if __name__ == '__main__':
    A0 = 1
    t0 = 10.0
    dt = 1.0
    initial_state = lambda t: A0 * np.exp(-0.5*t**2/t0**2)

    state = initial_state(x)
    for i in range(100):
        state = propagate(state, dt)
        scipy.linalg.norm(state)

        plt.plot(x, state)
        plt.xlim([-50,50])
    plt.savefig('splitstep.png')