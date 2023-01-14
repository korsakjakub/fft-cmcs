import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq, fftshift, ifft
import scipy.linalg

n = 2048
L = 40
dx = L / (n - 1)
x = np.linspace(-L / 2, L / 2, n)


def propagate(psi, dt):
    k_max = 1 / (2 * dx)
    k = fftshift(np.linspace(-k_max, k_max, n))
    psi_hat = fft(psi)
    return ifft(psi_hat * np.exp(1j * dt * k ** 2))


if __name__ == '__main__':
    A0 = 1
    t0 = 1
    dt = 1e-1
    state = (lambda t: A0 * np.exp(-t ** 2 / t0 ** 2))(x)
    plt.plot(x, state)
    for i in range(300):
        state = propagate(state, dt)
        scipy.linalg.norm(state)

        plt.plot(x, state)
        plt.xlim([-10, 10])
    plt.savefig('figures/splitstep.png')
