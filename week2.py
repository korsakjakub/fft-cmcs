import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq, fftshift, ifft
import scipy.linalg
import imageio.v2
import imageio
import os

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
    dt = 3e-1
    state = (lambda t: A0 * np.exp(-t ** 2 / t0 ** 2))(x)
    plt.plot(x, state)
    fnames = []
    for i in range(150):
        state = propagate(state, dt)
        scipy.linalg.norm(state)

        plt.plot(x, abs(state), label=f'|ψ|(t = {round(i*dt, 3)})')
        plt.plot(x, state.real, label=f'real part')
        plt.plot(x, state.imag, label=f'imaginary part')
        plt.xlim([-10, 10])
        plt.ylim([-0.3, 1.1])
        plt.legend()
        fnames.append(f'figures/splitstep/{i}.png')
        plt.savefig(fnames[-1])
        plt.close()

    with imageio.get_writer('figures/splitstep.gif', mode='I') as writer:
        for filename in fnames:
            image = imageio.v2.imread(filename)
            writer.append_data(image)

    for filename in set(fnames):
        os.remove(filename)