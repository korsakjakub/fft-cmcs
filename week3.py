import os

import imageio
import matplotlib.pyplot as plt
import numpy as np

n = 2048
L = 40
dx = L / (n - 1)
z = np.linspace(-L / 2, L / 2, n)

if __name__ == '__main__':
    A0 = 1
    dt = 1e-1
    t0 = 1.0
    beta1 = 0.5
    beta2 = 0.0
    state = lambda t: A0 * t0 / np.sqrt(t0 ** 2 - 1j * beta2 * z) * np.exp(
        -0.5 * (t - beta1 * z) ** 2 / (t0 ** 2 - 1j * beta2 * z))
    fnames = []
    for i in range(100):
        plt.plot(z, abs(state(t0 + i * dt)), label=f'analytic A(z, t = {round(t0 + i * dt, 3)})')
        plt.xlim([-20, 20])
        plt.ylim([-0.25, 1.2])
        plt.legend()
        fnames.append(f'figures/pulse/{i}.png')
        plt.savefig(fnames[-1])
        plt.close()

    with imageio.get_writer('figures/pulse-analytic.gif', mode='I') as writer:
        for filename in fnames:
            image = imageio.v2.imread(filename)
            writer.append_data(image)

    for filename in set(fnames):
        os.remove(filename)
