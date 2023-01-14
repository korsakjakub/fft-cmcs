import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq, fftshift, ifft

n = 2048
L = 10
dx = L / (n - 1)
k_max = 1 / (2 * dx)
x = np.linspace(-L / 2, L / 2, n)


def check_fft(f, f2, name):
    fhat = fftshift(fft(f))
    k = np.linspace(-k_max, k_max, len(fhat))
    plt.plot(k, abs(fhat * dx), label="FFT")
    plt.plot(k, abs(f2), color='c', label="formula")
    plt.xlim([-5, 5])
    plt.legend()
    plt.savefig(f"figures/{name}")
    plt.close()


def plot_dfFFT(f, f2=None, name=""):
    kp = fftshift((2 * np.pi / L) * np.arange(-n / 2, n / 2))
    fhat = fft(f)
    dfhat = kp * fhat * 1j
    dfFFT = ifft(dfhat).real

    plt.plot(x, dfFFT, label="FFT")
    if f2 is not None:
        plt.plot(x, f2, color='c', label="formula")
    plt.legend()
    plt.savefig(f"figures/{name}")
    plt.close()


def task1():
    k = np.linspace(-k_max, k_max, 2048)
    gauss = np.exp(-x ** 2)
    gauss_hat = np.sqrt(np.pi) * np.exp(-np.pi ** 2 * k ** 2)
    lorentzian = np.exp(-np.abs(x))
    lorentzian_hat = 2 / (4 * np.pi ** 2 * k ** 2 + 1)
    s = lambda v: 1.0 if -1 < v < 1 else 0.0
    step = np.array([s(v) for v in x])
    step_hat = 1j * np.exp(2j * np.pi * k) * (-1 + np.exp(-4 * 1j * np.pi * k)) / (2 * np.pi * k)
    c = lambda v: s(v) * np.cos(v)
    conv = np.array([c(v) for v in x])
    conv_hat = -1 / (2 * (4 * np.pi ** 2 * k ** 2 - 1)) * 1j * np.exp(-1j * (2 * np.pi * k + 1)) * (
            -2 * np.pi * k + np.exp(4j * np.pi * k + 2j) * (2 * np.pi * k - 1) + np.exp(4j * np.pi * k) * (
            2 * np.pi * k + 1) - np.exp(2j) * (2 * np.pi * k + 1) + 1)
    check_fft(gauss, gauss_hat, "gauss.png")
    check_fft(lorentzian, lorentzian_hat, "lorentz.png")
    check_fft(step, step_hat, "step.png")
    check_fft(conv, conv_hat, "conv.png")


def task2():
    gauss = np.exp(-np.power(x, 2))
    dgauss = - 2 * x * np.exp(-x ** 2)
    lorentzian = np.exp(-np.abs(x))
    dlorentzian = -np.exp(-abs(x)) * np.sign(x)
    s = lambda v: 1.0 if -1 < v < 1 else 0.0
    step = np.array([s(v) for v in x])
    c = lambda v: s(v) * np.cos(2 * v)
    conv = np.array([c(v) for v in x])
    dconv = -2 * step * np.sin(2 * x)
    plot_dfFFT(gauss, dgauss, 'dgauss.png')
    plot_dfFFT(lorentzian, dlorentzian, 'dlorentz.png')
    plot_dfFFT(step, name='dstep.png')
    plot_dfFFT(conv, dconv, 'dconv.png')


if __name__ == '__main__':
    task1()
    task2()
