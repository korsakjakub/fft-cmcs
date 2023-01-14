import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq, fftshift, ifft

n = 1000
L = 50
dx = L/n
x = np.arange(-L/2, L/2, dx, dtype="complex_")
k = (2*np.pi/L)*np.arange(-n/2, n/2)


def check_fft(f, f2, name):
    f = np.pad(f, (500,500))
    f2 = np.pad(f2, (500, 500))
    fhat = fftshift(fft(f))*dx
    k = np.linspace(-2*np.pi/L*n/2, 2*np.pi/L*n/2, len(fhat))
    plt.plot(k, abs(fhat), label="FFT")
    plt.plot(k, abs(f2), color='c', label="formula")
    plt.xlim([-20,20])
    plt.legend()
    plt.savefig(name)
    plt.close()

def plot_dfFFT(f, f2, name):
    fhat = fft(f)
    dfhat = k*fhat*(1j)    
    dfFFT = np.real(ifft(dfhat))

    plt.plot(x, dfFFT.real, label="FFT")
    plt.plot(x, f2, color='c', label="formula")
    plt.legend()
    plt.savefig(name)
    plt.close()

def task1():
    gauss = np.exp(-x ** 2)
    gauss_hat = np.exp(-0.25 * k**2)
    lorentzian = np.exp(-np.abs(x))
    lorentzian_hat = 1/(1+k**2)
    s = lambda v: 1.0 if -5 < v < 5 else 0.0
    step = np.array([s(v) for v in x])
    step_hat = 1/k * np.sin(5* k)
    c = lambda v: s(v) * np.cos(5.0 * v)
    conv = np.array([c(v) for v in x])
    conv_hat = (np.sin(5*(2-k))/(2-k) + np.sin(5*(2+k))/(2+k))/np.sqrt(2)
    check_fft(gauss, gauss_hat, "gauss.png")
    check_fft(lorentzian, lorentzian_hat, "lorentz.png")
    check_fft(step, step_hat, "step.png")
    check_fft(conv, conv_hat, "conv.png")

def task2():
    gauss = np.exp(-np.power(x, 2)/25)
    dgauss = - 2/25 * x * np.exp(-x**2 / 25)
    lorentzian = np.exp(-np.abs(x))
    dlorentzian = -np.exp(-abs(x))*np.sign(x)
    s = lambda v: 1.0 if -5 < v < 5 else 0.0
    step = np.array([s(v) for v in x])
    c = lambda v: s(v) * np.cos(2*v)
    conv = np.array([c(v) for v in x])
    plot_dfFFT(gauss, dgauss, 'dgauss.png')
    plot_dfFFT(lorentzian, dlorentzian, 'dlorentz.png')
    # plot_dfFFT(step, 'dstep.png')
    # plot_dfFFT(conv, 'dconv.png')

if __name__ == '__main__':
    task1()
    task2()
