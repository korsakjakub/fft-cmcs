import numpy as np
from numpy.fft import fft, fftshift
import matplotlib.pyplot as plt

def gauss(t, sigma):
    gauss = pow(2*sigma/np.pi, 1/4)*np.exp(-t**2*sigma)
    def gaussF_analitical(k):
        return pow(2*np.pi/sigma, 1/4)*np.exp(-np.pi**2*k**2/sigma)
    gaussD = -2 * pow(2*(sigma**5)/np.pi, 1/4) * t * np.exp(-sigma*t**2)
    return gauss, gaussF_analitical, gaussD

def step(t, dt, A):
    step =  np.piecewise(t, [np.abs(t) > dt/2, np.abs(t) <= dt/2], [0, A])
    def stepF(k):
        return A*np.sin(np.pi*k*dt)/(np.pi*k)
    return step, stepF

def pulse(t, dt, A, omega): 
    fp = np.piecewise(t, [np.abs(t) > dt/2, np.abs(t) <= dt/2], [0, A])*np.cos(omega*t)
    def fpF(k):
        return A*(np.sin(2*np.pi*(omega-k)*dt/2)/(omega-k)+np.sin(2*np.pi*(omega+k)*dt/2)/(omega+k))/(2*np.pi)
    return fp, fpF

def cos_sin_cos(t, w1, w2, w3):
    fun = np.cos(t*w1) + np.sin(t*w2) + np.cos(t*w3)
    
    funD = -w1* np.sin(t*w1) + w2* np.cos(t*w2) + w3 * -np.sin(t*w3)
    return fun, fun, funD 



def compare_analitical_with_fft(f):
    ft = np.pad(f[0], (512, 512), 'constant')
    fft_func = fftshift(fft(ft))
    k = np.linspace(-kmax, kmax, len(fft_func))
    f_analitic = f[1](k)
    
    plt.plot(k, abs(f_analitic), label='Analytic')
    plt.plot(k, abs(fft_func * dT), label='Fourier')
    plt.xlim(-5, 5)
    plt.legend()
    plt.savefig('test1.png')
    plt.close()

def derrivative(f):
    fh = np.fft.fft(f[0])
    kp = (2*np.pi/L) * np.arange(-N/2,N/2)
    kp = np.fft.fftshift(kp)
    dfhat = kp * fh * 1j
    df = np.real(np.fft.ifft(dfhat))

    if f[2] is not None:
        plt.plot(t,f[2])
    plt.plot(t, df)

    plt.savefig('test2.png')
    plt.close()


#main code
L, N = 10, 2048
dT = L / (N - 1)
t = np.linspace(-L/2, L/2, N)
kmax = 1 / (2 * dT)

derrivative(gauss(t, 2))
compare_analitical_with_fft(gauss(t, 2))