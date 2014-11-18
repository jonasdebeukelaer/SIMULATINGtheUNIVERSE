import pyfftw
import numpy

a = pyfftw.n_byte_align_empty(128, 16, 'complex128')
b = pyfftw.n_byte_align_empty(128, 16, 'complex128')
c = pyfftw.n_byte_align_empty(128, 16, 'complex128')

fft_object = pyfftw.FFTW(a, b)
ifft_object = pyfftw.FFTW(b, c, direction='FFTW_BACKWARD')

ar, ai = numpy.random.randn(2, 128)
a[:] = ar + 1j*ai

fft_a = fft_object()
ifft_b = ifft_object()