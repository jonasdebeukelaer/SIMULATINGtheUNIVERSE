import pyfftw
import numpy as np
from PIL import Image

lena = Image.open("lena_bw.bmp")
lenaArray = np.array(lena)

for i in range(0, 10000000):
	lenaFFT = pyfftw.builders.fft(lenaArray)
	transformLenaArray = lenaFFT()

	inverseLenaFFT = pyfftw.builders.ifft(transformLenaArray)
	lenaArray = inverseLenaFFT()

smooshFaceLena = Image.fromarray(np.uint8(lenaArray))
smooshFaceLena.save("Smoosh.bmp")