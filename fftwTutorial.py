import pyfftw
import math
import numpy as np
from PIL import Image

def xyGaussianArray(shape, xSigma, ySigma):

	xSize = shape[1]
	ySize = shape[0]

	xCentre = (xSize / 2.0) - 0.5
	yCentre = (ySize / 2.0) - 0.5

	gaussianArray = np.zeros(shape)

	for x in range(0, xSize):
		for y in range(0, ySize):
			gaussianArray[x][y] = math.exp(-(((x-xCentre)**2/2.0*xSigma**2)+((y-yCentre)**2/2.0*ySigma**2)))
	
	print gaussianArray
	return gaussianArray

lena = Image.open("lena_bw.bmp")
lenaArray = np.array(lena)

gaussianSmoother = xyGaussianArray(lenaArray.shape, 0.00001, 0.00001)

lenaFFT = pyfftw.builders.fft(lenaArray)
lenaTransform = lenaFFT()

gaussianFFT = pyfftw.builders.fft(gaussianSmoother)
gaussianTransform = gaussianFFT()

fourierSmoothedLena = np.multiply(gaussianTransform, lenaTransform)

inverseLenaFFT = pyfftw.builders.ifft(fourierSmoothedLena)
smoothedLena = inverseLenaFFT()

smooshFaceLena = Image.fromarray(np.uint8(smoothedLena))
smooshFaceLena.save("Smoosh.bmp")