import pyfftw
import math
import numpy as np
from scipy import signal
from PIL import Image

def xyGaussianArray(shape, xSigma, ySigma):

	xSize = shape[1]
	ySize = shape[0]

	xCentre = (xSize / 2.0) - 0.5
	yCentre = (ySize / 2.0) - 0.5

	gaussianArray = np.zeros(shape)

	for x in range(0, xSize):
		for y in range(0, ySize):
			xNumerator = (x - xCentre)**2
			yNumerator = (y - yCentre)**2

			xDenominator = 2.0 * xSigma**2
			yDenominator = 2.0 * ySigma**2

			exponent = (xNumerator / xDenominator) + (yNumerator / yDenominator)
			gaussianArray[x][y] = math.exp(-exponent)
	
	return gaussianArray

lena = Image.open("lena_bw.bmp")
lenaArray = np.array(lena)

gaussianSmoother = xyGaussianArray(lenaArray.shape, 5, 5)

lenaFFT     = pyfftw.builders.rfftn(lenaArray)
gaussianFFT = pyfftw.builders.rfftn(gaussianSmoother)

gaussianTransform = gaussianFFT()
lenaTransform = lenaFFT()
print gaussianTransform.shape
print lenaTransform.shape
fourierSmoothedLena = (gaussianTransform.real * lenaTransform.real) + (gaussianTransform.imag * lenaTransform.imag) * 1j

inverseFFT = pyfftw.builders.irfftn((fourierSmoothedLena))
inverseLenaFFT = inverseFFT()

maxAmplitude = np.amax(inverseLenaFFT)
inverseLenaFFT *= (256.0 / maxAmplitude)

smooshFaceLena = Image.fromarray(np.uint8(inverseLenaFFT))
smooshFaceLena.save("Smoosh.bmp")


