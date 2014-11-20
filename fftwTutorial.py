import pyfftw
import math
import numpy as np
from PIL import Image
from scipy import signal

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



lenaFFT = pyfftw.builders.rfft2(lenaArray)
lenaTransform = lenaFFT()

gaussianFFT = pyfftw.builders.rfft2(gaussianSmoother)
gaussianTransform = gaussianFFT()

fourierSmoothedLena = lenaTransform * gaussianTransform
#fourierSmoothedLena = np.zeros(shape=(gaussianTransform.shape))

#Center is at 256 128
fourierSmoothedLena[255][127] = 1.0
#fourierSmoothedLena[255][129] = 1.0
#fourierSmoothedLena[257][127] = 1.0
#fourierSmoothedLena[257][129] = 1.0

inverseLenaFFT = pyfftw.builders.irfft2(fourierSmoothedLena)
smoothedLena = inverseLenaFFT()

maxVal = np.amax(smoothedLena)
smoothedLena *= (256 / maxVal)


smooshFaceLena = Image.fromarray(np.uint8(smoothedLena))
print smooshFaceLena

smooshFaceLena.save("Smoosh.bmp")