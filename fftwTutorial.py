import pyfftw
import math
import numpy as np
from scipy import signal
from PIL import Image
import time

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

gaussianSmoother = xyGaussianArray(lenaArray.shape, 20, 2)

lenaFFT     = pyfftw.builders.rfftn(lenaArray)
gaussianFFT = pyfftw.builders.rfftn(gaussianSmoother)

gaussianTransform = gaussianFFT()
lenaTransform = lenaFFT()

fourierSmoothedLena = gaussianTransform * lenaTransform

fourierSmoothedLena = np.full(lenaTransform.shape, 0)
fourierSmoothedLena[0][1] = 200
print fourierSmoothedLena

inverseFFT = pyfftw.builders.ifftn(fourierSmoothedLena)
inverseLenaFFT = inverseFFT()
#inverseLenaFFT = np.fft.fftshift(inverseLenaFFT)

maxAmplitude = np.amax(inverseLenaFFT)
inverseLenaFFT *= (256.0 / maxAmplitude)

smooshFaceLena = Image.fromarray(inverseLenaFFT.astype('uint8'))
smooshFaceLena.save("pic_archive/Smoosh_%s.bmp" % (time.strftime("%c")))
smooshFaceLena.save("Smoosh.bmp")


