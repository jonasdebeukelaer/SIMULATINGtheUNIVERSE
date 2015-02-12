import pyfftw
import math
import numpy as np
from scipy import signal
from PIL import Image
import time

def MakeSphere(shape, radius):

	xSize = shape[1]
	ySize = shape[0]
	zSize = shape[2] #CHECK THIS asdkflbsaflabsflkjsabdfjsadlhfbsajdhb

	xCentre = (xSize / 2.0) - 0.5
	yCentre = (ySize / 2.0) - 0.5
	zCentre = (zSize / 2.0) - 0.5

	radiusSquared = radius**2
	sphere = np.zeros((shape))

	for x in range(int(xCentre - radius), int(xCentre + radius)):
		for y in range(int(yCentre - radius), int(yCentre + radius)):
			for z in range(int(zCentre - radius), int(zCentre + radius)):
				rSquared = (x - xCentre)**2 + (y - yCentre)**2 + (z - zCentre)**2

				if rSquared < radiusSquared:
					sphere[x][y][z] = 0.1
	return sphere


#-----------------------------------------------------------------------------------------------------

def xyzGaussianArray(shape, xSigma, ySigma, zSigma):

	xSize = shape[1]
	ySize = shape[0]
	zSize = shape[2] #CHECK THIS asdkflbsaflabsflkjsabdfjsadlhfbsajdhb

	xCentre = (xSize / 2.0) - 0.5
	yCentre = (ySize / 2.0) - 0.5
	zCentre = (zSize / 2.0) - 0.5

	gaussianArray = np.zeros(shape)

	for x in range(0, xSize):
		for y in range(0, ySize):
			for z in range(0, zSize):
				xNumerator = (x - xCentre)**2
				yNumerator = (y - yCentre)**2
				zNumerator = (z - zCentre)**2

				xDenominator = 2.0 * xSigma**2
				yDenominator = 2.0 * ySigma**2
				zDenominator = 2.0 * zSigma**2

				exponent = (xNumerator / xDenominator) + (yNumerator / yDenominator) + (zNumerator / zDenominator)
				gaussianArray[x][y][z] = math.exp(-exponent)
	
	return gaussianArray

#------------------------------------------------------------------------------------------------------

volume = [100, 100, 100]
sphereRadius = 20

#sphere = MakeSphere(volume, sphereRadius)
sphere = np.zeros((volume))
sphere[50][20][30] = 10
gaussianSmoother = xyzGaussianArray(volume, 4, 4, 4)

sphereFFT     = pyfftw.builders.rfftn(sphere)
gaussianFFT = pyfftw.builders.rfftn(gaussianSmoother)

gaussianTransform = gaussianFFT()
sphereTransform = sphereFFT()

fourierSmoothedSphere = gaussianTransform * sphereTransform

inverseFFT = pyfftw.builders.irfftn(fourierSmoothedSphere)
inverseSphereFFT = inverseFFT()
inverseSphereFFT = np.fft.fftshift(inverseSphereFFT)

maxAmplitude = np.amax(inverseSphereFFT)
inverseSphereFFT *= (1.0 / maxAmplitude)

f = open('gaussianSmoothedResults/3d_smoothing_0.3D', 'w')
f.write('x\ty\tz\tAmplitude\n')
f.write('%i\t%i\t%i\t%f\n' % (volume[0], volume[1], volume[2], 1))
f.write('%i\t%i\t%i\t%f\n' % (0, 0, 0, 1))
for x in range(volume[0]/2, volume[0]/2 +1):
	for y in range(0, volume[1]):
		for z in range(0, volume[2]):
			if sphere[x][y][z] != 0:
				f.write('%i\t%i\t%i\t%f\n' % (x, y, z, sphere[x][y][z]))

for x in range(0, volume[0]):
	for y in range(0, volume[1]):
		for z in range(volume[2]/2, volume[2]/2+1):
			if sphere[x][y][z] != 0:
				f.write('%i\t%i\t%i\t%f\n' % (x, y, z, sphere[x][y][z]))

for x in range(0, volume[0]):
	for y in range(volume[1]/2, volume[1]/2+1):
		for z in range(0, volume[2]):
			if sphere[x][y][z] != 0:
				f.write('%i\t%i\t%i\t%f\n' % (x, y, z, sphere[x][y][z]))

f.close()

f = open('gaussianSmoothedResults/3d_smoothing_1.3D', 'w')
f.write('x\ty\tz\tAmplitude\n')
f.write('%i\t%i\t%i\t%f\n' % (volume[0], volume[1], volume[2], 1))
f.write('%i\t%i\t%i\t%f\n' % (0, 0, 0, 1))
for x in range(volume[0]/2, volume[0]/2 + 1):
	for y in range(0, volume[1]):
		for z in range(0, volume[2]):
			if gaussianSmoother[x][y][z]  > 0.01:
				f.write('%i\t%i\t%i\t%f\n' % (x, y, z, gaussianSmoother[x][y][z]))

for x in range(0, volume[0]):
	for y in range(0, volume[1]):
		for z in range(volume[2]/2, volume[2]/2+1):
			if gaussianSmoother[x][y][z]  > 0.01:
				f.write('%i\t%i\t%i\t%f\n' % (x, y, z, gaussianSmoother[x][y][z]))

for x in range(0, volume[0]):
	for y in range(volume[1]/2, volume[1]/2+1):
		for z in range(0, volume[2]):
			if gaussianSmoother[x][y][z]  > 0.01:
				f.write('%i\t%i\t%i\t%f\n' % (x, y, z, gaussianSmoother[x][y][z]))

f.close()

f = open('gaussianSmoothedResults/3d_smoothing_2.3D', 'w')
f.write('x\ty\tz\tAmplitude\n')
f.write('%i\t%i\t%i\t%f\n' % (volume[0], volume[1], volume[2], 1))
f.write('%i\t%i\t%i\t%f\n' % (0, 0, 0, 1))
for x in range(volume[0]/2, volume[0]/2 +1):
	for y in range(0, volume[1]):
		for z in range(0, volume[2]):
			if inverseSphereFFT[x][y][z] > 0.01:
				f.write('%i\t%i\t%i\t%f\n' % (x, y, z, inverseSphereFFT[x][y][z]))

for x in range(0, volume[0]):
	for y in range(0, volume[1]):
		for z in range(volume[2]/2, volume[2]/2+1):
			if inverseSphereFFT[x][y][z]  > 0.01:
				f.write('%i\t%i\t%i\t%f\n' % (x, y, z, inverseSphereFFT[x][y][z]))

for x in range(0, volume[0]):
	for y in range(volume[1]/2, volume[1]/2+1):
		for z in range(0, volume[2]):
			if inverseSphereFFT[x][y][z]  > 0.01:
				f.write('%i\t%i\t%i\t%f\n' % (x, y, z, inverseSphereFFT[x][y][z]))

f.close()


