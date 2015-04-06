import pmHelpers as helpers

import numpy as np
import math
import sys
import pyfftw
import matplotlib.pyplot as plt

def FindMeshIndex(position, nGrid):
	index = math.floor(position + 0.5) + ((nGrid / 2) - 1)
	if index == -1:
		index = nGrid - 1
	return index

def CalculateDensityField(nGrid, particleList):
	densityFieldMesh = np.empty([nGrid, nGrid, nGrid])
	densityFieldMesh.fill(-1)

	for particle in particleList:
		xMesh = FindMeshIndex(particle.position[0], nGrid)
		yMesh = FindMeshIndex(particle.position[1], nGrid)
		zMesh = FindMeshIndex(particle.position[2], nGrid)

		densityFieldMesh[xMesh][yMesh][zMesh] += 1

	return densityFieldMesh

def CreateGreensFunction(nGrid):
	shape       = (nGrid, nGrid, nGrid / 2 + 1)	
	greensArray = np.zeros((shape))

	constant = 1

	for l in range(0, shape[0]):
		if l < (shape[0] / 2):
			kx = 2 * math.pi * l / (shape[0])
		else:
			kx = 2 * math.pi * (l - shape[0]) / (shape[0])

		for m in range(0, shape[1]):
			if m < (shape[1] / 2):
				ky = 2 * math.pi * m / (shape[1])
			else:
				ky = 2 * math.pi * (m - shape[1]) / (shape[1])

			for n in range(0, shape[2]):
				kz = math.pi * n / (shape[2])

				if l != 0 or m != 0 or n != 0:
					greensArray[l][m][n] = - constant / ((math.sin(kx * 0.5))**2 + (math.sin(ky * 0.5))**2 + (math.sin(kz * 0.5))**2)

	return greensArray
		
def InlineGreensConvolution(densityFFT, a):
	shape = densityFFT.shape
	convolutedDensity = np.zeros((shape), dtype='complex128')

	constant = 3 / (8 * a)

	for l in range(0, shape[0]):
		if l < (shape[0] / 2):
			kx = 2 * math.pi * l / (shape[0])
		else:
			kx = 2 * math.pi * (l - shape[0]) / (shape[0])

		for m in range(0, shape[1]):
			if m < (shape[1] / 2):
				ky = 2 * math.pi * m / (shape[1])
			else:
				ky = 2 * math.pi * (m - shape[1]) / (shape[1])

			for n in range(0, shape[2]):
				kz = math.pi * n / (shape[2])

				if l != 0 or m != 0 or n != 0:
					greensValue = - constant / ((math.sin(kx * 0.5))**2 + (math.sin(ky * 0.5))**2 + (math.sin(kz * 0.5))**2)
					convolutedDensity[l][m][n] = greensValue * densityFFT[l][m][n]

	return convolutedDensity

def SolvePotential(densityField, a, greensFunction, preComputeGreens):
	densityFieldFFT = pyfftw.builders.rfftn(densityField, threads = helpers.GetNumberOfThreads())
	densityFFT = densityFieldFFT()

	if not preComputeGreens:
		densityFieldConvoluted = InlineGreensConvolution(densityFFT, a)
	else:
		scaledGreensFunction = np.multiply(3 / (8 * a), greensFunction)
		densityFieldConvoluted = np.multiply(scaledGreensFunction, densityFFT)

	potentialReal = pyfftw.builders.irfftn(densityFieldConvoluted, threads = helpers.GetNumberOfThreads())
	potentialField = potentialReal()
	#potentialField = np.fft.irfftn(densityFieldConvoluted)

	return potentialField

def FindPlusMinus(meshIndex, nGrid):
	if meshIndex == 0:
		meshPlus  = meshIndex + 1
		meshMinus = nGrid - 1
	elif meshIndex == nGrid - 1:
		meshPlus  = 0
		meshMinus = meshIndex - 1
	else:
		meshPlus  = meshIndex + 1
		meshMinus = meshIndex - 1

	return (meshPlus, meshMinus)

def CalculateParticleAcceleration(particle, potentialField):
	meshShape = potentialField.shape

	xMesh = FindMeshIndex(particle.position[0], meshShape[0])
	yMesh = FindMeshIndex(particle.position[1], meshShape[1])
	zMesh = FindMeshIndex(particle.position[2], meshShape[2])

	xNeighbours = FindPlusMinus(xMesh, meshShape[0])
	yNeighbours = FindPlusMinus(yMesh, meshShape[1])
	zNeighbours = FindPlusMinus(zMesh, meshShape[2])

	xAcceleration = - (potentialField[xNeighbours[0]][yMesh][zMesh] - potentialField[xNeighbours[1]][yMesh][zMesh]) / 2.0
	yAcceleration = - (potentialField[xMesh][yNeighbours[0]][zMesh] - potentialField[xMesh][yNeighbours[1]][zMesh]) / 2.0
	zAcceleration = - (potentialField[xMesh][yMesh][zNeighbours[0]] - potentialField[xMesh][yMesh][zNeighbours[1]]) / 2.0

	return (xAcceleration, yAcceleration, zAcceleration)

def CalculateParticleAccelerationPP(particle, i, particleList, numParticles):
	for otherIndex in range(i+1, numParticles):
		otherParticle = particleList[otherIndex]
		force = CalcForce(particle.position, otherParticle.position)
		particle.acceleration += force
		otherParticle.acceleration -= force

def CalcForce(pos1, pos2):
	epsilon = 1
	separation = np.subtract(pos2, pos1)
	separationMagnitudeSquared = sum(separation**2) + epsilon**2
	force = 1 / separationMagnitudeSquared**(1.5)
	force = np.multiply(separation, force)
	return force

def GetF(a, omega_m = 1, omega_k = 0, omega_lambda = 0):
	return (a**(-1)*(omega_m + omega_k * a + omega_lambda * a**3))**(-0.5)

def PositionCorrect(particle, nGrid):
	positionLimit = nGrid / 2
	for index, position in enumerate(particle.position):
		if position > positionLimit:
			while position > positionLimit:
				position = position - nGrid
			particle.position[index] = position
		elif position < (- positionLimit):
			while position < (- positionLimit):
				position = position + nGrid
			particle.position[index] = position

def FindLocalDensity(particle, densityField):
	meshShape = densityField.shape

	xMesh = FindMeshIndex(particle.position[0], meshShape[0])
	yMesh = FindMeshIndex(particle.position[1], meshShape[1])
	zMesh = FindMeshIndex(particle.position[2], meshShape[2])

	return densityField[xMesh][yMesh][zMesh]

def CalculatePowerSpectrum(densityField, nGrid, lBox, dk):
	densityFFT = pyfftw.builders.rfftn(densityField, threads = helpers.GetNumberOfThreads())
	densityFourier = densityFFT()
	centeredDensityFourier = np.fft.fftshift(densityFourier, axes=(0,1,))

	kShape = centeredDensityFourier.shape
	kRadii = [kShape[0]/2, kShape[1]/2, kShape[2]]

	fourierRadius = int((1.0*nGrid/2.0))
	fourierBins   = int(fourierRadius/dk)
	binnedPowerSpectrum = np.zeros(fourierBins)
	numInBin = [0] * fourierBins
	bin = 0

	for i in range(0, kShape[0]):
		kiSquare = (i-kRadii[0])**2
		for j in range(0, kShape[1]):
			kjSquare = (j-kRadii[1])**2
			for k in range(0, kShape[2]):
				kkSquare = k**2
				kr = math.sqrt(kiSquare + kjSquare + kkSquare)
				if kr != 0 and kr < fourierRadius:
					bin = int(kr / dk)
					numInBin[bin] += 1
					binnedPowerSpectrum[bin] += 2*((centeredDensityFourier[i][j][k].real)**2 + (centeredDensityFourier[i][j][k].imag)**2)

	for i in range(1, fourierBins):
		if numInBin[i] != 0:
			binnedPowerSpectrum[i] /= numInBin[i]

	#print binnedPowerSpectrum

	xScale = [((float(i) - 0.5)*dk/nGrid) for i in range(1, len(binnedPowerSpectrum)+1)]

	#plt.plot(np.array(xScale), np.array(binnedPowerSpectrum))
	#plt.show()

	return list(binnedPowerSpectrum)

def OutputPowerSpectrum(initialPowerSpectrum, middlePowerSpectrum, finalPowerSpectrum, startingA, nGrid, lBox, dk):
	xScale = [((float(i) - 0.5)*nGrid*dk/lBox) for i in range(1, len(finalPowerSpectrum)+1)]

	for index, i in enumerate(finalPowerSpectrum):
		print i

	plt.plot(np.array(xScale), np.array(finalPowerSpectrum))
	plt.plot(np.array(xScale), np.array(initialPowerSpectrum))
	plt.plot(np.array(xScale), np.array(middlePowerSpectrum))
	plt.xlabel('Fourier mode')
	plt.ylabel('Amplitude')
	plt.title('Dimensionless Density Contrast Power Spectra')
	plt.legend([r'$a=1.00$', r'$a=%.2f$' % ((startingA + 1) / 2.), r'$a=%.2f$' % startingA])
	plt.show()
