import pmHelpers as helpers

import numpy as np
import math
import sys
import pyfftw

def FindMeshIndex(position, gridResolution, gridSize):
	index = round((position / gridResolution) + (gridResolution / 2)) + ((gridSize / 2) - 1)
	if index == -1:
		index = gridSize - 1
	return index

def CalculateDensityField(volume, gridResolution, particleList, populateArray = True):
	meshShape = [volume[0] / gridResolution, volume[1] / gridResolution, volume[2] / gridResolution]
	if meshShape[0] != int(meshShape[0]) or meshShape[1] != int(meshShape[1]) or meshShape[2] != int(meshShape[2]):
		sys.exit("Error: non-integer cell number defined pleuz fix")

	densityFieldMesh = np.zeros((meshShape))

	if populateArray:
		for particle in particleList:
			xMesh = FindMeshIndex(particle.position[0], gridResolution, meshShape[0])
			yMesh = FindMeshIndex(particle.position[1], gridResolution, meshShape[1])
			zMesh = FindMeshIndex(particle.position[2], gridResolution, meshShape[2])

			densityFieldMesh[xMesh][yMesh][zMesh] += particle.mass

		densityFieldMesh /= (gridResolution**3)

	return densityFieldMesh

def CreateGreensFunction(unalteredShape):
	shape       = (unalteredShape[0], unalteredShape[1], unalteredShape[2] / 2 + 1)	
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

	potentialField = np.fft.irfftn(densityFieldConvoluted)
	return potentialField

def FindPlusMinus(meshIndex, axisSize):
	if meshIndex == 0:
		meshPlus  = meshIndex + 1
		meshMinus = axisSize - 1
	elif meshIndex == axisSize - 1:
		meshPlus  = 0
		meshMinus = meshIndex - 1
	else:
		meshPlus  = meshIndex + 1
		meshMinus = meshIndex - 1

	return (meshPlus, meshMinus)

def CalculateParticleAcceleration(particle, potentialField, gridResolution):
	meshShape = potentialField.shape

	xMesh = FindMeshIndex(particle.position[0], gridResolution, meshShape[0])
	yMesh = FindMeshIndex(particle.position[1], gridResolution, meshShape[1])
	zMesh = FindMeshIndex(particle.position[2], gridResolution, meshShape[2])

	xNeighbours = FindPlusMinus(xMesh, meshShape[0])
	yNeighbours = FindPlusMinus(yMesh, meshShape[1])
	zNeighbours = FindPlusMinus(zMesh, meshShape[2])

	xAcceleration = - (potentialField[xNeighbours[0]][yMesh][zMesh] - potentialField[xNeighbours[1]][yMesh][zMesh]) / (2.0 * particle.mass)
	yAcceleration = - (potentialField[xMesh][yNeighbours[0]][zMesh] - potentialField[xMesh][yNeighbours[1]][zMesh]) / (2.0 * particle.mass)
	zAcceleration = - (potentialField[xMesh][yMesh][zNeighbours[0]] - potentialField[xMesh][yMesh][zNeighbours[1]]) / (2.0 * particle.mass)

	return (xAcceleration, yAcceleration, zAcceleration)

def GetF(a, omega_m = 1, omega_k = 0, omega_lambda = 0):
	return (a**(-1)*(omega_m + omega_k * a + omega_lambda * a**3))**(-0.5)

def PositionCorrect(particle, volumeLimits):
	positionLimits = np.array(volumeLimits) / 2
	for index, position in enumerate(particle.position):
		if position > positionLimits[index]:
			particle.position[index] = position - volumeLimits[index]
		elif position < (- positionLimits[index]):
			particle.position[index] = position + volumeLimits[index]