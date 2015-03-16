import pmHelpers as helpers

import numpy as np
import math
import sys
import pyfftw

def FindMeshIndex(position, nGrid):
	index = math.floor(position + 0.5) + ((nGrid / 2) - 1)
	if index == -1:
		index = nGrid - 1
	return index

def CalculateDensityField(nGrid, particleList):
	densityFieldMesh = np.zeros([nGrid, nGrid, nGrid])

	for particle in particleList:
		xMesh = FindMeshIndex(particle.position[0], nGrid)
		yMesh = FindMeshIndex(particle.position[1], nGrid)
		zMesh = FindMeshIndex(particle.position[2], nGrid)

		densityFieldMesh[xMesh][yMesh][zMesh] += particle.mass

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

	potentialField = np.fft.irfftn(densityFieldConvoluted)
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

	xAcceleration = - (potentialField[xNeighbours[0]][yMesh][zMesh] - potentialField[xNeighbours[1]][yMesh][zMesh]) / (2.0 * particle.mass)
	yAcceleration = - (potentialField[xMesh][yNeighbours[0]][zMesh] - potentialField[xMesh][yNeighbours[1]][zMesh]) / (2.0 * particle.mass)
	zAcceleration = - (potentialField[xMesh][yMesh][zNeighbours[0]] - potentialField[xMesh][yMesh][zNeighbours[1]]) / (2.0 * particle.mass)

	return (xAcceleration, yAcceleration, zAcceleration)

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