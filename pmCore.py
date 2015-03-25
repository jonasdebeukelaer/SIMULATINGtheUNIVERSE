import pmHelpers as helpers

import numpy as np
import math
import sys
import pyfftw

from operator import itemgetter

import plotly as py
from plotly.graph_objs import *
py.tools.set_credentials_file(username='jonasdb', api_key='2r0kiy0eq1', stream_ids=['fs43hhoxgj', '86nsdjyshj'])

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

def CalculatePowerSpectrum(densityField, nGrid, lBox):
	densityFFT = pyfftw.builders.rfftn(densityField, threads = helpers.GetNumberOfThreads())
	densityFourier = densityFFT()
	centeredDensityFourier = np.fft.fftshift(densityFourier)

	kShape = centeredDensityFourier.shape
	kRadii = [kShape[0]/2, kShape[1]/2, kShape[2]]
	numGridBoxes = kShape[0] * kShape[1] * kShape[2]
	centeredDensityFourier = centeredDensityFourier / nGrid**3

	distancesfromOrigin = list()

	for i in range(0, kShape[0]):
		kiSquare = (i-kRadii[0])**2
		for j in range(0, kShape[1]):
			kjSquare = (j-kRadii[1])**2
			for k in range(0, kShape[2]):
				kkSquare = k**2
				kr = math.sqrt(kiSquare + kjSquare + kkSquare)
				physicalK = kr / nGrid * lBox
				if kr != 0:
					distancesfromOrigin.append([physicalK, math.sqrt((centeredDensityFourier[i][j][k].real)**2 + (centeredDensityFourier[i][j][k].imag)**2)])

	distancesfromOrigin = sorted(distancesfromOrigin, key=itemgetter(0))

	resolution = 5
	numberFourierModes = int((lBox/resolution)+1)
	binnedPowerSpectrum = np.zeros(numberFourierModes)
	bin = 0
	topK = 0
	NumInBin = 0
	for i in distancesfromOrigin:
		physicalK = i[0]
		if physicalK > topK:
			if binnedPowerSpectrum[bin] != 0:
				binnedPowerSpectrum[bin] = float(binnedPowerSpectrum[bin]) / NumInBin
			topK += resolution
			bin += 1
			NumInBin = 0
		NumInBin += 1
		binnedPowerSpectrum[bin] += i[1]

	normalisedBinnedPowerSpectrum = binnedPowerSpectrum / max(binnedPowerSpectrum)

	return list(binnedPowerSpectrum)

def OutputPowerSpectrumHeatMap(powerSpectrumHeatMap, aArray, nGrid, lBox):
	numberFourierModes = int(math.sqrt(3)*(nGrid/2+1)+1)

	heatmap  = Heatmap(z=powerSpectrumHeatMap, y=aArray, x=[float(i)*lBox/nGrid for i in range(1, numberFourierModes)], colorscale='Greens')
	layout   = Layout(title='Power Spectrum Evolution In Time (nGrid=%s, lBox=%s)' % (nGrid, lBox), xaxis=XAxis(title='$Fourier Mode /Mpcs^{-1}$'), yaxis=YAxis(title='Scale Factor (a)'))
	data     = Data([heatmap])
	fig      = Figure(data=data, layout=layout)
	plot_url = py.plotly.plot(fig, filename='Power Spectrum Evolution (nGrid=%s, lBox=%s)' % (nGrid, lBox), world_readable=True)

def OutputPowerSpectrum(initialPowerSpectrum, finalPowerSpectrum, startingA, nGrid, lBox):
	resolution = 5
	chart1   = Scatter(x=[float(i)*resolution for i in range(1, len(finalPowerSpectrum))], y=initialPowerSpectrum, name='Initial (a=%d)' % (startingA), mode='line+marker', line=Line(shape='spline'))
	chart2   = Scatter(x=[float(i)*resolution for i in range(1, len(finalPowerSpectrum))], y=finalPowerSpectrum, name='Final (a=1)',  mode='line+marker', line=Line(shape='spline'))
	layout   = Layout(title='Power Spectrum of a Universe Simulation (nGrid=%s, lBox=%s)' % (nGrid, lBox), xaxis=XAxis(title='$Fourier Mode /Mpcs^{-1}$', rangemode='tozero'), barmode='group')
	data     = Data([chart1, chart2])
	fig      = Figure(data=data, layout=layout)
	plot_url = py.plotly.plot(fig, filename='InitialFinal Power Spectrums (nGrid=%s, lBox=%s)' % (nGrid, lBox), world_readable=True)
