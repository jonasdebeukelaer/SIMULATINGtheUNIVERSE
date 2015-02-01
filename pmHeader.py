import numpy as np
import math
import os
import sys
import random
import time
import pyfftw
from pync import Notifier
from enum import Enum

G = 1

class Particle:
	def __init__(self, position, velocity, mass):
		self.position     = position
		self.velocity     = velocity
		self.mass         = mass
		self.acceleration = [0.0, 0.0, 0.0]

class PositionDist(Enum):
	random       = 1
	randomShell  = 2
	evenShell    = 3
	correctField = 4

class VelocityDist(Enum):
	random       = 1
	zero         = 2
	correctField = 3

def InitialiseParticles(volume, initialisationResolution, numParticles, positionDistribution, velocityDistribution, maxVelocity):
	particleList = []
	
	if positionDistribution == PositionDist.correctField:
		check = numParticles**(1/3)
		if check != int(check):
			print 'Not a cubed number of particles...\n'
			check += 0.5
			check = round(check)
			numParticles = check**3
			print 'New numParticles = %d' % (numParticles)

	for i in range(0, numParticles):
		newParticle = Particle([0., 0., 0.] ,[0., 0., 0.], 1)
		particleList.append(newParticle)

	# Random coordinates distribution
	if positionDistribution == PositionDist.random:
		for particle in particleList:
			x                 = random.randrange((- volume[0] / 2) / initialisationResolution, (volume[0] / 2) / initialisationResolution) * initialisationResolution
			y                 = random.randrange((- volume[1] / 2) / initialisationResolution, (volume[1] / 2) / initialisationResolution) * initialisationResolution
			z                 = random.randrange((- volume[2] / 2) / initialisationResolution, (volume[2] / 2) / initialisationResolution) * initialisationResolution
			particle.position = [x, y, z]

	# Random shell of particles distribution
	elif positionDistribution == PositionDist.randomShell:
		r = min([volume[0] / 2, volume[1] / 2, volume[2] / 2])
		for particle in particleList:
			theta = math.acos(random.uniform(-1., 1.))
			phi   = random.uniform(0., 2 * math.pi)
			
			x     = r * math.sin(theta) * math.cos(phi)
			y     = r * math.sin(theta) * math.sin(phi)
			z     = r * math.cos(theta)

			particle.position = [x, y, z]

	# Even shell of particles distribution
	elif positionDistribution == PositionDist.evenShell:
		r            = min([volume[0] / 4, volume[1] / 4, volume[2] / 4])
		
		golden_angle = np.pi * (3 - np.sqrt(5))
		theta        = golden_angle * np.arange(numParticles)
		z            = np.linspace(r - 1.0 / numParticles, 1.0 / numParticles - r, numParticles)
		radius       = np.sqrt((r**2 - z * z))
		 
		points       = np.zeros((numParticles, 3))
		points[:,0]  = radius * np.cos(theta)
		points[:,1]  = radius * np.sin(theta)
		points[:,2]  = z

		i = 0
		for particle in particleList:
			particle.position = points[i, :]
			i += 1

	elif positionDistribution == 3:
		numPerSide = numParticles**(1/3)
		particleNumber = 0
		for i in range(0, numPerSide):
			x = volume[0]/2 - i*volume[0]/numPerSide
			for j in range(0, numPerSide):
				y = volume[1]/2 - i*volume[1]/numPerSide
				for k in range(0, numPerSide):
					z = volume[2]/2 - i*volume[2]/numPerSide
					particleList[particleNumber].position = [x,y,z]
					particleNumber+=1




	else:
		print 'Invalid position distribution selected'

	# Random velocity distribution
	if velocityDistribution == VelocityDist.random:
		for particle in particleList:
			xVelocity         = random.randrange(- maxVelocity / initialisationResolution, maxVelocity / initialisationResolution) * initialisationResolution
			yVelocity         = random.randrange(- maxVelocity / initialisationResolution, maxVelocity / initialisationResolution) * initialisationResolution
			zVelocity         = random.randrange(- maxVelocity / initialisationResolution, maxVelocity / initialisationResolution) * initialisationResolution
			particle.velocity = [xVelocity, yVelocity, zVelocity]

	# Zero velocity distribution
	elif velocityDistribution == VelocityDist.zero:
		for particle in particleList:
			particle.velocity = [0, 0, 0]

	else:
		print 'Invalid velocity distribution selected'

	return particleList

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
	shape       = (unalteredShape[0], unalteredShape[1], unalteredShape[2])# / 2) + 1)
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

def GetNumberOfThreads():
	user = os.getlogin()
	if user == "oliclipsham":
		threads = 8
	elif user == "jonasdebeukelaer":
		threads = 4
	else:
		print "Imposter!"
		threads = 1

	return threads

def SolvePotential(densityField, greensFunction):
	densityFieldFFT        = pyfftw.builders.rfftn(densityField, threads = GetNumberOfThreads())
	densityFieldConvoluted = np.multiply(greensFunction, densityFieldFFT())
	potentialFieldJumbled  = pyfftw.builders.irfftn(densityFieldConvoluted, threads = GetNumberOfThreads())
	potentialField         = np.fft.fftshift(potentialFieldJumbled())
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

def PositionCorrect(particle, volumeLimits):
	positionLimits = np.array(volumeLimits) / 2
	for index, position in enumerate(particle.position):
		if position > positionLimits[index]:
			particle.position[index] = position - volumeLimits[index]
		elif position < (- positionLimits[index]):
			particle.position[index] = position + volumeLimits[index]

def OutputPercentage(timeStep, numTimeSteps):
	i = (float(timeStep + 1) / numTimeSteps) * 100
	sys.stdout.write("\r%.2f%%" % i)
	sys.stdout.flush()

def OutputPotentialFieldXY(potentialField, particleList, volume, timeStep, gridResolution):
	g = open("PotentialResults/potential_frame%d.3D" % (timeStep), "w")
	g.write("x y z Potential\n")
	
	for i in range(0, volume[0]):
		for j in range(0, volume[1]):
			potentialValue = potentialField[i][j][49]
			if potentialValue > 0.01 and i % 2 == 0 and j % 2 == 0:
				g.write("%f %f 0 %f\n" % (i*gridResolution-((volume[0]/2)-1), j*gridResolution-((volume[1]/2)-1), potentialValue))

	for particle in particleList:
		g.write("%f %f %f %f\n" % (particle.position[0], particle.position[1], particle.position[2], 0.1))

	#g.write("%f %f %f %f\n%f %f %f %f\n" % (volume[0] / 2, volume[1] / 2, volume[2] / 2, 0., - volume[0] / 2, - volume[1] / 2, - volume[2] / 2, 0.))
	g.close()