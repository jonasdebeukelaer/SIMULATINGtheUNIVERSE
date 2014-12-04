import numpy as np
import math
import sys
import random
import time
import pyfftw

G = 1

class Particle:
	def __init__(self, position, velocity, mass):
		self.position = position
		self.velocity = velocity
		self.mass = mass
		self.acceleration = [0.0, 0.0, 0.0]

def InitialiseParticles(numParticles, initialisationResolution, maxCoordinate, maxVelocity, positionDistribution, velocityDistribution):
	particleList = []
	#maxCoordinates is randomMaxCoordinates in the main section
	
	for i in range(0, numParticles):
		newParticle = Particle([0., 0., 0.] ,[0., 0., 0.], 1)
		particleList.append(newParticle)

	#random coordinates dist
	if positionDistribution == 0:
		for particle in particleList:
			x = random.randrange(-maxCoordinate, maxCoordinate) * initialisationResolution
			y = random.randrange(-maxCoordinate, maxCoordinate) * initialisationResolution
			z = random.randrange(-maxCoordinate, maxCoordinate) * initialisationResolution
			particle.position = [x, y, z]

	#random shell of particles distribution
	elif positionDistribution == 1:
		r = maxCoordinate * initialisationResolution
		for particle in particleList:
			theta = math.acos(random.uniform(-1., 1.))
			phi = random.uniform(0., 2 * math.pi)
			
			x = r * math.sin(theta) * math.cos(phi)
			y = r * math.sin(theta) * math.sin(phi)
			z = r * math.cos(theta)

			particle.position = [x, y, z]

	#Even shell of particles distribution
	elif positionDistribution == 2:
		r = maxCoordinate * initialisationResolution
		
		golden_angle = np.pi * (3 - np.sqrt(5))
		theta = golden_angle * np.arange(numParticles)
		z = np.linspace(r - 1.0 / numParticles, 1.0 / numParticles - r, numParticles)
		radius = np.sqrt((r**2 - z * z))
		 
		points = np.zeros((numParticles, 3))
		points[:,0] = radius * np.cos(theta)
		points[:,1] = radius * np.sin(theta)
		points[:,2] = z

		i = 0
		for particle in particleList:
			particle.position = points[i, :]
			i+=1

	else:
		print 'Invalid position distribution selected'

	#random velocity dist
	if velocityDistribution == 0:
		for particle in particleList:
			xVelocity = random.randrange(-maxVelocity, maxVelocity) * initialisationResolution
			yVelocity = random.randrange(-maxVelocity, maxVelocity) * initialisationResolution
			zVelocity = random.randrange(-maxVelocity, maxVelocity) * initialisationResolution
			particle.velocity = [xVelocity, yVelocity, zVelocity]

	#zero velocity distribution
	elif velocityDistribution == 1:
		for particle in particleList:
			particle.velocity = [0, 0, 0]

	else:
		print 'Invalid velocity distribution selected'

	return particleList

def FindMeshIndex(position, gridResolution, gridSize): 
	return round((position / gridResolution) + (gridResolution / 2)) + ((gridSize / 2) - 1)

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
	shape = (unalteredShape[0], unalteredShape[1], (unalteredShape[2] / 2) + 1)
	greensArray = np.zeros((shape))

	constant = 1

	for l in range(0, shape[0]):
		if l <= shape[0] / 2:
			kx = 2 * math.pi * l / (shape[0] - 1)
		else:
			kx = 2 * math.pi * (l - shape[0]) / (shape[0] - 1)

		for m in range(0, shape[1]):
			if m <= shape[1] / 2:
				ky = 2 * math.pi * m / (shape[1] - 1)
			else:
				ky = 2 * math.pi * (m - shape[1]) / (shape[1] - 1)

			for n in range(0, shape[2]):
				if n <= shape[2] / 2:
					kz = 2 * math.pi * n / (shape[2] - 1)
				else:
					kz = 2 * math.pi * (n - shape[2]) / (shape[2] - 1)

				if l != 0 or m != 0 or n != 0:
					greensArray[l][m][n] = - constant / ((math.sin(kx * 0.5))**2 + (math.sin(ky * 0.5))**2 + (math.sin(kz * 0.5))**2)

	return greensArray

def SolvePotential(densityField, greensFunction):
	densityFieldFFT = pyfftw.builders.rfftn(densityField)
	densityFieldConvoluted = densityFieldFFT() * greensFunction
	potentialFieldJumbled = pyfftw.builders.irfftn(densityFieldConvoluted)
	potentialField = np.fft.fftshift(potentialFieldJumbled())
	return potentialField

def ConvertPotentialToForce(potentialField, gridResolution):
	volume = potentialField.shape
	forceField = np.zeros((volume[0], volume[1], volume[2], 3))

	for i in range(0, volume[0]):
		for j in range(0, volume[1]):
			for k in range(0, volume[2]):
				if i != 0 and i != (volume[0] - 1) and j != 0 and j != (volume[1] - 1) and k != 0 and k != (volume[2] - 1):
					xForce = - (potentialField[i+1][j][k] - potentialField[i-1][j][k]) / (2 * gridResolution)
					yForce = - (potentialField[i][j+1][k] - potentialField[i][j-1][k]) / (2 * gridResolution)
					zForce = - (potentialField[i][j][k+1] - potentialField[i][j][k-1]) / (2 * gridResolution)

					forceField[i][j][k][0] = xForce
					forceField[i][j][k][1] = yForce
					forceField[i][j][k][2] = zForce

	return forceField

def CalculateParticleAcceleration(particle, forceField):
	xMesh = int(particle.position[0] / gridResolution)
	yMesh = int(particle.position[1] / gridResolution)
	zMesh = int(particle.position[2] / gridResolution)

	particleAcceleration = (forceField[xMesh][yMesh][zMesh][0]/particle.mass, forceField[xMesh][yMesh][zMesh][1]/particle.mass, forceField[xMesh][yMesh][zMesh][2]/particle.mass)
	return particleAcceleration

start = time.time()
print "Seeding..."
random.seed(89321)
print "Done\n"

numParticles = 100
initialisationResolution = 0.1
maxCoordinate = 25
randomMaxCoordinates = maxCoordinate / initialisationResolution
maxVelocity = 1
positionDistribution = 0
velocityDistribution = 1
hasCenterParticle = False

print "Initialising particles..."
particleList = InitialiseParticles(numParticles, initialisationResolution, randomMaxCoordinates, maxVelocity, positionDistribution, velocityDistribution)
if hasCenterParticle:
	centreParticle = Particle([100., 100., 100.], [0., 0., 0.,], 20)
	particleList.append(centreParticle)
	numParticles += 1

particleFile = open("PotentialInvestigation0.3D", 'w')
particleFile.write("x\ty\tz\tValue\n")
particleFile.write("50\t50\t50\t0\n-50\t-50\t-50\t0\n")
for particle in particleList:
	particleFile.write("%f\t%f\t%f\t%f\n" % (particle.position[0], particle.position[1], particle.position[2], particle.mass))
particleFile.close()
print"Done\n"

timeStepSize = 0.01
numTimeSteps = 1
shootEvery = 100

volume = [100, 100, 100]
gridResolution = 1

print "Determining mesh shape..."
densityField = CalculateDensityField(volume, gridResolution, particleList, False)
print "Done\n"

print "Calculating Green's function..."
greensFunction = CreateGreensFunction(densityField.shape)
print "Done\n"

print "Iterating..."
for timeStep in range(0, numTimeSteps):
	#time += timeStepSize
	#shoot = True if (timeStep % shootEvery) == 0 else False

	#percentage counter
	i = (float(timeStep) / numTimeSteps) * 100
	sys.stdout.write("\r%.2f%%" % i)
	sys.stdout.flush()

	#if shoot:
		#f = open("Results/values_frame%d.3D" % (timeStep), "w")
		#f.write("x y z VelocityMagnitude\n")

	densityField   = CalculateDensityField(volume, gridResolution, particleList)
	densityFile = open("PotentialInvestigation1.3D", 'w')
	densityFile.write("x\ty\tz\tValue\n")
	densityFile.write("50\t50\t50\t0\n-50\t-50\t-50\t0\n")
	densityShape = densityField.shape
	for i in range(0, densityShape[0]):
		x = (i - ((densityShape[0] / 2) - 1)) * gridResolution
		for j in range(0, densityShape[1]):
			y = (j - ((densityShape[1] / 2) - 1)) * gridResolution
			for k in range(0, densityShape[2]):
				z = (k - ((densityShape[2] / 2) - 1)) * gridResolution
				if densityField[i][j][k] != 0:
					densityFile.write("%f\t%f\t%f\t%f\n" % (x, y, z, densityField[i][j][k]))
	densityFile.close()

	potentialField = SolvePotential(densityField, greensFunction)
	potentialFile = open("PotentialInvestigation2.3D", 'w')
	potentialFile.write("x\ty\tz\tValue\n")
	potentialFile.write("50\t50\t50\t0\n-50\t-50\t-50\t0\n")
	potentialShape = potentialField.shape
	if potentialShape != densityShape:
		sys.exit("Error: something definitely went wrong")
	for i in range(0, potentialShape[0]):
		x = (i - ((potentialShape[0] / 2) - 1)) * gridResolution
		for j in range(0, potentialShape[1]):
			y = (j - ((potentialShape[1] / 2) - 1)) * gridResolution
			for k in range(0, potentialShape[2]):
				z = (k - ((potentialShape[2] / 2) - 1)) * gridResolution
				if potentialField[i][j][k] != 0:
					potentialFile.write("%f\t%f\t%f\t%f\n" % (x, y, z, potentialField[i][j][k]))
	potentialFile.close()

	#forceField     = ConvertPotentialToForce(potentialField, gridResolution)

	#for particle in particleList:

		#particleAcceleration = CalculateParticleAcceleration(particle, forceField)
		
		#if timeStep == 0:
			#particle.acceleration = particleAcceleration

		#particle.position += np.multiply(timeStepSize, particle.velocity) + np.multiply(0.5 * timeStepSize**2, particle.acceleration)
		#particle.velocity += np.multiply((0.5 * timeStepSize), (np.add(particle.acceleration, particleAcceleration)))
		#velocityMagnitude = ((particle.velocity[0])**2 + (particle.velocity[1])**2 + (particle.velocity[2])**2)
		#particle.acceleration = particleAcceleration

		#if shoot:
			#f.write("%f %f %f %f\n" % (particle.position[0], particle.position[1], particle.position[2], velocityMagnitude))

	#if shoot:
		#f.write("%f %f %f %f\n%f %f %f %f\n" % (maxCoordinate, maxCoordinate, maxCoordinate, 0., -maxCoordinate, -maxCoordinate, -maxCoordinate, 0.))
		#f.close()

end = time.time()
sys.stdout.write("\n")
print end - start









