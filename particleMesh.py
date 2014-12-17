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

def SolvePotential(densityField, greensFunction, shoot, timeStep, outputFourierSpace=False):
	densityFieldFFT = pyfftw.builders.rfftn(densityField)
	densityFieldConvoluted = densityFieldFFT() * greensFunction

	if outputFourierSpace and shoot:
		#shiftedConvolution = np.fft.fftshift(densityFieldConvoluted)
		OutputFourierSpaceTo3D(densityFieldConvoluted, timeStep)
	
	potentialFieldJumbled = pyfftw.builders.irfftn(densityFieldConvoluted)
	
	potentialField = np.fft.fftshift(potentialFieldJumbled())
	return potentialFieldJumbled()

def OutputFourierSpaceTo3D(densityFieldConvoluted, timeStep):
	fourierFile = open('fourierResults/convolutedFourierSpace_%s.3D' % timeStep, 'w')
	fourierFile.write("x y z VelocityMagnitude\n")
	for kx in range(0, len(densityFieldConvoluted[:][1][1])):
		for ky in range(0, len(densityFieldConvoluted[1][:][1])):
			for kz in range(0, len(densityFieldConvoluted[1][1][:])):
				if densityFieldConvoluted[kx][ky][kz] > 0.9:
					fourierFile.write("%d %d %d %f\n" % (kx, ky, kz, math.log(densityFieldConvoluted[kx][ky][kz])))

	fourierFile.close()

def CalculateParticleAcceleration(particle, potentialField, gridResolution):
	meshShape = potentialField.shape

	xMesh = FindMeshIndex(particle.position[0], gridResolution, meshShape[0])
	yMesh = FindMeshIndex(particle.position[1], gridResolution, meshShape[1])
	zMesh = FindMeshIndex(particle.position[2], gridResolution, meshShape[2])

	xAcceleration = - (potentialField[xMesh + 1][yMesh][zMesh] - potentialField[xMesh - 1][yMesh][zMesh]) / 2.0
	yAcceleration = - (potentialField[xMesh][yMesh + 1][zMesh] - potentialField[xMesh][yMesh - 1][zMesh]) / 2.0
	zAcceleration = - (potentialField[xMesh][yMesh][zMesh + 1] - potentialField[xMesh][yMesh][zMesh - 1]) / 2.0

	return (xAcceleration, yAcceleration, zAcceleration)

start = time.time()
print "Seeding..."
random.seed(89321)
print "Done\n"

numParticles = 20
initialisationResolution = 0.1
maxCoordinate = 25
randomMaxCoordinates = maxCoordinate / initialisationResolution
maxVelocity = 1
positionDistribution = 0
velocityDistribution = 1
hasCenterParticle = False
printFourierSpace = True

print "Initialising particles..."
particleList = InitialiseParticles(numParticles, initialisationResolution, randomMaxCoordinates, maxVelocity, positionDistribution, velocityDistribution)
if hasCenterParticle:
	centreParticle = Particle([0., 0., 0.], [0., 0., 0.,], 20)
	particleList.append(centreParticle)
	numParticles += 1
print"Done\n"

timeStepSize = 0.01
numTimeSteps = 1
shootEvery = 200

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
	shoot = True if (timeStep % shootEvery) == 0 else False

	#percentage counter
	i = (float(timeStep) / numTimeSteps) * 100
	sys.stdout.write("\r%.2f%%" % i)
	sys.stdout.flush()

	if shoot:
		f = open("Results/values_frame%d.3D" % (timeStep), "w")
		f.write("x y z VelocityMagnitude\n")

	densityField   = CalculateDensityField(volume, gridResolution, particleList)
	potentialField = SolvePotential(densityField, greensFunction, shoot, timeStep, printFourierSpace)

	for particle in particleList:

		particleAcceleration = CalculateParticleAcceleration(particle, potentialField, gridResolution)
		
		if timeStep == 0:
			particle.acceleration = particleAcceleration

		particle.position += np.multiply(timeStepSize, particle.velocity) + np.multiply(0.5 * timeStepSize**2, particle.acceleration)
		particle.velocity += np.multiply((0.5 * timeStepSize), (np.add(particle.acceleration, particleAcceleration)))
		velocityMagnitude = ((particle.velocity[0])**2 + (particle.velocity[1])**2 + (particle.velocity[2])**2)
		particle.acceleration = particleAcceleration

		if shoot:
			f.write("%f %f %f %f\n" % (particle.position[0], particle.position[1], particle.position[2], velocityMagnitude))

	if shoot:
		f.write("%f %f %f %f\n%f %f %f %f\n" % (volume[0] / 2, volume[1] / 2, volume[2] / 2, 0., - volume[0] / 2, - volume[1] / 2, - volume[2] / 2, 0.))
		f.close()

end = time.time()
sys.stdout.write("\n")
print end - start









