import numpy as np
import math
import sys
import random
import time

G = 1

class Particle:
	def __init__(self, position, velocity, mass):
		self.position = position
		self.velocity = velocity
		self.mass = mass
		self.accumulatedForce = float(0)
		self.acceleration = [0.0, 0.0, 0.0]


def CalcForce(pos1, pos2, mass1, mass2):

	epsilon = 2
	separation = np.subtract(pos2, pos1)
	separationMagnitudeSquared = sum(separation**2) + epsilon**2
	force = G * mass1 * mass2 / separationMagnitudeSquared**1.5
	force = np.multiply(separation, force)
	return force

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

start = time.time()
random.seed(89321)

numParticles = 1
initialisationResolution = 0.1
maxCoordinate = 50
randomMaxCoordinates = maxCoordinate / initialisationResolution
maxVelocity = 1
positionDistribution = 2
velocityDistribution = 1
hasCenterParticle = True

particleList = InitialiseParticles(numParticles, initialisationResolution, randomMaxCoordinates, maxVelocity, positionDistribution, velocityDistribution)
if hasCenterParticle:
	centreParticle = Particle([0., 0., 0.], [0., 0., 0.,], 20)
	particleList.append(centreParticle)
	numParticles += 1

timeStepSize = 0.01
numTimeSteps = 100000
shootEvery = 100
#time = 0.0

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

	for i, particle in enumerate(particleList):

		for otherIndex in range(i + 1, numParticles):

			otherParticle = particleList[otherIndex]
			force = CalcForce(particle.position, otherParticle.position, particle.mass, otherParticle.mass)
			particle.accumulatedForce += force
			otherParticle.accumulatedForce -= force

		particleAcceleration = particle.accumulatedForce / particle.mass
		
		if timeStep == 0:
			particle.acceleration = particleAcceleration

		particle.position += np.multiply(timeStepSize, particle.velocity) + np.multiply(0.5 * timeStepSize**2, particle.acceleration)
		particle.velocity += np.multiply((0.5 * timeStepSize), (np.add(particle.acceleration, particleAcceleration)))
		velocityMagnitude = ((particle.velocity[0])**2 + (particle.velocity[1])**2 + (particle.velocity[2])**2)
		particle.acceleration = particleAcceleration
		particle.accumulatedForce = [0.0, 0.0, 0.0]

		if shoot:
			f.write("%f %f %f %f\n" % (particle.position[0], particle.position[1], particle.position[2], velocityMagnitude))

	if shoot:
		f.write("%f %f %f %f\n%f %f %f %f\n" % (maxCoordinate, maxCoordinate, maxCoordinate, 0., -maxCoordinate, -maxCoordinate, -maxCoordinate, 0.))
		f.close()

end = time.time()
sys.stdout.write("\n")
print end - start





