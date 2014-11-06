import numpy as np
import math
import sys

G = 1

class Particle:
	def __init__(self, position, velocity, mass):
		self.position = position
		self.velocity = velocity
		self.mass = mass
		self.accumulatedForce = float(0)
		self.acceleration = [0.0, 0.0, 0.0]


def CalcForce(pos1, pos2, mass1, mass2):

	separation = np.subtract(pos2, pos1)
	separationMagnitudeSquared = sum(separation**2)
	force = G * mass1 * mass2 / separationMagnitudeSquared
	unitVector = np.divide(separation, separationMagnitudeSquared**0.5)
	force = np.multiply(unitVector, force)
	return force

particleList = []

particle1 = Particle([50.0, 0.0, 0.0], [0.0, 1.4135, 0.0], 1)
particle2 = Particle([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 100)

particleList.append(particle1)
particleList.append(particle2)

numParticles = len(particleList)

timeStepSize = 0.01
numTimeSteps = 50000
time = 0.0

for timeStep in range(0, numTimeSteps):
	time += timeStepSize

	f = open("Results/values_frame%d.3D" % (timeStep), "wt")
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

		if particle.mass == 100:
			particle.position = [0.0, 0.0, 0.0]

		f.write("%f %f %f %f\n" % (particle.position[0], particle.position[1], particle.position[2], velocityMagnitude))

	f.close()










