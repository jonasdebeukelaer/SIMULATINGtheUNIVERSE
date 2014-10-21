import numpy as np
import cmath
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

numParticles = 2
particleList = []
timeStepSize = 0.1

particle1 = Particle([2.0, 0.0, 0.0], [0.0, 10.0, 0.0], 1)
particle2 = Particle([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 100)

particleList.append(particle1)
particleList.append(particle2)

file_ = open('results.txt', 'w')

time = 0.0

for timeStep in range(0, 10):
	time += timeStepSize
	for i, particle in enumerate(particleList):

		for otherIndex in range(i+1, numParticles):

			otherParticle = particleList[otherIndex]
			force = CalcForce(particle.position, otherParticle.position, particle.mass, otherParticle.mass)
			particle.accumulatedForce += force
			otherParticle.accumulatedForce -= force
			
			#print str(time) + "\t" + str(force[0]) + "\t" + str(particle.accumulatedForce[0]) + "\t" + str(otherParticle.accumulatedForce[0])

		particleAcceleration = particle.accumulatedForce / particle.mass
		#if particle.mass == 1:
			#print particleAcceleration
		
		if timeStep == 0:
			particle.acceleration = particleAcceleration
		particle.position += np.multiply(particle.velocity, timeStepSize) + np.multiply(0.5 * timeStepSize**2, particle.acceleration)
		particle.velocity += np.multiply((0.5 * timeStepSize), (np.add(particle.acceleration, particleAcceleration)))
		particle.acceleration = particleAcceleration
		particle.accumulatedForce = [0.0, 0.0, 0.0]

		if particle.mass == 1:
			print str(time) + "\t" + str(particle.position[0]) + "\t" + str(particle.velocity[0]) + "\t" + str(particle.acceleration[0])
			file_.write(str(time) + "\t" + str(particle.position[0]) + "\t" + str(particle.position[1]) + "\t" + str(particle.velocity[0]) + "\t" + str(particle.acceleration[0]) + '\n')
file_.close()














