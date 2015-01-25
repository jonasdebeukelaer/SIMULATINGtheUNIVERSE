import pmHeader as pm

import numpy as np
import math
import os
import sys
import random
import time
import pyfftw
from pync import Notifier

start = time.time()
print "Seeding..."
random.seed(89321)
print "Done\n"

#------------INITIALISATION PARAMETERS-----------#

volume                   = [2.5, 2.5, 2.5]
initialisationResolution = 1
gridResolution           = 0.01

numParticles             = 1
positionDistribution     = pm.PositionDist.random
velocityDistribution     = pm.VelocityDist.random

maxVelocity              = 1
hasCenterParticle        = True

numTimeSteps             = 100000
timeStepSize             = 0.001
shootEvery               = 1

#------------------------------------------------#

#------------INITIALISATION FUNCTIONS------------#

print "Initialising particles..."
#particleList = pm.InitialiseParticles(volume, initialisationResolution, numParticles, positionDistribution, velocityDistribution, maxVelocity)
particleList = []
particle1 = pm.Particle([0.25, 0., 0.], [0., 3., 0.], 1)
particleList.append(particle1)
if hasCenterParticle:
	centerParticle = pm.Particle([0., 0., 0.], [0., 0., 0.,], 20)
	particleList.append(centerParticle)
	numParticles += 1
print"Done\n"

print "Determining mesh shape..."
densityField = pm.CalculateDensityField(volume, gridResolution, particleList, False)
print densityField.shape
print "Done\n"

print "Calculating Green's function..."
greensFunction = pm.CreateGreensFunction(densityField.shape)
print "Done\n"

if os.path.exists("Results/values_frame0.3D"):
	raw_input("Delete yo motherflippin results from the last test, you simpleton! Or, if you'rereally sure, just hit enter I guess...\n")

#------------------------------------------------#

#-----------------ITERATION LOOP-----------------#

print "Iterating..."
timeStep = 0
while timeStep < numTimeSteps:

	shoot = True if (timeStep % shootEvery) == 0 else False
	if shoot:
		f = open("Results/values_frame%d.3D" % (timeStep), "w")
		f.write("x y z VelocityMagnitude\n")

	densityField   = pm.CalculateDensityField(volume, gridResolution, particleList)
	potentialField = pm.SolvePotential(densityField, greensFunction)

	for particle in particleList:

		particleAcceleration = pm.CalculateParticleAcceleration(particle, potentialField, gridResolution)
		if timeStep == 0:
			particle.acceleration = particleAcceleration

		particle.position     += np.multiply(timeStepSize, particle.velocity) + np.multiply(0.5 * timeStepSize**2, particle.acceleration)
		pm.PositionCorrect(particle, volume)
		particle.velocity     -= np.multiply((0.5 * timeStepSize), (np.add(particle.acceleration, particleAcceleration))) # THIS IS UBER WRONG
		velocityMagnitude     =  ((particle.velocity[0])**2 + (particle.velocity[1])**2 + (particle.velocity[2])**2)
		particle.acceleration =  particleAcceleration

		if particle.mass == 20:
			particle.position = [0, 0, 0]
			particle.velocity = [0, 0, 0]
		if shoot:
			f.write("%f %f %f %f\n" % (particle.position[0], particle.position[1], particle.position[2], velocityMagnitude))

	if shoot:
		f.write("%f %f %f %f\n%f %f %f %f\n" % (volume[0] / 2, volume[1] / 2, volume[2] / 2, 0., - volume[0] / 2, - volume[1] / 2, - volume[2] / 2, 0.))
		f.close()

	pm.OutputPercentage(timeStep, numTimeSteps)

	timeStep += 1

	if timeStep == numTimeSteps:
		Notifier.notify('%d time steps complete' % (numTimeSteps), title = 'User input required')
		moreSteps = raw_input("\n\nPlease check your VisIt output... would you like to add more time steps         (currently %d run)? (y/n): " % (numTimeSteps))
		if moreSteps == 'y':
			extraSteps = int(raw_input("\nInput desired number of extra steps: "))
			sys.stdout.write("\n")
			numTimeSteps += extraSteps

#------------------------------------------------#

end = time.time()
sys.stdout.write("\n")
sys.stdout.write("Simulation time = %fs\n\n" % (end - start))

Notifier.notify('The universe has been solved', title = 'Thanks to the finest minds of the 21st century...')