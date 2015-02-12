import pmHeader as pm

import numpy as np
import math
import os
import sys
import random
import time
import pyfftw
import glob
from pync import Notifier

start = time.time()
print "Seeding..."
random.seed(89321)
print "Done\n"

#------------INITIALISATION PARAMETERS-----------#

volume                 = [10, 10, 10]
gridResolution         = 1

numParticles           = 500
positionDistribution   = pm.PositionDist.zeldovich
velocityDistribution   = pm.VelocityDist.zero

maxVelocity            = 1
hasCenterParticle      = False

startingA              = 0.1
maxA                   = 1
stepSize               = 0.01
shootEvery             = 100

outputPotentialFieldXY = False
outputSystemEnergy     = False

#------------------------------------------------#

#------------INITIALISATION FUNCTIONS------------#

print "Initialising particles..."
#particleList = pm.InitialiseParticles(volume, initialisationResolution, numParticles, positionDistribution, velocityDistribution, maxVelocity)
particleList = []
particle2 = pm.Particle([-2., 0., 1.], [0., 0., 0.], 4)
particleList.append(particle2)

particle1 = pm.Particle([2., 0., -1.], [0., 0., 0.], 4)
particleList.append(particle1)


if hasCenterParticle:
	centerParticle = pm.Particle([0., 24.64, 0.], [0., -1., 0.,], 20)
	particleList.append(centerParticle)
	numParticles += 1
print"Done\n"

print "Determining mesh shape..."
densityField = pm.CalculateDensityField(volume, gridResolution, particleList, False)
print "Done\n"

print "Calculating Green's function..."
greensFunction = pm.CreateGreensFunction(densityField.shape)
print "Done\n"

if os.path.exists("Results/values_frame0.3D"):
	Notifier.notify('You should probably make them not exist...', title = 'Results still exist')
	deleteFiles = raw_input("Would you like to delete yo motherflippin results from the last test, you       simpleton? (y/n): ")
	if deleteFiles == 'y':
		fileList = glob.glob("Results/*.3D")
		potentialFileList = glob.glob("PotentialResults/*.3D")
		for resultFile in fileList:
			os.remove(resultFile)
		for potentialFile in potentialFileList:
			os.remove(potentialFile)

#------------------------------------------------#

#-----------------ITERATION LOOP-----------------#

print "\nIterating..."
a = startingA
while a < maxA:

	shoot = True if (timeStep % shootEvery) == 0 else False
	if shoot:
		f = open("Results/values_frame%d.3D" % (timeStep), "w")
		f.write("x y z VelocityMagnitude\n")

	densityField   = pm.CalculateDensityField(volume, gridResolution, particleList)
	potentialField = pm.SolvePotential(densityField, greensFunction)

	accumulatedEnergy = 0

	for i, particle in enumerate(particleList):

		particleAcceleration = pm.CalculateParticleAcceleration(particle, potentialField, gridResolution)
		if timeStep == 0:
			particle.acceleration = particleAcceleration

		particle.position     += np.multiply(timeStepSize, particle.velocity) + np.multiply(0.5 * timeStepSize**2, particle.acceleration)
		pm.PositionCorrect(particle, volume)
		particle.velocity     += np.multiply((0.5 * timeStepSize), (np.add(particle.acceleration, particleAcceleration)))
		velocityMagnitude      = (((particle.velocity[0])**2 + (particle.velocity[1])**2 + (particle.velocity[2])**2))**0.5
		particle.acceleration  = particleAcceleration

		if shoot:
			f.write("%f %f %f %f\n" % (particle.position[0], particle.position[1], particle.position[2], 0))

			if outputSystemEnergy:
				accumulatedEnergy += pm.OutputTotalEnergy(i, particle, particleList, velocityMagnitude)

	
	if shoot:		
		f.write("0 %d 0 %f\n" % (volume[2]/2, accumulatedEnergy))
		print "\t", accumulatedEnergy

		f.write("%f %f %f %f\n%f %f %f %f\n" % (volume[0] / 2, volume[1] / 2, volume[2] / 2, 0., - volume[0] / 2, - volume[1] / 2, - volume[2] / 2, 0.))
		f.close()

		if outputPotentialFieldXY:
			pm.OutputPotentialFieldXY(potentialField, particleList, volume, timeStep, gridResolution)

	pm.OutputPercentage(timeStep, numTimeSteps)

	a += stepSize

	if timeStep == numTimeSteps:
		Notifier.notify('%d time steps complete' % (numTimeSteps), title = 'User input required')
		moreSteps = raw_input("\n\nPlease check your VisIt output... would you like to add more time steps (currently %d run)? (y/n): " % (numTimeSteps))
		if moreSteps == 'y':
			extraSteps = int(raw_input("\nInput desired number of extra steps: "))
			sys.stdout.write("\n")
			numTimeSteps += extraSteps

#------------------------------------------------#

end = time.time()
sys.stdout.write("\n")
sys.stdout.write("Simulation time = %fs\n\n" % (end - start))

Notifier.notify('The universe has been solved', title = 'Thanks to the finest minds of the 21st century...')