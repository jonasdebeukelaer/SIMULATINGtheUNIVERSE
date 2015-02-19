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

volume                 = [20, 20, 20]
gridResolution         = 1

numParticles           = 0
positionDistribution   = pm.PositionDist.zeldovich
velocityDistribution   = pm.VelocityDist.zeldovich

maxVelocity            = 1
hasCenterParticle      = False

startingA              = 0.1000
maxA                   = 1.0000
stepSize               = 0.0001

shootEvery             = 50

outputPotentialFieldXY = False
outputSystemEnergy     = False

#------------------------------------------------#

#------------INITIALISATION FUNCTIONS------------#

print "Initialising particles..."
particleList = pm.InitialiseParticles(volume, gridResolution, numParticles, positionDistribution, velocityDistribution, maxVelocity, startingA, stepSize)
if positionDistribution == pm.PositionDist.zeldovich:
	numParticles = len(particleList)

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
	print '\n'
	if deleteFiles == 'y':
		fileList = glob.glob("Results/*.3D")
		potentialFileList = glob.glob("PotentialResults/*.3D")
		for resultFile in fileList:
			os.remove(resultFile)
		for potentialFile in potentialFileList:
			os.remove(potentialFile)

#------------------------------------------------#

#-----------------ITERATION LOOP-----------------#

initial = open("Results/values_initialframe.3D", "w")
initial.write("x y z MomentumMagnitude\n")
for particle in particleList:
	initial.write("%f %f %f %f\n" % (particle.position[0], particle.position[1], particle.position[2], 0))
initial.close()

print "Iterating..."
a = startingA
frameNo = 0
maxFrameNo = int((maxA - startingA) / stepSize)
while frameNo < maxFrameNo:

	shoot = True if (frameNo % shootEvery) == 0 else False
	if shoot:
		f = open("Results/values_frame%d.3D" % (frameNo), "w")
		f.write("x y z MomentumMagnitude\n")

	densityField   = pm.CalculateDensityField(volume, gridResolution, particleList)
	potentialField = pm.SolvePotential(densityField, greensFunction, a)

	accumulatedEnergy = 0

	for i, particle in enumerate(particleList):

		particle.acceleration      = pm.CalculateParticleAcceleration(particle, potentialField, gridResolution)
		particle.halfStepMomentum += np.multiply(pm.GetF(a - stepSize) * stepSize, particle.acceleration)
		momentumMagnitude          = math.sqrt(particle.halfStepMomentum[0]**2 + particle.halfStepMomentum[1]**2 + particle.halfStepMomentum[2]**2)
		particle.position         += np.multiply((a - (stepSize / 2))**(-2) * pm.GetF(a - (stepSize / 2)) * stepSize, particle.halfStepMomentum)
		pm.PositionCorrect(particle, volume)

		if shoot:
			f.write("%f %f %f %f\n" % (particle.position[0], particle.position[1], particle.position[2], momentumMagnitude))

			if outputSystemEnergy:
				accumulatedEnergy += pm.OutputTotalEnergy(i, particle, particleList, momentumMagnitude, a)

	
	if shoot:	
		if outputSystemEnergy:	
			f.write("0 %d 0 %f\n" % (volume[2]/2, accumulatedEnergy))
			print "\t", accumulatedEnergy

		f.write("%f %f %f %f\n%f %f %f %f\n" % (volume[0] / 2, volume[1] / 2, volume[2] / 2, 0., - volume[0] / 2, - volume[1] / 2, - volume[2] / 2, 0.))
		f.close()

		if outputPotentialFieldXY:
			pm.OutputPotentialFieldXY(potentialField, particleList, volume, frameNo, gridResolution)

	pm.OutputPercentage(frameNo, (maxA - startingA) / stepSize, time.time() - start)

	a       += stepSize
	frameNo += 1
	if frameNo == maxFrameNo:
		Notifier.notify('a = %.3f reached' % (a), title = 'User input required')
		moreSteps = raw_input("\n\nPlease check your VisIt output... would you like to increase max A (currently at a = %.3f)? (y/n): " % (a))
		if moreSteps == 'y':
			maxA = float(raw_input("\nInput new max a: "))
			maxFrameNo = int(float(maxA - startingA) / float(stepSize))
			sys.stdout.write("\n")

#------------------------------------------------#

end = time.time()
sys.stdout.write("\n")
sys.stdout.write("Simulation time = %fs\n\n" % (end - start))

Notifier.notify('The universe has been solved', title = 'Thanks to the finest minds of the 21st century...')