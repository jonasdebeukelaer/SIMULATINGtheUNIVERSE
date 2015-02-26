import pmClass
import pmInitialisation as initialisation
import pmCore as core
import pmDebug as debug
import pmHelpers as helpers

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
random.seed(30091874)
print "Done\n"

#------------INITIALISATION PARAMETERS-----------#

volume                 = [20, 20, 20]
gridResolution         = 0.05
Lbox 				   = 14000

numParticles           = 0
positionDistribution   = pmClass.PositionDist.sineWave1D
velocityDistribution   = pmClass.VelocityDist.sineWave1D

preComputeGreens       = True

maxVelocity            = 1
hasCenterParticle      = False

startingA              = 0.010
maxA                   = 1.000
stepSize               = 0.001

shootEvery             = 2

#----------------DEBUG PARAMETERS----------------#

outputPotentialFieldXY = False
outputSystemEnergy     = False
outputDensityField     = False

#------------------------------------------------#

#------------INITIALISATION FUNCTIONS------------#

print "Initialising particles..."
particleList = initialisation.InitialiseParticles(volume, numParticles, positionDistribution, velocityDistribution, maxVelocity, startingA, stepSize, Lbox)
if positionDistribution == pmClass.PositionDist.zeldovich:
	numParticles = len(particleList)

if hasCenterParticle:
	centerParticle = pmClass.Particle([0., 24.64, 0.], [0., -1., 0.,], 20)
	particleList.append(centerParticle)
	numParticles += 1
print"Done\n"

print "Determining mesh shape..."
densityField = core.CalculateDensityField(volume, gridResolution, particleList, False)
print "Done\n"

if preComputeGreens:
	print "Calculating Green's function..."
	greensFunction = core.CreateGreensFunction(densityField.shape)
	print "Done\n"
else:
	greensFunction = 0

if os.path.exists("Results/values_frame0.3D"):
	Notifier.notify('You should probably make them not exist...', title = 'Results still exist')
	deleteFiles = raw_input("Would you like to delete yo motherflippin results from the last test, you       simpleton? (y/n): ")
	if deleteFiles == 'y':
		fileList          = glob.glob("Results/*.3D")
		potentialFileList = glob.glob("PotentialResults/*.3D")
		densityFileList   = glob.glob("DensityResults/*.3D")
		for resultFile in fileList:
			os.remove(resultFile)
		for potentialFile in potentialFileList:
			os.remove(potentialFile)
		for densityFile in densityFileList:
			os.remove(densityFile)

#------------------------------------------------#

#-----------------ITERATION LOOP-----------------#

initial = open("Results/values_frame0.3D", "w")
initial.write("x y z MomentumMagnitude\n")
for particle in particleList:
	initial.write("%f %f %f %f\n" % (particle.position[0], particle.position[1], particle.position[2], 0))
initial.close()

print "Iterating..."
a = startingA
iterationStart = time.time()
frameNo = 1
maxFrameNo = int((maxA - startingA) / stepSize) + 1
while frameNo < maxFrameNo:

	shoot = True if (frameNo % shootEvery) == 0 else False
	if shoot:
		f = open("Results/values_frame%d.3D" % (frameNo), "w")
		f.write("x y z MomentumMagnitude\n")

	densityField   = core.CalculateDensityField(volume, gridResolution, particleList)
	potentialField = core.SolvePotential(densityField, a, greensFunction, preComputeGreens)

	accumulatedEnergy = 0

	for i, particle in enumerate(particleList):

		particle.acceleration      = core.CalculateParticleAcceleration(particle, potentialField, gridResolution)
		particle.halfStepMomentum += np.multiply(core.GetF(a - stepSize) * stepSize, particle.acceleration)
		momentumMagnitude          = math.sqrt(particle.halfStepMomentum[0]**2 + particle.halfStepMomentum[1]**2 + particle.halfStepMomentum[2]**2)
		particle.position         += np.multiply((a - (stepSize / 2))**(-2) * core.GetF(a - (stepSize / 2)) * stepSize, particle.halfStepMomentum)
		core.PositionCorrect(particle, volume)

		if shoot:
			f.write("%f %f %f %f\n" % (particle.position[0], particle.position[1], particle.position[2], momentumMagnitude))
	
	if shoot:
		if outputSystemEnergy:	
			accumulatedEnergy = debug.OutputTotalEnergy(particleList, potentialField, a, stepSize, volume)
			f.write("0 %d 0 %f\n" % (volume[2]/2, accumulatedEnergy))
			print "\t", accumulatedEnergy

		if outputDensityField:
			debug.OutputDensityField(volume, densityField)

		if outputPotentialFieldXY:
			debug.OutputPotentialFieldXY(potentialField, particleList, volume, frameNo, gridResolution)

		f.write("%f %f %f %f\n%f %f %f %f\n" % (volume[0] / 2, volume[1] / 2, volume[2] / 2, 0., - volume[0] / 2, - volume[1] / 2, - volume[2] / 2, 0.))
		f.close()

	helpers.OutputPercentage(frameNo, (maxA - startingA) / stepSize, time.time()-start)

	a       += stepSize
	frameNo += 1
	if frameNo == maxFrameNo:

		end = time.time()
		sys.stdout.write("\n")
		sys.stdout.write("Simulation time = %fs\n\n" % (end - iterationStart))

		Notifier.notify('a = %.3f reached' % (a), title = 'User input required')
		moreSteps = raw_input("\n\nPlease check your VisIt output... would you like to increase max A (currently at a = %.3f)? (y/n): " % (a))
		if moreSteps == 'y':
			maxA = float(raw_input("\nInput new max a: "))
			maxFrameNo = int(float(maxA - startingA) / float(stepSize))
			sys.stdout.write("\n")

#------------------------------------------------#


Notifier.notify('The universe has been solved', title = 'Thanks to the finest minds of the 21st century...')