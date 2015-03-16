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

print "Seeding..."
random.seed(30091874)
print "Done\n"

#------------INITIALISATION PARAMETERS-----------#

nGrid                  = 40
lBox 				   = 100

numParticles           = 0
positionDistribution   = pmClass.PositionDist.zeldovich
velocityDistribution   = pmClass.VelocityDist.zeldovich

preComputeGreens       = True

maxVelocity            = 1
hasCenterParticle      = False

startingA              = 0.100
maxA                   = 1.000
stepSize               = 0.001

shootEvery             = 2
outputAsSphere         = False

#----------------DEBUG PARAMETERS----------------#

outputPotentialFieldXY = False
outputSystemEnergy     = False
outputDensityField     = False

#------------------------------------------------#

#------------INITIALISATION FUNCTIONS------------#

print "Initialising particles..."
particleList = initialisation.InitialiseParticles(nGrid, numParticles, positionDistribution, velocityDistribution, maxVelocity, startingA, stepSize, lBox)
if positionDistribution == pmClass.PositionDist.zeldovich:
	numParticles = len(particleList)

if hasCenterParticle:
	centerParticle = pmClass.Particle([0., 24.64, 0.], [0., -1., 0.,], 20)
	particleList.append(centerParticle)
	numParticles += 1
print"Done\n"

print "Determining mesh shape..."
densityField = core.CalculateDensityField(nGrid, particleList)
print "Done\n"

if preComputeGreens:
	print "Calculating Green's function..."
	greensFunction = core.CreateGreensFunction(nGrid)
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
initial.write("x y z LocalDensity\n")
for particle in particleList:
	localDensity = core.FindLocalDensity(particle, densityField)
	initial.write("%f %f %f %f\n" % (particle.position[0], particle.position[1], particle.position[2], math.log(localDensity)))
initial.write("%f %f %f %f\n%f %f %f %f\n" % (nGrid / 2, nGrid / 2, nGrid / 2, 0., - nGrid / 2, - nGrid / 2, - nGrid / 2, 0.))
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
		f.write("x y z LocalDensity\n")
		if outputAsSphere:
			sphereF = open("Results/sphereValues_frame%d.3D" % (frameNo), "w")
			sphereF.write("x y z LocalDensity\n")

	densityField   = core.CalculateDensityField(nGrid, particleList)
	potentialField = core.SolvePotential(densityField, a, greensFunction, preComputeGreens)

	accumulatedEnergy = 0

	for i, particle in enumerate(particleList):

		particle.acceleration      = core.CalculateParticleAcceleration(particle, potentialField)
		particle.halfStepMomentum += np.multiply(core.GetF(a - stepSize) * stepSize, particle.acceleration)
		localDensity               = core.FindLocalDensity(particle, densityField)
		particle.position         += np.multiply((a - (stepSize / 2))**(-2) * core.GetF(a - (stepSize / 2)) * stepSize, particle.halfStepMomentum)
		core.PositionCorrect(particle, nGrid)

		if shoot:
			f.write("%f %f %f %f\n" % (particle.position[0], particle.position[1], particle.position[2], math.log(localDensity)))
			if outputAsSphere:
				r = math.sqrt(particle.position[0]**2 + particle.position[1]**2 + particle.position[2]**2)
				if r <= volume[0]/2:
					sphereF.write("%f %f %f %f\n" % (particle.position[0], particle.position[1], particle.position[2], math.log(localDensity)))
				
	if shoot:
		if outputSystemEnergy:	
			accumulatedEnergy = debug.OutputTotalEnergy(particleList, potentialField, a, stepSize, nGrid)
			f.write("0 %d 0 %f\n" % (nGrid/2, accumulatedEnergy))
			print "\t", accumulatedEnergy

		if outputDensityField:
			debug.OutputDensityField(nGrid, densityField, frameNo)

		if outputPotentialFieldXY:
			debug.OutputPotentialFieldXY(potentialField, particleList, nGrid, frameNo, gridResolution)

		f.write("%f %f %f %f\n%f %f %f %f\n" % (nGrid / 2, nGrid / 2, nGrid / 2, 0., - nGrid / 2, - nGrid / 2, - nGrid / 2, 0.))
		f.close()

	helpers.OutputPercentage(frameNo, (maxA - startingA) / stepSize, time.time() - iterationStart)

	a       += stepSize
	frameNo += 1

#------------------------------------------------#

Notifier.notify('The universe has been solved', title = 'Thanks to the finest minds of the 21st century...')