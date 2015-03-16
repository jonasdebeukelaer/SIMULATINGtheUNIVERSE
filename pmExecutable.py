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
random.seed(30091874) #30091875)		#
print "Done\n"

#------------INITIALISATION PARAMETERS-----------#

volume                 = [10, 10, 10]
gridResolution         = 1
Lbox 				   = 10

numParticles           = 0
positionDistribution   = pmClass.PositionDist.zeldovich
velocityDistribution   = pmClass.VelocityDist.zeldovich

preComputeGreens       = True

maxVelocity            = 1
hasCenterParticle      = False

startingA              = 0.10
maxA                   = 1.000
stepSize               = 0.001

shootEvery             = 2
outputAsSphereOnly     = False

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
#particleList = []
#particleList.append(pmClass.Particle([-5, 0, 0], [1, 0, 0], 100))
#particleList.append(pmClass.Particle([-3, 1, 0], [0, -0.2, 0], 100))
#numParticles = 2


if hasCenterParticle:
	centerParticle = pmClass.Particle([0., 24.64, 0.], [0., -1., 0.,], 20)
	particleList.append(centerParticle)
	numParticles += 1
print"Done\n"

print "Determining mesh shape..."
densityField = core.CalculateDensityField(volume, gridResolution, particleList)
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
initial.write("x y z LocalDensity\n")
for particle in particleList:
	localDensity = core.FindLocalDensity(particle, densityField, gridResolution)
	initial.write("%f %f %f %f\n" % (particle.position[0], particle.position[1], particle.position[2], math.log(localDensity)))
initial.write("%f %f %f %f\n%f %f %f %f\n" % (volume[0] / 2, volume[1] / 2, volume[2] / 2, 1., - volume[0] / 2, - volume[1] / 2, - volume[2] / 2, 1.))
initial.close()


energyFile = open("energyResults.txt", "w")
energyFile.write("a\tTotal\tpotential\tkinetic\t%% kinetic off\t%% error in energy\n")
startE = 0

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

	densityField   = core.CalculateDensityField(volume, gridResolution, particleList)
	potentialField = core.SolvePotential(densityField, a, greensFunction, preComputeGreens)

	accumulatedEnergy = 0

	for i, particle in enumerate(particleList):

		particle.acceleration      = core.CalculateParticleAcceleration(particle, potentialField, gridResolution)
		particle.halfStepMomentum += np.multiply(core.GetF(a - stepSize) * stepSize, particle.acceleration)
		localDensity               = core.FindLocalDensity(particle, densityField, gridResolution)
		particle.position         += np.multiply((a - (stepSize / 2))**(-2) * core.GetF(a - (stepSize / 2)) * stepSize, particle.halfStepMomentum)
		core.PositionCorrect(particle, volume)

		if shoot:
			if positionDistribution == pmClass.PositionDist.sineWave1D:
				f.write("%f %f %f %f\n" % (particle.position[0], particle.position[1], particle.position[2], particle.halfStepMomentum[0]))
			else:
				if outputAsSphereOnly:
					r = math.sqrt(particle.position[0]**2 + particle.position[1]**2 + particle.position[2]**2)
					if r <= volume[0]/2:
						f.write("%f %f %f %f\n" % (particle.position[0], particle.position[1], particle.position[2], math.log(localDensity)))
				else:
					f.write("%f %f %f %f\n" % (particle.position[0], particle.position[1], particle.position[2], localDensity))

	if shoot:
		if outputSystemEnergy:	
			energyResults = debug.OutputTotalEnergy(particleList, potentialField, a, stepSize, volume)
			accumulatedEnergy = energyResults[0]
			potentialE = energyResults[1]
			ke = energyResults[2]
			f.write("0 %d 0 %f\n" % (volume[2]/2, accumulatedEnergy))
			if startE == 0:
				startE = accumulatedEnergy
			energyFile.write("%f\t%f\t%f\t%f\t%f\t%f\n" % (a, accumulatedEnergy, potentialE, ke, (startE - accumulatedEnergy) / ke, (startE - accumulatedEnergy)/startE * 100))
			print "\t", accumulatedEnergy, '\t potential=', potentialE, ' ke=', ke, '*ke off=', (startE - potentialE) / ke

		if outputDensityField:
			debug.OutputDensityField(volume, densityField)

		if outputPotentialFieldXY:
			debug.OutputPotentialFieldXY(potentialField, particleList, volume, frameNo, gridResolution)

		f.write("%f %f %f %f\n%f %f %f %f\n" % (volume[0] / 2, volume[1] / 2, volume[2] / 2, 1., - volume[0] / 2, - volume[1] / 2, - volume[2] / 2, 1.))
		f.close()

	#core.OutputPowerSpectrum(densityField)
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

energyFile.close()

if outputSystemEnergy:
	print 'startE=', startE, ' endE=', accumulatedEnergy, '\t difference percent=', (accumulatedEnergy-startE)/ startE * 100

#------------------------------------------------#






Notifier.notify('The universe has been solved', title = 'Thanks to the finest minds of the 21st century...')