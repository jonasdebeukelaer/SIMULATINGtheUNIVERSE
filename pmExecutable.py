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

nGrid                  = 32
lBox 				   = 32
ns                     = 1

numParticles           = 0
positionDistribution   = pmClass.PositionDist.zeldovich
velocityDistribution   = pmClass.VelocityDist.zeldovich

usePPmethod            = False
preComputeGreens       = True

maxVelocity            = 1
hasCenterParticle      = False

startingA              = 'auto'
maxA                   = 1.000
stepSize               = 0.001

shootEvery             = 2
outputAsSphere         = False

outputPowerSpectrum    = True
powerSpectrumRes       = 0.5

sendEmail              = False

#----------------DEBUG PARAMETERS----------------#

outputPotentialFieldXY = False
outputSystemEnergy     = False
outputDensityField     = False

#------------------------------------------------#

#------------INITIALISATION FUNCTIONS------------#

particleList = []
if numParticles == 'grafic':
	graficConditions = open('p3m_32.dat', 'r')

	params     = graficConditions.readline()
	paramsList = params.split()
	startingA  = float("{0:.3f}".format(float(paramsList[6])))
	#startingA  = 1.0/float("{0:.3f}".format(float(paramsList[6])))
	factor     = startingA - (stepSize / 2)

	for pV in graficConditions:
		values = pV.split()
		particle = pmClass.Particle([float(values[0]) - 16., float(values[1]) - 16., float(values[2]) - 16.], [float(values[3]) * factor, float(values[4]) * factor, float(values[5]) * factor])
		particleList.append(particle)
else:
	print "Initialising particles..."
	startingA, particleList = initialisation.InitialiseParticles(nGrid, numParticles, positionDistribution, velocityDistribution, maxVelocity, startingA, stepSize, lBox, ns)

print '\n'
print startingA
#sys.exit('Everything good?')

if positionDistribution == pmClass.PositionDist.zeldovich or numParticles == 'grafic':
	numParticles = len(particleList)

if hasCenterParticle:
	centerParticle = pmClass.Particle([0., 24.64, 0.], [0., -1., 0.,])
	particleList.append(centerParticle)
	numParticles += 1
print"\nDone\n"

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
	initial.write("%f %f %f %f\n" % (particle.position[0], particle.position[1], particle.position[2], localDensity+1))
initial.write("%f %f %f %f\n%f %f %f %f\n" % (nGrid / 2, nGrid / 2, nGrid / 2, 1., - nGrid / 2, - nGrid / 2, - nGrid / 2, 1.))
initial.close()

print "Iterating..."
a = startingA
iterationStart = time.time()
frameNo = 1
maxFrameNo = int(round((maxA - startingA) / stepSize)) + 1

if outputSystemEnergy:
	energyFile = open("energyResults.txt", "w")
	energyFile.write("a\tTotal\tpotential\tkinetic\t%% kinetic off\t%% error in energy\n")
	startE = 0

if outputPowerSpectrum:
	initialPowerSpectrum = core.CalculatePowerSpectrum(densityField, nGrid, lBox, powerSpectrumRes)

while frameNo < maxFrameNo:
	shoot = True if (frameNo % shootEvery) == 0 else False
	if shoot:
		f = open("Results/values_frame%d.3D" % (frameNo), "w")
		f.write("x y z LocalDensity\n")
		if outputAsSphere:
			sphereF = open("Results/sphereValues_frame%d.3D" % (frameNo), "w")
			sphereF.write("x y z LocalDensity\n")

	densityField   = core.CalculateDensityField(nGrid, particleList)
	if usePPmethod == False:
		potentialField = core.SolvePotential(densityField, a, greensFunction, preComputeGreens)

	accumulatedEnergy = 0

	if frameNo == maxFrameNo/2:
		middlePowerSpectrum = core.CalculatePowerSpectrum(densityField, nGrid, lBox, powerSpectrumRes)

	for i, particle in enumerate(particleList):

		if usePPmethod == False:
			particle.acceleration  = core.CalculateParticleAcceleration(particle, potentialField)
		else:
			core.CalculateParticleAccelerationPP(particle, i, particleList, numParticles)    
		particle.halfStepMomentum += np.multiply(core.GetF(a - stepSize) * stepSize, particle.acceleration)
		particle.position         += np.multiply((a - (stepSize / 2))**(-2) * core.GetF(a - (stepSize / 2)) * stepSize, particle.halfStepMomentum)
		core.PositionCorrect(particle, nGrid)
		localDensity               = core.FindLocalDensity(particle, densityField)

		if usePPmethod == False:
			particle.accumulatedForce = [0.0, 0.0, 0.0]

		if shoot:
			f.write("%f %f %f %f\n" % (particle.position[0], particle.position[1], particle.position[2], localDensity+1))
			if outputAsSphere:
				r = math.sqrt(particle.position[0]**2 + particle.position[1]**2 + particle.position[2]**2)
				if r <= nGrid/2:
					sphereF.write("%f %f %f %f\n" % (particle.position[0], particle.position[1], particle.position[2], localDensity+1))
	if shoot:
		if outputSystemEnergy:	
			energyResults = debug.OutputTotalEnergy(particleList, potentialField, a, stepSize, nGrid)
			accumulatedEnergy = energyResults[0]
			potentialE = energyResults[1]
			ke = energyResults[2]
			f.write("0 %d 0 %f\n" % (nGrid/2, accumulatedEnergy))
			if startE == 0:
				startE = accumulatedEnergy
			energyFile.write("%f\t%f\t%f\t%f\t%f\t%f\n" % (a, accumulatedEnergy, potentialE, ke, (startE - accumulatedEnergy) / ke, (startE - accumulatedEnergy)/startE * 100))
			print "\t", accumulatedEnergy, '\t potential=', potentialE, ' ke=', ke, '*ke off=', (startE - potentialE) / ke

		if outputDensityField:
			debug.OutputDensityField(nGrid, densityField, frameNo)

		if outputPotentialFieldXY:
			debug.OutputPotentialFieldXY(potentialField, particleList, nGrid, frameNo)

		f.write("%f %f %f %f\n%f %f %f %f\n" % (nGrid / 2, nGrid / 2, nGrid / 2, 1., - nGrid / 2, - nGrid / 2, - nGrid / 2, 1.))
		f.close()

	helpers.OutputPercentage(frameNo, (maxA - startingA) / stepSize, time.time() - iterationStart)

	a       += stepSize
	frameNo += 1

if outputSystemEnergy:
	print 'startE=', startE, ' endE=', accumulatedEnergy, '\t difference percent=', (accumulatedEnergy-startE)/ startE * 100
	energyFile.close()

#------------------------------------------------#

print '\n', time.time() - iterationStart
Notifier.notify('The universe has been solved', title = 'Thanks to the finest minds of the 21st century...')
if sendEmail: 
	helpers.SendEmail()

if outputPowerSpectrum:
	finalPowerSpectrum = core.CalculatePowerSpectrum(densityField, nGrid, lBox, powerSpectrumRes)
	core.OutputPowerSpectrum(initialPowerSpectrum, middlePowerSpectrum, finalPowerSpectrum, startingA, nGrid, lBox, powerSpectrumRes)
