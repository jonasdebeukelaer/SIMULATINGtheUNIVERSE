import pmClass
import pmCore as core
import pmHelpers as helpers

import numpy as np
import math
import cmath
import random
import time
import sys

def CalculateStartingA(x, y, z):
	print "\nCalculating starting A..."
	maxList = [np.amax(np.absolute(x)), np.amax(np.absolute(y)), np.amax(np.absolute(z))]
	maxDisplacement = max(maxList)
	print maxDisplacement
	#displacementAmplitudes = x**2 + y**2 + z**2
	#print np.amax(displacementAmplitudes)
	#startingA = math.sqrt(3 * 0.5**2) / math.sqrt(np.amax(displacementAmplitudes))
	startingA = 0.5 / maxDisplacement
	if startingA < 0.001:
		sys.exit("You literally can't start soon enough... starting a calculated as %f" % startingA)
	formattedA = float("{0:.3f}".format(startingA))
	if formattedA >= 1.:
		sys.exit("Unable to execute simulating, starting a calculated as %f" % formattedA)
	else:
		print ("startingA = %f" % formattedA)
	return formattedA

def ComputeDisplacementVectors(nGrid, lBox, ns):
	xDisplacementFourier = np.zeros([nGrid, nGrid, (nGrid / 2) + 1], dtype = 'complex128')
	yDisplacementFourier = np.zeros([nGrid, nGrid, (nGrid / 2) + 1], dtype = 'complex128')
	zDisplacementFourier = np.zeros([nGrid, nGrid, (nGrid / 2) + 1], dtype = 'complex128')

	print 'Calculating displacements...'
	displacementStart = time.time()
	displacementsCalculated = 0
	for l in range(0, nGrid):
		if l < (nGrid / 2):
			kx = 2 * math.pi * l / nGrid
		else:
			kx = 2 * math.pi * (l - nGrid) / nGrid

		for m in range(0, nGrid):
			if m < (nGrid / 2):
				ky = 2 * math.pi * m / nGrid
			else:
				ky = 2 * math.pi * (m - nGrid) / nGrid

			for n in range(0, (nGrid / 2) + 1):
				kz = 2 * math.pi * n / nGrid
				displacementsCalculated += 1

				kSquare = float(kx**2 + ky**2 + kz**2)
				k       = math.sqrt(kSquare)
				if kSquare != 0:
					powerValue = 8 * math.pi**4 * 2 * 10**(-9) * ((k * nGrid) / (lBox * 0.05))**ns
					ak         = math.sqrt(powerValue) * random.gauss(0., 1.) / ((k * nGrid / lBox)**2 * math.sqrt(2))
					bk         = math.sqrt(powerValue) * random.gauss(0., 1.) / ((k * nGrid / lBox)**2 * math.sqrt(2))
				else:
					ak = 0
					bk = 0
				ck = (ak - bk * 1j)

				xDisplacementFourier[l][m][n] = ck * kx
				yDisplacementFourier[l][m][n] = ck * ky
				zDisplacementFourier[l][m][n] = ck * kz

				helpers.OutputPercentage(displacementsCalculated, nGrid**2 * ((nGrid / 2) + 1), time.time() - displacementStart)

	xDisplacementReal = np.fft.irfftn(xDisplacementFourier) * nGrid**3
	yDisplacementReal = np.fft.irfftn(yDisplacementFourier) * nGrid**3
	zDisplacementReal = np.fft.irfftn(zDisplacementFourier) * nGrid**3

	return (xDisplacementReal, yDisplacementReal, zDisplacementReal)

def InitialiseParticles(nGrid, numParticles, positionDistribution, velocityDistribution, maxVelocity, a, deltaA, lBox, ns):
	particleList = []

	if positionDistribution == pmClass.PositionDist.zeldovich:

		xDisplacements, yDisplacements, zDisplacements = ComputeDisplacementVectors(nGrid, lBox, ns)

		maxContrast = 0
		aIncrease = 0.001
		while maxContrast <= 2:
			if a == 'auto':
				a = CalculateStartingA(xDisplacements, yDisplacements, zDisplacements)
			a += aIncrease

			print 'Creating particles...'
			creationStart = time.time()
			particlesCreated = 0
			gridX = - nGrid / 2
			for i in range(0, nGrid):
				gridX += 1
				gridY = - nGrid / 2
				for j in range(0, nGrid):
					gridY += 1
					gridZ = - nGrid / 2
					for k in range(0, nGrid): 
						gridZ += 1
						particlesCreated += 1
					
						x = gridX - a * (xDisplacements[i][j][k])
						y = gridY - a * (yDisplacements[i][j][k])
						z = gridZ - a * (zDisplacements[i][j][k])

						xMomentum = - (xDisplacements[i][j][k]) * (a - (deltaA / 2))**2
						yMomentum = - (yDisplacements[i][j][k]) * (a - (deltaA / 2))**2
						zMomentum = - (zDisplacements[i][j][k]) * (a - (deltaA / 2))**2

						newParticle = pmClass.Particle([x, y, z], [xMomentum, yMomentum, zMomentum])
						core.PositionCorrect(newParticle, nGrid)
						particleList.append(newParticle)
						helpers.OutputPercentage(particlesCreated, nGrid**3, time.time() - creationStart)

			densityField = core.CalculateDensityField(nGrid, particleList)
			maxContrast = int(np.amax(densityField))
			print maxContrast, a, aIncrease


	elif positionDistribution == pmClass.PositionDist.sineWave1D:
		wavelength = nGrid
		kBox = 2.0 * math.pi / wavelength
		aCross = 10.0 * a
		waveAmplitude = 1.0 / (aCross * kBox)

		qx = - nGrid / 2
		for i in range(0, nGrid):
			qx +=  (wavelength / nGrid)

			print str(a * waveAmplitude * math.sin(kBox * qx))

			qy = - nGrid / 2
			for j in range(0, nGrid):
				qy += 1
				qz = - nGrid / 2
				for k in range(0, nGrid):
					qz += 1

					x = qx - a * waveAmplitude * math.sin(kBox * qx)
					y = qy
					z = qz
					
					xMomentum = - waveAmplitude * math.sin(kBox * qx)
					yMomentum = 0
					zMomentum = 0

					newParticle = pmClass.Particle([x, y, z], [xMomentum, yMomentum, zMomentum])
					core.PositionCorrect(newParticle, nGrid)
					particleList.append(newParticle)

	else:
		for i in range(0, numParticles):
			newParticle = pmClass.Particle([0., 0., 0.], [0., 0., 0.])
			particleList.append(newParticle)

		# Random coordinates distribution
		if positionDistribution == pmClass.PositionDist.random:
			for particle in particleList:
				x                 = random.randrange((- nGrid / 2) / 0.001, (nGrid / 2) / 0.001) * 0.001
				y                 = random.randrange((- nGrid / 2) / 0.001, (nGrid / 2) / 0.001) * 0.001
				z                 = random.randrange((- nGrid / 2) / 0.001, (nGrid / 2) / 0.001) * 0.001
				particle.position = [x, y, z]

		# Random shell of particles distribution
		elif positionDistribution == pmClass.PositionDist.randomShell:
			r = nGrid / 2
			for particle in particleList:
				theta = math.acos(random.uniform(-1., 1.))
				phi   = random.uniform(0., 2 * math.pi)
				
				x     = r * math.sin(theta) * math.cos(phi)
				y     = r * math.sin(theta) * math.sin(phi)
				z     = r * math.cos(theta)

				particle.position = [x, y, z]

		# Even shell of particles distribution
		elif positionDistribution == pmClass.PositionDist.evenShell:
			r            = nGrid / 2
			
			golden_angle = np.pi * (3 - np.sqrt(5))
			theta        = golden_angle * np.arange(numParticles)
			z            = np.linspace(r - 1.0 / numParticles, 1.0 / numParticles - r, numParticles)
			radius       = np.sqrt((r**2 - z * z))
			 
			points       = np.zeros((numParticles, 3))
			points[:,0]  = radius * np.cos(theta)
			points[:,1]  = radius * np.sin(theta)
			points[:,2]  = z

			i = 0
			for particle in particleList:
				particle.position = points[i, :]
				i += 1

		else:
			print 'Invalid position distribution selected'

	# Random velocity distribution
	if velocityDistribution == pmClass.VelocityDist.random:
		for particle in particleList:
			xVelocity                 = random.randrange(- maxVelocity / 0.001, maxVelocity / 0.001) * 0.001
			yVelocity                 = random.randrange(- maxVelocity / 0.001, maxVelocity / 0.001) * 0.001
			zVelocity                 = random.randrange(- maxVelocity / 0.001, maxVelocity / 0.001) * 0.001
			particle.halfStepMomentum = [xVelocity, yVelocity, zVelocity]

	# Zero velocity distribution
	elif velocityDistribution == pmClass.VelocityDist.zero:
		for particle in particleList:
			particle.halfStepMomentum = [0, 0, 0]

	elif velocityDistribution == pmClass.VelocityDist.zeldovich or velocityDistribution == pmClass.VelocityDist.sineWave1D:
		pointless = 0

	else:
		print 'Invalid velocity distribution selected'

	return [a, particleList]