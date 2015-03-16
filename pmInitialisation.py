import pmClass
import pmCore as core

import numpy as np
import math
import cmath
import random



def ComputeDisplacementVectors(nGrid, lBox, a):
	xDisplacementFourier = np.zeros([nGrid, nGrid, (nGrid / 2) + 1], dtype = 'complex128')
	yDisplacementFourier = np.zeros([nGrid, nGrid, (nGrid / 2) + 1], dtype = 'complex128')
	zDisplacementFourier = np.zeros([nGrid, nGrid, (nGrid / 2) + 1], dtype = 'complex128')

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

				kSquare = float(kx**2 + ky**2 + kz**2) #* (2 * math.pi / Lbox)**2
				k       = math.sqrt(kSquare)
				if kSquare != 0:
					powerValue = 10**(-4.5) * ((k * nGrid) / (lBox * 0.05))
					ak         = math.sqrt(powerValue) * random.gauss(0., 1.) / (k * nGrid / lBox)**2
					bk         = math.sqrt(powerValue) * random.gauss(0., 1.) / (k * nGrid / lBox)**2
				else:
					ak = 0
					bk = 0
				ck = (ak - bk * 1j) / 2

				xDisplacementFourier[l][m][n] = ck * kx
				yDisplacementFourier[l][m][n] = ck * ky
				zDisplacementFourier[l][m][n] = ck * kz

	xDisplacementReal = np.fft.irfftn(xDisplacementFourier) * nGrid**3
	yDisplacementReal = np.fft.irfftn(yDisplacementFourier) * nGrid**3
	zDisplacementReal = np.fft.irfftn(zDisplacementFourier) * nGrid**3

	return (xDisplacementReal, yDisplacementReal, zDisplacementReal)

def ComputeDisplacementAlternative(nGrid, lBox, a):
	xDisplacementFourier = np.zeros([nGrid, nGrid, (nGrid / 2) + 1], dtype = 'complex128')
	yDisplacementFourier = np.zeros([nGrid, nGrid, (nGrid / 2) + 1], dtype = 'complex128')
	zDisplacementFourier = np.zeros([nGrid, nGrid, (nGrid / 2) + 1], dtype = 'complex128')

	dk = 2 * math.pi / lBox

	for nx in range(0, nGrid):
		kx = 2 * math.pi * nx / nGrid
		for ny in range(0, nGrid):
			ky = 2 * math.pi * ny / nGrid
			for nz in range(0, (nGrid / 2) + 1):
				kz = 2 * math.pi * nz / nGrid

				k          = math.sqrt((nx**2 + ny**2 + nz**2) * dk**2)
				powerValue = (2 * math.pi**2) * k * (1. / 70.)**4 * (10**7)**2 * a**2

				Ax = float(math.sqrt(-math.log(1. - random.random()) * powerValue))
				Ay = float(math.sqrt(-math.log(1. - random.random()) * powerValue))
				Az = float(math.sqrt(-math.log(1. - random.random()) * powerValue))

				thetaX = float(random.uniform(0. , 2 * math.pi))
				thetaY = float(random.uniform(0. , 2 * math.pi))
				thetaZ = float(random.uniform(0. , 2 * math.pi))

				deltaKx = Ax * cmath.exp(thetaX * 1j)
				deltaKy = Ay * cmath.exp(thetaY * 1j)
				deltaKz = Az * cmath.exp(thetaZ * 1j)

				xDisplacementFourier[nx][ny][nz] = deltaKx
				yDisplacementFourier[nx][ny][nz] = deltaKy
				zDisplacementFourier[nx][ny][nz] = deltaKz

	xDisplacementReal = np.fft.irfftn(xDisplacementFourier) * nGrid**3
	yDisplacementReal = np.fft.irfftn(yDisplacementFourier) * nGrid**3
	zDisplacementReal = np.fft.irfftn(zDisplacementFourier) * nGrid**3

	return (xDisplacementReal, yDisplacementReal, zDisplacementReal)

def InitialiseParticles(nGrid, numParticles, positionDistribution, velocityDistribution, maxVelocity, a, deltaA, lBox):
	particleList = []

	if positionDistribution == pmClass.PositionDist.zeldovich:

		displacementVectors = ComputeDisplacementVectors(nGrid, lBox, a)
		xDisplacements = (displacementVectors[0])
		yDisplacements = (displacementVectors[1])
		zDisplacements = (displacementVectors[2])

		gridX = - nGrid / 2
		for i in range(0, nGrid):
			gridX += 1
			gridY = - nGrid / 2
			for j in range(0, nGrid):
				gridY += 1
				gridZ = - nGrid / 2
				for k in range(0, nGrid): 
					gridZ += 1
				
					x = gridX - a * (xDisplacements[i][j][k])
					y = gridY - a * (yDisplacements[i][j][k])
					z = gridZ - a * (zDisplacements[i][j][k])

					xMomentum = - (a - (deltaA / 2))**2 * (xDisplacements[i][j][k])
					yMomentum = - (a - (deltaA / 2))**2 * (yDisplacements[i][j][k])
					zMomentum = - (a - (deltaA / 2))**2 * (zDisplacements[i][j][k])

					newParticle = pmClass.Particle([x, y, z], [xMomentum, yMomentum, zMomentum], 1)
					core.PositionCorrect(newParticle, nGrid)
					particleList.append(newParticle)

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

					newParticle = pmClass.Particle([x, y, z], [xMomentum, yMomentum, zMomentum], 1)
					core.PositionCorrect(newParticle, nGrid)
					particleList.append(newParticle)

	else:
		for i in range(0, numParticles):
			newParticle = pmClass.Particle([0., 0., 0.], [0., 0., 0.], 1)
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

	return particleList