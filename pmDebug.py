import pmCore as core

import numpy as np
import math

def OutputPotentialFieldXY(potentialField, particleList, nGrid, timeStep):
	g = open("PotentialResults/potential_frame%d.3D" % (timeStep), "w")
	g.write("x y z Potential\n")
	
	for i in range(0, nGrid):
		for j in range(0, nGrid):
			for k in range(0, nGrid):
				potentialValue = potentialField[i][j][k]
				g.write("%f %f %f %f\n" % (i - (nGrid / 2), j - (nGrid / 2), k - (nGrid / 2), potentialValue))

	sliceOutput = open("PotentialResults/potential_slice%d.3D" % timeStep, "w")
	sliceOutput.write("x y z Potential\n")

	for i in range(0, nGrid):
		for j in range(0, nGrid):
			potentialValue = potentialField[i][j][0]	
			sliceOutput.write("%f %f 0 %f\n" % (i - (nGrid / 2), j - (nGrid / 2), potentialValue))

	#for particle in particleList:
	#	sliceOutput.write("%f %f %f %f\n" % (particle.position[0], particle.position[1], particle.position[2], -0.1))

	g.write("%f %f %f %f\n%f %f %f %f\n" % (nGrid / 2, nGrid / 2, nGrid / 2, 0., - nGrid / 2, - nGrid / 2, - nGrid / 2, 0.))
	g.close()

	sliceOutput.write("%f %f %f %f\n%f %f %f %f\n" % (nGrid / 2, nGrid / 2, nGrid / 2, 0., - nGrid / 2, - nGrid / 2, - nGrid / 2, 0.))
	sliceOutput.close()


def OutputTotalEnergy(particleList, potentialField, a, stepSize, nGrid):
	
	potential = 0
	for particle in particleList:
		xIndex = core.FindMeshIndex(particle.position[0], nGrid)
		yIndex = core.FindMeshIndex(particle.position[1], nGrid)
		zIndex = core.FindMeshIndex(particle.position[2], nGrid)

		potential += potentialField[xIndex][yIndex][zIndex] * particle.mass

	potentialE  = 0.5 * potential * a

	totalMoms = [(part.halfStepMomentum[0]**2 + part.halfStepMomentum[1]**2 + part.halfStepMomentum[2]**2) / part.mass for part in particleList]
	kinetic = sum(np.multiply(0.5, totalMoms))

	return [potentialE + kinetic, potentialE, kinetic]

def OutputDensityField(nGrid, densityField, frameNo):
	densityFile = open("Densityresults/density_frame%d.3D" % (frameNo), "w")
	densityFile.write("x y z Density\n")

	for i in range(0, nGrid):
		for j in range(0, nGrid):
			for k in range(0, nGrid):
				cellDensity = densityField[i][j][k]
				if cellDensity != 0:
					densityFile.write("%f %f %f %f\n" % (i-(nGrid/2-1), j-(nGrid/2-1), k-(nGrid/2-1), abs(cellDensity)))

	densityFile.close()