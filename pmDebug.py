import pmCore as core

import numpy as np
import math

def OutputPotentialFieldXY(potentialField, particleList, volume, timeStep, gridResolution):
	g = open("PotentialResults/potential_frame%d.3D" % (timeStep), "w")
	g.write("x y z Potential\n")
	
	for i in range(0, int(volume[0]/gridResolution)):
		for j in range(0, int(volume[1]/gridResolution)):
			for k in range(0, int(volume[2]/gridResolution)):
				potentialValue = potentialField[i][j][k]
				if i % 2 == 0 and j % 2 == 0 and k % 2 ==0:
					g.write("%f %f %f %f\n" % (i*gridResolution-((volume[0]/2)), j*gridResolution-((volume[1]/2)), k*gridResolution-((volume[2]/2)-1), potentialValue))

	sliceOutput = open("PotentialResults/potential_slice%d.3D" % timeStep, "w")
	sliceOutput.write("x y z Potential\n")

	for i in range(0, int(volume[0]/gridResolution)):
		for j in range(0, int(volume[1]/gridResolution)):
			potentialValue = potentialField[i][j][0]
			if i % 4 == 0 and j % 4 == 0:	
				sliceOutput.write("%f %f 0 %f\n" % (i*gridResolution-((volume[0]/2)), j*gridResolution-((volume[1]/2)), potentialValue))


	for particle in particleList:
		sliceOutput.write("%f %f %f %f\n" % (particle.position[0], particle.position[1], particle.position[2], -0.1))

	g.write("%f %f %f %f\n%f %f %f %f\n" % (volume[0] / 2, volume[1] / 2, volume[2] / 2, 0., - volume[0] / 2, - volume[1] / 2, - volume[2] / 2, 0.))
	g.close()

	sliceOutput.write("%f %f %f %f\n%f %f %f %f\n" % (volume[0] / 2, volume[1] / 2, volume[2] / 2, 0., - volume[0] / 2, - volume[1] / 2, - volume[2] / 2, 0.))
	sliceOutput.close()


def OutputTotalEnergy(particleList, potentialField, a, stepSize, volume):
	
	potential = 0
	for particle in particleList:
		xIndex = core.FindMeshIndex(particle.position[0], 1, volume[0])
		yIndex = core.FindMeshIndex(particle.position[1], 1, volume[1])
		zIndex = core.FindMeshIndex(particle.position[2], 1, volume[2])

		potential += potentialField[xIndex][yIndex][zIndex] * particle.mass

	potentialE  = 0.5 * potential * a

	totalMoms = [(part.halfStepMomentum[0]**2 + part.halfStepMomentum[1]**2 + part.halfStepMomentum[2]**2) / part.mass for part in particleList]
	kinetic = sum(np.multiply(0.5, totalMoms))

	return [potentialE + kinetic, potentialE, kinetic]

def OutputDensityField(volume, densityField):
	densityFile = open("Densityresults/density_frame%d.3D" % (frameNo), "w")
	densityFile.write("x y z Density\n")

	for i in range(0, volume[0]):
		for j in range(0, volume[1]):
			for k in range(0, volume[2]):
				cellDensity = densityField[i][j][k]
				if cellDensity != 0:
					cellDensity -= 990 if cellDensity >= 1000 else cellDensity 
					if cellDensity != 0:
						densityFile.write("%f %f %f %f\n" % (i-(volume[0]/2-1), j-(volume[1]/2-1), k-(volume[2]/2-1), math.log(abs(cellDensity))))

	densityFile.close()