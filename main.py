import numPy as np
import cmath

class Particle:
	def __init__(self, position, velocity):
		self.position = position
		self.velocity = velocity



def CalcForce(pos1[3], pos2[3]):
	separationSquared = sum((pos1 - pos2)^2)
	force = 1 / (separationSquared)



	return force

