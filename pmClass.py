from enum import Enum

class Particle:
	def __init__(self, position, halfStepMomentum, mass):
		self.position         = position
		self.halfStepMomentum = halfStepMomentum
		self.mass             = mass
		self.acceleration     = [0.0, 0.0, 0.0]

class PositionDist(Enum):
	random      = 1
	randomShell = 2
	evenShell   = 3
	zeldovich   = 4
	sineWave1D  = 5

class VelocityDist(Enum):
	random    	= 1
	zero      	= 2
	zeldovich 	= 3
	sineWave1D	= 4