import math
import os
import sys

def GetNumberOfThreads():
	user = os.getlogin()
	if user == "oliclipsham":
		threads = 8
	elif user == "jonasdebeukelaer":
		threads = 4
	else:
		print "Imposter!"
		threads = 1

	return threads

def OutputPercentage(timeStep, numTimeSteps, timeElapsed):
	i = (float(timeStep + 1) / numTimeSteps) * 100
	timeLeft = timeElapsed * ((100. / i) - 1)
	hours = math.floor(timeLeft / 3600.)
	minutes = math.floor((timeLeft % 3600) / 60.)
	seconds = (timeLeft % 60.) + 0.5
	sys.stdout.write("\r%.2f%% - Estimated time remaining: %02d:%02d:%02d" % (i, hours, minutes, seconds))
	sys.stdout.flush()