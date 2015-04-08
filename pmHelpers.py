import math
import os
import sys

import smtplib

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
	i = (float(timeStep) / numTimeSteps) * 100
	if round(numTimeSteps) == timeStep:
		timeLeft = 0.
	else:
		timeLeft = timeElapsed * ((100. / i) - 1.)
	hours = math.floor(timeLeft / 3600.)
	minutes = math.floor((timeLeft % 3600) / 60.)
	seconds = (timeLeft % 60.) + 0.5
	sys.stdout.write("\r%.2f%% - Estimated time remaining: %02d:%02d:%02d" % (i, hours, minutes, seconds))
	sys.stdout.flush()

def SendEmail():
	sender = 'oliverclipsham@yahoo.com'
	receivers = ['oliverclipsham@yahoo.com']

	message = """From: From Person <oliverclipsham@yahoo.com>
	To: To Person <oliverclipsham@yahoo.com>
	Subject: Simulation complete!

	Simulation complete muthafucka!
	"""

	smtpObj = smtplib.SMTP('localhost')
	smtpObj.sendmail(sender, receivers, message)         
	print "Successfully sent email"