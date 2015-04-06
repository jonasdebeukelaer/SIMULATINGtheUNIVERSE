import numpy as np
import matplotlib.pyplot as plt
import statistics
import math

processingResults = []

source = open('timingResults.txt', 'r')
for line in source:
	values = line.split()
	appended = False

	for entry in processingResults:
		if values[0] == entry[0] and values[1] == entry[1]:
			entry.append(values[2])
			appended = True

	if appended == False:
		processingResults.append(values)
		appended = True

ppX     = []
ppY     = []
ppError = []
pmX     = []
pmY     = []
pmError = []

for dataSet in processingResults:
	timings  = map(float, dataSet[2:])
	discount = max(timings)
	valid    = []
	for time in timings:
		if time != discount:
			valid.append(time)

	average = statistics.mean(valid)
	error   = statistics.stdev(valid)

	if dataSet[1] == 'True':
		ppX.append(float(dataSet[0]))
		ppY.append(average)
		ppError.append(error)
	elif dataSet[1] == 'False':
		pmX.append(float(dataSet[0]))
		pmY.append(average)
		pmError.append(error)

ppSlope, ppIntercept = np.polyfit(np.log(ppX), np.log(ppY), 1)
pmSlope, pmIntercept = np.polyfit(np.log(pmX[2:]), np.log(pmY[2:]), 1)

plt.plot(ppX, ppY, 'k-o', pmX, pmY, 'k--o', linewidth=1.0)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$N_g$', fontsize=20.0, y=1.01)
plt.ylabel(r'Time Elapsed / s')

plt.title(r'Execution Time for One Time Step vs. $N_g$', y=1.01)
plt.legend(['Direct Method', 'Particle Mesh'], loc=4)
plt.text(2, 100, r'$\frac{d(log(t))}{d(log(N_g))}=%.2f$' % ppSlope, fontsize = 20)
plt.text(55, 1, r'$\frac{d(log(t))}{d(log(N_g))}=%.2f$' % pmSlope, fontsize = 20)
plt.grid(True)

plt.gcf().subplots_adjust(bottom=0.15)

plt.show()