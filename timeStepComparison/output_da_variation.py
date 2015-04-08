import matplotlib.pyplot as plt
import numpy as np

rawdata = open('rawDaPowerSpectrumResultsngrid32.txt', 'r')

titles = rawdata.readline()

powerSpectra  = list()
x = list()


for index, i in enumerate(rawdata):
	line = i.split()
	line = line[1:]
	powerSpectra.append(line)
	x.append((float(index+1) + 0.5)*0.5)

powerSpectra = np.transpose(powerSpectra)

fig = plt.figure()

ax = fig.add_subplot(111)
#fig.subplots_adjust(top=0.85)
ax.set_title('Effect of Time Step Size on Final $P(k)$ ($a=1$)', y=1.01)

ax.set_xlabel('$k /Mpc^{-1}$', y=1.01)
ax.set_ylabel('$P(k)$ /$10^3$')


for spectrum in powerSpectra:
	for i, value in enumerate(spectrum):
		spectrum[i] = float(value) / 1000
	plt.plot(x, spectrum)

ax.legend(['$a=0.0005$', '$a=0.001$', '$a=0.002$', '$a=0.005$', '$a=0.01$', '$a=0.02$', '$a=0.05$', '$a=0.1$'], loc=2, fontsize=17)
plt.gcf().subplots_adjust(bottom=0.15)
plt.show()

rawdata.close()


#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

scaledata = open('dapowerspectrumResults.txt')
titles = scaledata.readline()

x = list()
rCorrelation = list()

for index, i in enumerate(scaledata):
	line = i.split()
	x.append(line[0])
	rCorrelation.append(line[1])

fig = plt.figure()

ax = fig.add_subplot(111)
#fig.subplots_adjust(top=0.85)
ax.set_title('r Correlation Scaling With Time Step Size')

ax.set_xlabel('Time Step')
ax.set_ylabel('r Correlation')

#ax.axis([0, 0.1, 0.99, 1])
ax.set_xscale('log')

plt.plot(x, rCorrelation, 'bo', x, rCorrelation)
plt.show()
