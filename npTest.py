import numpy as np
from numpy import random

testInput = random.rand(10, 10, 10)

intermediateStage = np.fft.rfftn(testInput)

print intermediateStage[1][1][1]