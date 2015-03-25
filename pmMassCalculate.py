nGrid = float(raw_input("Please input nGrid: "))
lBox  = float(raw_input("Please input lBox: "))

h = 0.7

totalMass = h**2 * 2.79 * 10**11 * lBox**3
particleMass = totalMass / (nGrid**3)

print ("The mass of each individual particle is %.2e solar masses" % particleMass)