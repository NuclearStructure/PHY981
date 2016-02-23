import numpy as np 
import time

# expectation value for the one body part, Harmonic oscillator in three dimensions
def onebody(i, n, l):
	homega = 10
	return homega*(2*n[i] + l[i] + 1.5)

if __name__ == '__main__':
	
	""" Read quantum numbers from file """
        index = []
	n = []
        l = []
        j = []	
        mj = []
        tz = []
	spOrbitals = 0
	with open("nucleispnumbers.dat", "r") as qnumfile:
		for line in qnumfile:
			nums = line.split()
			if len(nums) != 0:
				index.append(int(nums[0]))
				n.append(int(nums[1]))
				l.append(int(nums[2]))
				j.append(int(nums[3]))
				mj.append(int(nums[4]))
				tz.append(int(nums[5]))
				spOrbitals += 1

	""" Read two-nucleon interaction elements (integrals) from file """
	nninteraction = np.zeros([spOrbitals, spOrbitals, spOrbitals, spOrbitals])
	with open("nucleitwobody.dat", "r") as infile:
		for line in infile:
			number = line.split()
			a = int(number[0]) - 1
			b = int(number[1]) - 1
			c = int(number[2]) - 1
			d = int(number[3]) - 1
			#print a, b, c, d, float(l[4])
			nninteraction[a][b][c][d] = float(number[4])
	""" Set up single-particle integral """
	singleparticleH = np.zeros(spOrbitals)
	for i in range(spOrbitals):
		singleparticleH[i] = onebody(i, n, l)
	
	""" Compute binding energy """
        spOrbitals = 80
        """ Computing the binding energy without Hartree-Fock """
        sum = 0.0
        for i in range(spOrbitals):
            for j in range(spOrbitals):
                sum += nninteraction[i][j][i][j]
#            energy += sum # + singleparticleH[i]

        print " Energy = "
        print sum*0.5            
