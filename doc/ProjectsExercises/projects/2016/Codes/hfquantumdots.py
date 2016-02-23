import numpy as np 
import time

# expectation value for the one body part
def onebody(i, n, ml):
	hbar = 1
	omega = 1
	return hbar*omega*(2*n[i] + abs(ml[i]) + 1)

if __name__ == '__main__':
	
	Norbitals = 2 # number of electrons
	

	""" Read quantum numbers from file """
	n = []
	ml = []
	spOrbitals = 0
	with open("qdspnumbers.dat", "r") as qnumfile:
		for line in qnumfile:
			nums = line.split()
			if len(nums) != 0:
				n.append(int(nums[0]))
				ml.append(int(nums[1]))
				spOrbitals += 1
	print spOrbitals, " spOrbitals"

	""" Read two-electron integrals from file """
	two_electron_I = np.zeros([spOrbitals, spOrbitals, spOrbitals, spOrbitals])
	with open("qdtwobody.dat", "r") as infile:
		for line in infile:
			l = line.split()
			a = int(l[0]) - 1
			b = int(l[1]) - 1
			c = int(l[2]) - 1
			d = int(l[3]) - 1
			#print a, b, c, d, float(l[4])
			two_electron_I[a][b][c][d] = float(l[4])
	#raw_input("") 

	""" Calculate one-electron integral """
	singleparticleH = np.zeros(spOrbitals)
	for i in range(spOrbitals):
		singleparticleH[i] = onebody(i, n, ml)
		print singleparticleH[i]

	""" Run HF-iterations """
        Nparticles = 2
	C = np.ones([spOrbitals, spOrbitals]) # HF coefficients
        maxHFiter = 100
        epsilon =  1.0e-10 
        diff = 1.0
	hf_count = 0
	oldenergies = np.zeros(spOrbitals)
	newenergies = np.zeros(spOrbitals)
	while diff > epsilon and hf_count < maxHFiter:
		print "############### Iteration %i ###############" % hf_count
		h = np.zeros([spOrbitals, spOrbitals])
   	        HFmatrix = np.zeros([spOrbitals,spOrbitals])		
		for alpha in range(spOrbitals):
			for beta in range(spOrbitals):
                            if alpha == beta:   HFmatrix[alpha][beta] += singleparticleH[alpha]
     		            FockTermSum = 0.0
                            for p in range(Nparticles):
				for gamma in range(spOrbitals):
                                    for delta in range(spOrbitals):
                                        FockTermSum += C[p][gamma]*C[p][delta]*(two_electron_I[alpha][gamma][beta][delta]-two_electron_I[alpha][gamma][delta][beta])
                        print " Fock term"
                        print alpha, beta, FockTermSum
                        HFmatrix[alpha][beta] += FockTermSum
  		   
                print " Hartree-Fock Matrix" 
                print HFmatrix
		spenergies, C = np.linalg.eigh(HFmatrix)
		newenergies = spenergies
		diff = sum(abs(newenergies-oldenergies))/spOrbitals
                oldenergies = newenergies
                print "Single-particle energies"
		print oldenergies
		hf_count += 1

