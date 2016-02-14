import numpy as np 
import time

# expectation value for the one body part, Harmonic oscillator in three dimensions
def onebody(i, n, ml):
	hbar = 1
	omega = 1
	return hbar*omega*(2*n[i] + l[i] + 1.5)

# get energy from eigenvalues
def get_energy(epsilon):
	return np.sum(epsilon)

# map nucleon number to state 
def state(alpha):
	return np.floor(alpha/2.)


if __name__ == '__main__':
	
	Norbitals = 2 # number of nucleons
	

	""" Read quantum numbers from file """
	n = []
	ml = []
	states = 0
	with open("qnumbers.dat", "r") as qnumfile:
		for line in qnumfile:
			nums = line.split()
			if len(nums) != 0:
				n.append(int(nums[0]))
				ml.append(int(nums[1]))
				states += 1
	print states, " states"

	""" Read two-electron integrals from file """
	two_electron_I = np.zeros([states, states, states, states])
	with open("matrixelements.dat", "r") as infile:
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
	one_electron_I = np.zeros(states)
	for i in range(states):
		one_electron_I[i] = onebody(i, n, ml)
		#print one_electron_I[i]
	#raw_input("")
	
	#startTime = time.time()
	""" Run HF-iterations """
	hf_count = 0
	E_old = 0
	E_new = 1
	C = np.ones([Norbitals, Norbitals]) # HF coefficients

	while hf_count < 5:
		print "############### Iteration %i ###############" % hf_count

		h = np.zeros([Norbitals, Norbitals])
		F = np.zeros([Norbitals, Norbitals])
		
		for i in range(Norbitals):
			h[i][i] = C[i][i]*one_electron_I[state(i)] 
			for j in range(Norbitals):
				
				# F[i][j] += C[i][j]*(two_electron_I[state(i)][state(j)][state(i)][state(j)] - \
				# 							two_electron_I[state(i)][state(j)][state(j)][state(i)])

				alpha = i
				gamma = j

				for p in range(Norbitals):
					for beta in range(Norbitals):
						for delta in range(Norbitals):
							a, b, d, g = state(alpha), state(beta), state(delta), state(gamma)

							F[alpha][gamma] += C[p][beta]*C[p][delta]* \
									(two_electron_I[a][b][g][d] - two_electron_I[a][b][d][g])



		print "F:"
		print F
		print "h:"
		print h
		
		epsilon, C = np.linalg.eigh(h + F)
		E_old = E_new
		E_new = get_energy(epsilon)

		print "new eps"
		print epsilon
		print "new C"
		print C
		print 
		print "dE = ", abs(E_new - E_old)
		print "Energy = ", E_new
		print
		print

		hf_count += 1

