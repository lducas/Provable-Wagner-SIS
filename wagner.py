from math import sqrt, log, pi, exp
from scipy.special import erf, comb

def find_w(d, log2_N):
	for w in range(d):
		S = 2**w * comb(d, w)
		if log2_N < log(S, 2):
			return w
	raise ValueError("No solution")

def log2_success_one_sample(d, s, b):
	x = erf(b/(sqrt(2)*s))
	return d * log(x, 2)

def success(n, m, q, b, log2_N, quantizing=True):
	""" 
	Return None if the success probability of Wagner's algorithm for SIS^∞ with parameters n,m,q,b 
	and with N samples is less than 1/2, following the heuristic analysis of section. 
	Otherwise, return some details about the attack parameters. 
	""" 
	assert(m > n)
	assert(2*b < q)
    
    # Initializing sparse ternary vectors 
	w = find_w(m-n, log2_N) # Initial Hamming weight 
	s_0 = sqrt(w / (m-n)) # Initial standard deviation sigma_0 

	s_i = s_0
	n_ = n
	Q = sqrt(1./(2 * pi * exp(1))) if quantizing else sqrt(1./12) 
	i = 0
	while n_ > 0:
		i += 1 
		p_i = Q * q / s_i # Rounding or quantization method determines the factor Q  
		if p_i < 1:
			return None

		b_i = log2_N / log(p_i, 2) # In this estimation, p_i and b_i may not be integers as they should (see 'Fractional parameters')
		n_ -= b_i
		s_i *= sqrt(2) 

        # The attack is considered successful if N*p>1/2 (see 'Central Gaussian Heuristic')
        # It aborts as soon as this condition is reached (see 'Early Abort')
		log2_proba_gauss = log2_success_one_sample(m - n_, s_i, b)
		log2_proba_unif = n_ * log(2*b/q, 2)
		log2_proba = log2_proba_gauss + log2_proba_unif
		if log2_N + log2_proba > -1: 
			return {"w": w, "r": i, "s_0": s_0, "s_r": s_i, "l": n_, "log2_pps": log2_proba, "log2_N": log2_N}  

def find_attack(n, m, q, b, quantizing=True, step=.1):
	log2_N = 5
	while True:
		res = success(n, m, q, b, log2_N, quantizing=quantizing)
		if res is not None:
			return res
		log2_N += step

# Parameters of SIS instances in Dilithium: Table 1 in our section 'Concrete Analysis Against Dilithium'
def SISparameters(q, gamma1, gamma2, betaD, d, tau, k, l):
    """ 
	Given Dilithium parameters, return the corresponding parameters of the corresponding SIS instance.  
	"""
    n,m = 256*k, 256*(k+l+1) # The MSIS ring in Dilithium is of dimension 256 
    zeta = max(gamma1-betaD, 2*gamma2+1+2**(d-1)*tau) # Equation (7) in [1]
    zeta_= max(2*(gamma1-betaD),4*gamma2+2) # Equation (8) in [1]
    beta = int(min(zeta, zeta_))
    return n,m,q,beta

# Attack parameters for NIST security level 2, 3, and 5: Table 2 in our section 'Concrete Analysis Against Dilithium'
print("DILITHIUM 2")
q = 8380417; gamma1 = 2**17; gamma2 = (q-1)/88; betaD = 78;  d = 13; tau = 39; k = 4; l = 4 # Notation from [1]
n,m,q,beta = SISparameters(q, gamma1, gamma2, betaD, d, tau, k, l)
print("n = 256*", n//256, "; m = 256*", m//256, "; q =", q, "; beta =", beta)
#print(find_attack(n,m,q,beta, quantizing=False))
print(find_attack(n,m,q,beta, quantizing=True))

print("\nDILITHIUM 3")
q = 8380417 ; gamma1 = 2**19; gamma2 = (q-1)/32; betaD = 196; d = 13; tau = 49; k = 6; l = 5
n,m,q,beta = SISparameters(q, gamma1, gamma2, betaD, d, tau, k, l)
print("n = 256*", n//256, "; m = 256*", m//256, "; q =", q, "; beta =", beta)
#print(find_attack(n,m,q,beta, quantizing=False))
print(find_attack(n,m,q,beta, quantizing=True))

print("\nDILITHIUM 5")
q = 8380417 ; gamma1 = 2**19; gamma2 = (q-1)/32; betaD = 120; d = 13; tau = 60; k = 8; l = 7
n,m,q,beta = SISparameters(q, gamma1, gamma2, betaD, d, tau, k, l)
print("n = 256*", n//256, "; m = 256*", m//256, "; q =", q, "; beta =", beta)
#print(find_attack(n,m,q,beta, quantizing=False))
print(find_attack(n,m,q,beta, quantizing=True))

print("\nWhen Would Wagner Shine?")
n = 500 ; m = 600 ; q = 1000 ; beta = q//4
print("n =", n, "; m =", m, "; q =", q, "; beta =", beta)
#print(find_attack(n,m,q,beta, quantizing=False))
print(find_attack(n,m,q,beta, quantizing=True))

# References: 
# [1]: CRYSTALS-Dilithium specification https://pq-crystals.org/dilithium/data/dilithium-specification-round3-20210208.pdf