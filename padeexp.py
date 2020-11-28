from math import exp, ceil, factorial as fac
from sys import argv

# calculating e^x using Padé approximants

# https://mathoverflow.net/questions/41226/pade-approximant-to-exponential-function

def Npq_exp(p, q):
	coeffs = [ 0 ] * (p + 1)

	for j in range(p + 1):
		n = fac(p + q - j) * fac(p)
		d = fac(p + q) * fac(j) * fac(p - j)
		coeffs[j] = n / d

	return coeffs

def Dpq_exp(p, q):
	coeffs = [ 0 ] * (q + 1)

	for j in range(q + 1):
		n = fac(p + q - j) * fac(q)
		d = fac(p + q) * fac(j) * fac(q - j)
		coeffs[j] = n / d

	return coeffs


# https://en.wikipedia.org/wiki/Horner%27s_method

def horner(x, coeffs):
	p = 0

	for c in coeffs[::-1]:
		p = c + x * p

	return p

def exp_pade(x, Npq, Dpq):
	# deal with negative argument

	if x < 0:
		return 1 / exp_pade(-x, Npq, Dpq)

	# Range reduction logic stolen from
	# https://www.pseudorandom.com/implementing-exp

	log2 = 0.6931471805599453
	k = int(x / log2)
	p = 1 << k
	r = x - (k * log2)

	N = horner(r, Npq)
	D = horner(-r, Dpq)

	print("  ", r, k, p)


	return p * (N / D)

P = 4
Q = 4

NPQ = Npq_exp(P, Q)
DPQ = Dpq_exp(P, Q)

x = 3.14159265359

print(NPQ)
print(DPQ)

epx = exp_pade(x, NPQ, DPQ)
ex = exp(x)

print(x, epx, ex, abs(ex - epx) / ex)
