import math
import numpy as np
from sympy import *
from sympy.combinatorics.named_groups import CyclicGroup
from sympy.combinatorics.perm_groups import PermutationGroup
from sympy.ntheory.modular import crt
import time
import matplotlib.pyplot as plt
from decimal import *

getcontext().prec = 1000000

def get_order(g, p):
	e = 1

	for d in range(1, p):
		if ((int(g)**d) % p) == e:
			return d

def get_inverse(g, group):
	m = group.order()

	for b in range(m):
		if (g*b) % m == 1:
			return b

def extendedGCD(a, b):
	if a == 0:
		return b, 0, 1

	gcd, x1, y1 = extendedGCD(b % a, a)

	x = y1 - (b // a) * x1
	y = x1

	return gcd, x, y


#Shanks' algorithm
def shanks_algo(g, h, group):
	p = group.order()
	n = get_order(g, p)

	m = math.floor(math.sqrt(n)) + 1

	u = pow(g, m, p)
	b = extendedGCD(g, p)[1]

	l1 = [(j, pow(u, j, p)) for j in range(m)]
	l2 = [(i, h * (b**i) % p) for i in range(m)]
	
	l1.sort(key=lambda tup: tup[1])
	l2.sort(key=lambda tup: tup[1])

	i, j = 0, 0

	while (i < m) and (j < m):
		if l1[j][1] == l2[i][1]:
			return m * l1[j][0] + l2[i][0] % p
		elif abs(l1[j][1]) > abs(l2[i][1]):
			i += 1
		else:
			j += 1


#Pohlig-Hellman algorithm
def pohlig_hellman_algo(g, h, group):

	n = group.order()
	primes = ntheory.factorint(n)

	xlist = []
	plist = []
	
	for p, e in primes.items():
		gi = (g**int((n/(p**e)))) % n
		hi = (h**int((n/(p**e)))) % n
		xi = ntheory.residue_ntheory.discrete_log(n, int(hi), int(gi))
		xlist.append(xi)
		plist.append(p**e)

	res = crt(xlist, plist)
	return res[1]

def main_calc(g, h, p):

	print('Calculating {}^x = {} % {}'.format(g, h, p))
	G = CyclicGroup(p)
	print('...')

	start2 = time.time()
	x2 = pohlig_hellman_algo(g, h, G)
	end2 = time.time()
	print('Pohlig-Hellman result: {}'.format(x2))
	print('Pohlig-Hellman time: {0:.4f}'.format(end2 - start2))

	if end2 - start2 > 0.001:
		return math.inf, end2 - start2

	start1 = time.time()
	x1 = shanks_algo(g, h, G)
	end1 = time.time()
	print('Shanks result {}'.format(x1))
	print('Shanks time: {0:.4f}'.format(end1 - start1))

	print('====================')
	return end1 - start1, end2 - start2


if __name__ == '__main__':

	a1, b1 = main_calc(6, 7531, 8101)
	a2, b2 = main_calc(7, 777, 14947)
	a3, b3 = main_calc(23, 9689, 11251)
	a4, b4 = main_calc(5701, 9227, 23719)
	a5, b5 = main_calc(397, 9173, 32341)
	a6, b6 = main_calc(1741, 36899, 13043)
	a7, b7 = main_calc(5581, 5399, 27509)
	a8, b8 = main_calc(4651, 2837, 13709)
	a9, b9 = main_calc(1277, 7951, 14843)
	a10, b10 = main_calc(7393, 8101, 20249)
	a11, b11 = main_calc(75109, 781661, 72678533)
	a12, b12 = main_calc(38569, 5525087, 51562439)
	
	
	shanks = np.mean([a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12])
	pohlig_hellman = np.mean([b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12])

	print(shanks)
	print(pohlig_hellman)


	plt.plot([1,2,3,4,5,6,7,8,9,10, 11, 12], np.log([a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12]), linestyle='-', label='Shanks', color='g')
	plt.plot([1,2,3,4,5,6,7,8,9,10, 11, 12], np.log([b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12]), linestyle='--', label='Pohlig-Hellman', color='r')
	plt.legend()
	plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
	plt.ylabel('Logarithmic time (s)')
	plt.xlabel('Discrete log g^x = h % p')
	plt.show()

	
