#!/home/ubeuser/anaconda3/bin/python

#####################################################################
# WEEK 5 - LAB
#####################################################################

import math
import numpy as np
#import shipy as sp
#import matplotlib as plt
from qutip import *


qb0 = basis(2,0)
qb1 = basis(2,1)

qb00 = tensor(qb0, qb0)
qb01 = tensor(qb0, qb1)
qb10 = tensor(qb1, qb0)
qb11 = tensor(qb1, qb1)

qb000 = tensor(qb0, qb0, qb0)
qb001 = tensor(qb0, qb0, qb1)
qb010 = tensor(qb0, qb1, qb0)
qb011 = tensor(qb0, qb1, qb1)
qb100 = tensor(qb1, qb0, qb0)
qb101 = tensor(qb1, qb0, qb1)
qb110 = tensor(qb1, qb1, qb0)
qb111 = tensor(qb1, qb1, qb1)

qbHp = (1.0 / math.sqrt(2)) * (qb0 + qb1)
qbHm = (1.0 / math.sqrt(2)) * (qb0 - qb1)

qbY0 = (1.0 / math.sqrt(2)) * (qb0 + 1j * qb1)
qbY1 = (1.0 / math.sqrt(2)) * (qb0 - 1j * qb1)


def normalizer(N):
	return 1.0 / math.sqrt(N)

def state2dm(psi):
	return psi * psi.dag()

def qprint(msg, q):
	print("{}:\n{}\n".format(msg, q))


def test01():
  print("== TEST 01 ==")

def test02():
	print("== TEST 02 ==")
	
def test03():
	print("== TEST 03 ==")


def calc_Hmin(rhoX):
	qprint('rhoX', rhoX)
	evalues = rhoX.eigenenergies()
	print("evalues={}".format(evalues))
	lmax = max(evalues)
	Hmin = (-1.0) * math.log(lmax,2)
	print("lmax = %g   ==>  Hmin = %g" % (lmax, Hmin))
	return Hmin


def calc_pGuess(psiA, plusKet, minusKet):
	rhoA = ket2dm(psiA)
	rhoPlus = ket2dm(plusKet)
	rhoMinus = ket2dm(minusKet)
	pPlus = (rhoA * rhoPlus).tr().real
	pMinus = (rhoA * rhoMinus).tr().real
	pGuess = max(pPlus, pMinus)
	print("p+=%g p-=%g pGuess=%g" % (pPlus, pMinus, pGuess))
	return pGuess

def get_H_Golay():
	H = np.matrix( '1 0 0 1 1 1 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0;	\
					1 0 1 0 1 1 0 1 1 0 0 1 0 1 0 0 0 0 0 0 0 0 0;	\
					1 0 1 1 0 1 1 0 1 0 1 0 0 0 1 0 0 0 0 0 0 0 0;	\
					1 0 1 1 1 0 1 1 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0;	\
					1 1 0 0 1 1 1 0 1 1 0 0 0 0 0 0 1 0 0 0 0 0 0;	\
					1 1 0 1 0 1 1 1 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0;	\
					1 1 0 1 1 0 0 1 1 0 1 0 0 0 0 0 0 0 1 0 0 0 0;	\
					1 1 1 0 0 1 0 1 0 1 1 0 0 0 0 0 0 0 0 1 0 0 0;	\
					1 1 1 0 1 0 1 0 0 0 1 1 0 0 0 0 0 0 0 0 1 0 0;	\
					1 1 1 1 0 0 0 0 1 1 0 1 0 0 0 0 0 0 0 0 0 1 0;	\
					0 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 1')
	return H

def get_H_Hamming():
	H = np.matrix(	'1 0 1 0 1 0 1; 	\
					 0 1 1 0 0 1 1; 	\
					 0 0 0 1 1 1 1');
	return H
	
	

def exercise_1a():
	print("== EXERCISE 1A ==")
	H_Hamming = get_H_Hamming()
	H_Golay = get_H_Golay()
	print("H-rank={} shape={}".format(np.linalg.matrix_rank(H_Hamming), H_Hamming.shape))
	print("G-rank={} shape={}".format(np.linalg.matrix_rank(H_Golay), H_Golay.shape))

	


def exercise_1d():
	print("== EXERCISE 1D ==")
	H_Hamming = get_H_Hamming()
	H_Golay = get_H_Golay()
	qprint('HH', H_Hamming)
	qprint('HG', H_Golay)

	v0 = np.matrix('1 0 1 1 0 1 1').T
	v1 = np.matrix('1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0').T
	qprint('v0T', v0.T)
	qprint('v1T', v1.T)

	v0E = H_Hamming * v0
	v1E = H_Golay * v1
	qprint('v0ET', v0E.T)
	qprint('v1ET', v1E.T)
	

def exercise_2():
	print("== EXERCISE 1D ==")
	H_Hamming = get_H_Hamming()
	H_Golay = get_H_Golay()
	v = np.matrix('1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0').T
	s0 = mp.matrix('1 1 1 0 0 0 0 1 1 1 0').T
	s1 = mp.matrix('1 0 1 0 0 1 1 1 0 1 0').T
	


	


def main():
	#test01()
	#test02()
	#test03()
	exercise_1a()
	exercise_1d()

if __name__ == "__main__":
	main()


