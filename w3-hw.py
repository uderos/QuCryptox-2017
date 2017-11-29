#!/home/ubeuser/anaconda3/bin/python

#####################################################################
# WEEK3 - HOMEWORK
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

	

def problem_3_2():
	print("== PAGE 3 PROBLEM 2 ==")
	alpha = 1.0 / math.sqrt(2 + math.sqrt(2))
	psiA = alpha * (qb0 + qbHp)
	psiPlus1 = qb0
	psiMinus1 = qb1
	psiPlus2 = qbHp
	psiMinus2 = qbHm
	psiPlus3 = qbY0
	psiMinus3 = qbY1
	pGuess1 = calc_pGuess(psiA, psiPlus1, psiMinus1)
	pGuess2 = calc_pGuess(psiA, psiPlus2, psiMinus2)
	pGuess3 = calc_pGuess(psiA, psiPlus3, psiMinus3)
	pGuess = (pGuess1 + pGuess2 + pGuess3) / 3.0
	print("pGuess = %r" % pGuess)


	


def problem_1_2():
	print("== PAGE 3 PROBLEM 2 ==")

	rhoX = ket2dm(qb00)
	calc_Hmin(rhoX)

	rhoX = (0.5 * ket2dm(qb00)) + (0.5 * ket2dm(qb11))
	calc_Hmin(rhoX)

	rhoX = ((3.0/4.0) * ket2dm(qb0)) + ((1.0/4.0) * ket2dm(qb1))
	calc_Hmin(rhoX)

	rhoX = ((3.0/4.0) * ket2dm(qbHp)) + ((1.0/4.0) * ket2dm(qbHm))
	calc_Hmin(rhoX)

	epsilon = 0.01
	rho00 = 0.25 * ket2dm(qb00)
	rho11 = 0.25 * ket2dm(qb11)
	rho01 = (0.25 - epsilon) * ket2dm(qb01)
	rho10 = (0.25 + epsilon) * ket2dm(qb10)
	rhoX = rho00 + rho01 + rho10 + rho11
	calc_Hmin(rhoX)





def main():
	#test01()
	#test02()
	#test03()
	problem_1_2()
	problem_3_2()

if __name__ == "__main__":
	main()


