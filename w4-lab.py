#!/home/ubeuser/anaconda3/bin/python

#####################################################################
# WEEK4 - LAB
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

	

def exercise_2a():
	print("== EXERCISE 2A ==")
	rhoExt = ket2dm(qb0)
	rhoU = (1.0/8.0)*qeye(2)
	qprint('rhoExt', rhoExt)
	qprint('rhoU', rhoU)
	td = tracedist(rhoExt, rhoU)
	print("D(rhoExt,rhoU)={}".format(td))
	

def exercise_2b():
	print("== EXERCISE 2B ==")
	rhoExt = 0.25*ket2dm(qb0) + 0.75*ket2dm(qb1)
	rhoU = 0.5*qeye(2)
	qprint('rhoExt', rhoExt)
	qprint('rhoU', rhoU)
	td = tracedist(rhoExt, rhoU)
	print("tr(rhoExt)={} trace(rhoU)={}".format(rhoExt.tr(), rhoU.tr()))
	print("D(rhoExt,rhoU)={}".format(td))
	

	


def main():
	#test01()
	#test02()
	#test03()
	exercise_2a()
	exercise_2b()

if __name__ == "__main__":
	main()


