#!/home/ubeuser/anaconda3/bin/python

#####################################################################
# WEEK 6 - LAB
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

qbY0 = (1.0 / math.sqrt(2)) * (qb0 + (1j * qb1))
qbY1 = (1.0 / math.sqrt(2)) * (qb0 - (1j * qb1))


def normalizer(N):
	return 1.0 / math.sqrt(N)

def state2dm(psi):
	return psi * psi.dag()

def qprint(msg, q):
	print("{}:\n{}\n".format(msg, q))


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
	puess1 = max(pPlus, pMinus)
	print("p+=%g p-=%g puess1=%g" % (pPlus, pMinus, puess1))
	return puess1


def problem_1a_6S():
	print("== PROBLEM 1A 6S ==")
	pTheta0 = 0.5
	pTheta1 = 0.5
	rho0 = (1.0/3.0) * ( ket2dm(qb0) + ket2dm(qbHp) + ket2dm(qbY0) )
	qprint('rho0', rho0)
	rho1 = (1.0/3.0) * ( ket2dm(qb1) + ket2dm(qbHm) + ket2dm(qbY1) )
	qprint('rho1', rho1)
	td = tracedist(rho0, rho1)
	pguess2 = 0.5 + (0.5 * td)
	print("td={} pguess2={}".format(td, pguess2))
	
def problem_1a_BB84():
	print("== PROBLEM 1A BB84 ==")
	pTheta0 = 0.5
	pTheta1 = 0.5
	rho0 = (0.5 * ket2dm(qb0)) + (0.5 * ket2dm(qbHp))
	qprint('rho0', rho0)
	rho1 = (0.5 * ket2dm(qb1)) + (0.5 * ket2dm(qbHm))
	qprint('rho1', rho1)
	td = tracedist(rho0, rho1)
	pguess2 = 0.5 + (0.5 * td)
	print("td={} pguess2={}".format(td, pguess2))
	

def problem_1b_BB84():
	print("== PROBLEM 1B BB84 ==")
	pTheta0 = 0.5
	pTheta1 = 0.5
	rho0 = 0.5 * (ket2dm(qb0) + ket2dm(qbHp))
	rho1 = 0.5 * (ket2dm(qb1) + ket2dm(qbHm))
	rho01 = rho0 - rho1
	e = rho01.eigenstates()
	qprint('e', e)
	eValues = e[0]
	eVects = e[1]
	m0 = qeye(2)
	m1 = qeye(2)
	for i in range(len(eValues)):
		if eValues[i] >= 0:
			m0 = m0 + ket2dm(eVects[i])
		else:
			m1 = m1 + ket2dm(eVects[i])
	m0 = m0 / m0.tr()
	m1 = m1 / m1.tr()
	qprint('m0', m0)
	qprint('m1', m1)
	qprint('m0+m1', m0+m1)
	rho = (pTheta0 * rho0) + (pTheta1 * rho1)
	rhoPost0 = (m0 * rho * m0.dag()) / (m0.dag()*m0*rho).tr()
	qprint('rhoPost0', rhoPost0)
	rhoPost1 = (m1 * rho * m1.dag()) / (m1.dag()*m1*rho).tr()
	qprint('rhoPost1', rhoPost1)
	td0 = tracedist(rho, rhoPost0)
	td1 = tracedist(rho, rhoPost1)
	print("td0={} td1={}".format(td0, td1))
	

def problem_1b_6S():
	print("\n\n== PROBLEM 1B 6S ==")
	pTheta0 = 0.5
	pTheta1 = 0.5
	rho0 = (1.0/3.0) * ( ket2dm(qb0) + ket2dm(qbHp) + ket2dm(qbY0) )
	qprint('rho0', rho0)
	rho1 = (1.0/3.0) * ( ket2dm(qb1) + ket2dm(qbHm) + ket2dm(qbY1) )
	rho01 = rho0 - rho1
	e = rho01.eigenstates()
	qprint('e', e)
	eValues = e[0]
	eVects = e[1]
	m0 = qeye(2)
	m1 = qeye(2)
	for i in range(len(eValues)):
		if eValues[i] >= 0:
			m0 = m0 + ket2dm(eVects[i])
		else:
			m1 = m1 + ket2dm(eVects[i])
	m0 = m0 / m0.tr()
	m1 = m1 / m1.tr()
	qprint('m0', m0)
	qprint('m1', m1)
	qprint('m0+m1', m0+m1)
	rho = (pTheta0 * rho0) + (pTheta1 * rho1)
	rhoPost0 = (m0 * rho * m0.dag()) / (m0.dag()*m0*rho).tr()
	qprint('rhoPost0', rhoPost0)
	rhoPost1 = (m1 * rho * m1.dag()) / (m1.dag()*m1*rho).tr()
	qprint('rhoPost1', rhoPost1)
	td0 = tracedist(rho, rhoPost0)
	td1 = tracedist(rho, rhoPost1)
	print("td0={} td1={}".format(td0, td1))
	

def test_dephasing_noise(p, rho, msg):
	X = sigmax()
	rho_with_noise = ((1 - p)*rho) + (p * X * rho * X)
	td = tracedist(rho, rho_with_noise)
	#qprint('rho', rho)
	#qprint('rhoWN', rho_with_noise)
	print("DePhasing %s ==> td = %.3f" % (msg, td))

def test_depolarizing_noise(p, rho, msg):
	X = sigmax()
	rho_with_noise = ((1 - p)*rho) + (p * 0.5 * qeye(2))
	td = tracedist(rho, rho_with_noise)
	#qprint('rho', rho)
	#qprint('rhoWN', rho_with_noise)
	print("Depolarizing: %s ==> td = %.3f" % (msg, td))

def problem_3():
	test_dephasing_noise(0.2,  ket2dm(qb0), 	"1")
	test_dephasing_noise(0.3,  ket2dm(qb0), 	"2")
	test_dephasing_noise(0.2,  ket2dm(qb1), 	"3")
	test_dephasing_noise(0.3,  ket2dm(qb1), 	"4")
	test_dephasing_noise(0.2,  ket2dm(qbHp), 	"5")
	test_dephasing_noise(0.3,  ket2dm(qbHp), 	"6")
	test_dephasing_noise(0.2,  ket2dm(qbHm), 	"7")
	test_dephasing_noise(0.3,  ket2dm(qbHm), 	"8")
	
	test_depolarizing_noise(0.2,  ket2dm(qb0), 	 	 "9")
	test_depolarizing_noise(0.3,  ket2dm(qb0), 		"10")
	test_depolarizing_noise(0.2,  ket2dm(qb1), 		"11")
	test_depolarizing_noise(0.3,  ket2dm(qb1), 		"12")
	test_depolarizing_noise(0.2,  ket2dm(qbHp), 	"13")
	test_depolarizing_noise(0.3,  ket2dm(qbHp), 	"14")
	test_depolarizing_noise(0.2,  ket2dm(qbHm), 	"15")
	test_depolarizing_noise(0.3,  ket2dm(qbHm), 	"16")

def calc_p_flip_dephasing(rho, p):
	X = sigmax()
	rhoWN = ((1 - p)*rho) + (p * X * rho * X)
	td = tracedist(rho, rhoWN)
	return td

def calc_p_flip_depolarizing(rho, p):
	X = sigmax()
	rhoWN = ((1 - p)*rho) + (p * 0.5 * qeye(2))
	td = tracedist(rho, rhoWN)
	return td

def show_dephasing(rho, p):
	X = sigmax()
	rhoWN = ((1 - p)*rho) + (p * X * rho * X)
	noise = rho - rhoWN
	print("p={}".format(p))
	qprint('rho', rho)
	qprint('noise', noise)
	print("td={}".format(tracedist(rho, rhoWN)))

def problem_4_1():
	fname = "== PROBLEM 4.1 =="
	states = [ ket2dm(qb0), ket2dm(qb1), ket2dm(qbHp), ket2dm(qbHm) ]
	Pstates = [ 0.25, 0.25, 0.25, 0.25 ]
	p = 0.1
	pflip = 0
	for i in range(len(states)):
		pflip += Pstates[i] * calc_p_flip_dephasing(states[i], p)
	print("%s: p = %g   pflip = %.3f" % (fname, p, pflip))
	
	
def problem_4_2():
	fname = "== PROBLEM 4.2 =="
	states = [ ket2dm(qb0), ket2dm(qb1), ket2dm(qbHp), ket2dm(qbHm) ]
	Pstates = [ 0.25, 0.25, 0.25, 0.25 ]
	p = 0.1
	pflip = 0
	for i in range(len(states)):
		pflip += Pstates[i] * calc_p_flip_depolarizing(states[i], p)
	print("%s: p = %g   pflip = %.3f" % (fname, p, pflip))
	
	
def problem_4_3():
	fname = "== PROBLEM 4.3 =="
	states = [ ket2dm(qb0), ket2dm(qb1), ket2dm(qbHp), ket2dm(qbHm), ket2dm(qbY0), ket2dm(qbY1) ]
	Pstates = [ (1.0/6.0), (1.0/6.0), (1.0/6.0), (1.0/6.0), (1.0/6.0), (1.0/6.0) ]
	p = 0.1
	pflip = 0
	for i in range(len(states)):
		pflip += Pstates[i] * calc_p_flip_dephasing(states[i], p)
	print("%s: p = %g   pflip = %.3f" % (fname, p, pflip))
	
	
def problem_4_4():
	fname = "== PROBLEM 4.3 =="
	states = [ ket2dm(qb0), ket2dm(qb1), ket2dm(qbHp), ket2dm(qbHm), ket2dm(qbY0), ket2dm(qbY1) ]
	Pstates = [ (1.0/6.0), (1.0/6.0), (1.0/6.0), (1.0/6.0), (1.0/6.0), (1.0/6.0) ]
	p = 0.1
	pflip = 0
	for i in range(len(states)):
		pflip += Pstates[i] * calc_p_flip_depolarizing(states[i], p)
	print("%s: p = %g   pflip = %.3f" % (fname, p, pflip))
	
	


def test01():
	print("== TEST 01 ==")
	show_dephasing(ket2dm(qb0), 0.1)
	show_dephasing(ket2dm(qb1), 0.1)
	show_dephasing(ket2dm(qbHp), 0.1)

def test02():
	print("== TEST 02 ==")
	
def test03():
	print("== TEST 03 ==")


def main():
	#test01()
	#test02()
	#test03()
	#problem_1a_6S()
	#problem_1a_BB84()
	#problem_1b_BB84()
	#problem_1b_6S()
	#problem_3()
	problem_4_1()
	problem_4_2()
	problem_4_3()
	problem_4_4()

if __name__ == "__main__":
	main()


