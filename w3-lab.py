#!/home/ubeuser/anaconda3/bin/python

###########################################################
# WEEK 3 - LAB
###########################################################

import math
import numpy as np
import scipy
#import matplotlib as plt
from qutip import *

#################################################
# CONSTANTS
#################################################

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


#################################################
# UTILITY FUNCTIONS
#################################################

def normalizer(N):
	return 1.0 / math.sqrt(N)

def kets2dm(psi1, psi2):
	return psi1 * psi2.dag()

def qprint(msg, q):
	print("{}:\n{}\n".format(msg, q))

def makeWstate(N):
	q = []
	for i in range(N):
		l = []
		for j in range(N):
			if (i == j):
				l.append(qb1)
			else:
				l.append(qb0)
		sub_state = tensor(l)
		q.append(sub_state)
	wN = (1.0/math.sqrt(N)) * sum(q)
	return wN


def partial_trace_operator(N):
	ptOp = qzero(N)
	if N == 2:
		Tb = 0.25 * (state2dm(qb00, qb00) + state2dm(qb01, qb10) + state2dm(qb10, qb01) + state2dm(qb11, qb11))
		qprint('Tb', Tb)
		ptOp = tensor(qeye(2), Tb)
		qprint('ptOp', ptOp)
	return ptOp

def is_separable(rho):
	rho_pt = partial_transpose(rho, [0,1])
	evalues = rho_pt.eigenenergies()
	qprint('rho', rho)
	qprint('rho_pt', rho_pt)
	print("evalues={}".format(evalues))
	is_separable_flag = True
	for e in evalues:
		if e < 0:
			is_separable_flag = False
	if is_separable_flag:
		print("SEPARABLE\n")
	else:
		print("NOT separable\n")
	return is_separable_flag


def calc_shmidt_coeffs(rho):
	rhoA = rho.ptrace(0)
	#qprint('rhoA', rhoA)
	evalues = rhoA.eigenenergies()
	coeffs = [ math.sqrt(i) for i in evalues ]
	return coeffs

def calc_shmidt_rank(rho):
	epsilon = 1e-6
	shmidt_coeffs = calc_shmidt_coeffs(rho)
	rank = 0
	for coeff in shmidt_coeffs:
		if abs(coeff) > epsilon:
			rank = rank + 1
	return rank

def is_entangled(rho):
	shmidt_rank = calc_shmidt_rank(rho)
	is_entangled_flag = (shmidt_rank > 1)
	qprint('rho', rho)
	print("is_entangled? %r" % is_entangled_flag)
	return is_entangled_flag

def are_uhlmann_equivalent(rho1, rho2, idx=0):
	#print("Are Uhlmann Equivalent: idx=%d" % (idx))
	rho1Pt = rho1.ptrace(idx)
	rho2Pt = rho2.ptrace(idx)
	are_equivalent = (rho1Pt == rho2Pt)
	#qprint('rho1', rho1)
	#qprint('rho2', rho2)
	qprint('rho1Pt', rho1Pt)
	qprint('rho2Pt', rho2Pt)
	print("Are Uhlmann Equivalent: idx=%d %r" % (idx, are_equivalent))
	
def meas_outcome_d2(rhoAB, meas_A_ket, meas_B_ket):
	meas_ket = tensor(meas_A_ket, meas_B_ket)
	meas_dm = ket2dm(meas_ket)
	outcome = (rhoAB * meas_dm).tr()
	return outcome

def udr_minimize(func, x0):
	res = scipy.optimize.minimize(func, x0, 
			method='nelder-mead', 
			options={'xtol': 1e-8, 'disp': False})
	if (not res.status == 0):
		raise RuntimeError("minimize() failure: {}".format(res.status))
	minX = res.x[0]
	minY = res.fun
	return [minX, minY]


def udr_calc_Hmin(rhoAE):
	#qprint('rhoAE', rhoAE)
	Hmin0 = entropy_conditional(rhoAE, 2)
	Hmin1 = entropy_conditional(rhoAE, 2)
	print("Hmin1=%g  (Hmin0=%g)" % (Hmin1, Hmin0))
	return Hmin1

		
#################################################
# TESTS
#################################################

def testfun1(x):
	return 1.0 + ((x - 2.0) ** 2)

def testfun2(x):
	return 1.0 + 1.0 / (1.0 + (x - 2.0) ** 2)

def test01():
	print("\n\n== TEST 01 ==")
	x0 = 0
	res = scipy.optimize.minimize(testfun1, x0, 
			method='nelder-mead', 
			options={'xtol': 1e-8, 'disp': True})
	print("X0=%r res=%r" % (x0, res))
	print("minY=%g minX=%g" % (res.fun, res.x[0]))
	

def test02():
	print("\n\n== TEST 02 ==")
	x0 = 0
	[minX, minY] = udr_minimize(testfun1, x0)
	print("minY=%g minX=%g" % (minY, minX))

def test03():
	print("\n\n== TEST 03 ==")
	A = 0.5 * qeye(2)
	E = ket2dm(qb0)
	rhoAE = tensor(A, E)
	qprint('rhoAE', rhoAE)
	Hmin0 = entropy_conditional(rhoAE, 0, 2)
	Hmin1 = entropy_conditional(rhoAE, 1, 2)
	print("Hmin0=%g  Hmin1=%g" % (Hmin0, Hmin1))



#################################################
# LAB PROBLEMS
#################################################

def page1_question1():
	print("\n\n== PAGE 1 QUESTION 01 ==")
	rho0 = 0.5 * (ket2dm(qb0) + ket2dm(qbHp))
	rho1 = 0.5 * (ket2dm(qb1) + ket2dm(qbHm))
	qprint('rho0', rho0)
	qprint('rho1', rho1)

	# Standard basis ==> 0.75
	M0 = ket2dm(qb0)
	M1 = ket2dm(qb1)
	pguess = 0.5 * (M0*rho0).tr() + 0.5 * (M1*rho1).tr()
	qprint('M0', M0)
	qprint('M1', M1)
	print("Pguess = %r" % pguess)

	# My guess ==> 0.75
	M0 = rho0
	M1 = rho1
	pguess = 0.5 * (M0*rho0).tr() + 0.5 * (M1*rho1).tr()
	qprint('M0', M0)
	qprint('M1', M1)
	print("Pguess = %r" % pguess)

	# My guess ==> 0.25
	M0 = rho1
	M1 = rho0
	pguess = 0.5 * (M0*rho0).tr() + 0.5 * (M1*rho1).tr()
	qprint('M0', M0)
	qprint('M1', M1)
	print("Pguess = %r" % pguess)

	# My guess ==> 0.125
	M0 = (0.5 * rho0) + (0.5 * rho1)
	M1 = (0.5 * rho0) - (0.5 * rho1)
	pguess = 0.5 * (M0*rho0).tr() + 0.5 * (M1*rho1).tr()
	qprint('M0', M0)
	qprint('M1', M1)
	print("Pguess = %r" % pguess)

	# Hadamard basis ==> 0.75
	M0 = ket2dm(qbHp)
	M1 = ket2dm(qbHm)
	pguess = 0.5 * (M0*rho0).tr() + 0.5 * (M1*rho1).tr()
	qprint('M0', M0)
	qprint('M1', M1)
	print("Pguess = %r" % pguess)

	# Try this:
	N = 25
	for i in range(N):
		alpha = math.pi * (i / N)
		cc = math.cos(alpha)
		ss = math.sin(alpha)
		M0 = ket2dm(cc * qb0 + ss * qb1)
		M1 = ket2dm((-1.0) * ss * qb0 + cc * qb1)
		pguess = 0.5 * (M0*rho0).tr() + 0.5 * (M1*rho1).tr()
		print("alpha=%g Pguess = %r" % (alpha, pguess))

	# Search for the best
	print('\n Srearching ...')
	step = 0.01
	alphaL = 0
	alphaR =  alphaL + step
	pguessMax = 0.0
	while (alphaR - alphaL >= 1e-6):
		alpha = alphaR
		cc = math.cos(alpha)
		ss = math.sin(alpha)
		M0 = ket2dm(cc * qb0 + ss * qb1)
		M1 = ket2dm((-1.0) * ss * qb0 + cc * qb1)
		pguess = (0.5 * (M0*rho0).tr() + 0.5 * (M1*rho1).tr()).real
		print("alpha=%g = %gpi step=%e Pguess=%r" % (alpha, alpha, step, pguess))
		if (pguess > pguessMax):
			pguessMax = pguess
			alphaL = alphaR
			alphaR = alphaR + step
		else:
			step = step * 0.1
			alphaR = alphaL + step

			
def udr_tracedist1(rho0, rho1):
	A = rho0 - rho1
	#qprint('A', A)
	evalues = A.eigenenergies()
	#print("evals = {}".format(evalues))
	td = 0
	for ev in evalues:
		td = td + abs(ev)
	return 0.5*td

def udr_tracedist2(rho0, rho1):
	A = rho0 - rho1
	#qprint('A', A)
	evalues = A.eigenenergies()
	#print("evals = {}".format(evalues))
	td = 0
	for ev in evalues:
		if (ev > 0):
			td = td + ev
	return td

def show_trace_distance(rho0, rho1):
	#qprint('rho0', rho0)
	#qprint('rho1', rho1)
	td = tracedist(rho0, rho1)
	td_udr1 = udr_tracedist1(rho0, rho1)
	td_udr2 = udr_tracedist2(rho0, rho1)
	err1 = td_udr1 - td
	err2 = td_udr2 - td
	print("Trace Distance = %r err1=%e err2=%e" % (td, err1, err2))

def page2_question1():
	print("\n\n== PAGE 2 QUESTION 01 ==")
	show_trace_distance(ket2dm(qb00), ket2dm(qb00))
	show_trace_distance(ket2dm(qb00), ket2dm(qb11))

	phi0 = ket2dm(qb00)
	phi1 = 0.5 * (ket2dm(qb00 + qb11))
	show_trace_distance(phi0, phi1)

	show_trace_distance(ket2dm(qb0), ket2dm(qbHp))


def calc_Hmin(rhoE0, rhoE1):
	td = tracedist(rhoE0, rhoE1)
	Pguess = 0.5 + 0.5 * td
	Hmin = (-1.0) * math.log(Pguess, 2)
	print("td=%g pG=%g Hmin=%g" % (td, Pguess, Hmin))
	return Hmin

#def calcM0(alpha):
#	cc = math.cos(alpha)
#	ss = math.sin(alpha)
#	M0 = ket2dm(cc * qb0 + ss * qb1)
#	return M0
#	
#def calcPguessStep(rhoXE, p0, p1, M0)
#	M1 = eye(2) - M0
#	res = p0 * (M0 * rhoXE).tr() + p1 * (M1 * rhoXE).tr()
#	return res

	
def page3_question1():
	print("\n\n== PAGE 3 QUESTION 01 ==")

	rhoA = 0.5 * qeye(2)
	rhoE = ket2dm(qb0)
	#qprint('rhoA', rhoA)
	#qprint('rhoE', rhoE)
	rhoAE = tensor(rhoA, rhoE)
	Hmin = udr_calc_Hmin(rhoAE)
	print("Question1: Hmin=%g" % Hmin)

	psiAE = normalizer(2) * (qb00 + qb11)
	rhoAE = ket2dm(psiAE)
	Hmin = udr_calc_Hmin(rhoAE)
	print("Question2: Hmin=%g" % Hmin)

	rhoAE1 = 0.5 * tensor(ket2dm(qb0), ket2dm(qb0))
	rhoAE2 = 0.5 * tensor(ket2dm(qb1), ket2dm(qbHp))
	rhoAE = rhoAE1 + rhoAE2
	qprint('rhoAE', rhoAE)
	Hmin = udr_calc_Hmin(rhoAE)
	print("Question3: Hmin=%g" % Hmin)

def page3_question2():
	print("\n\n== PAGE 3 QUESTION 02 ==")
	p00 = 1.0 / 4.0
	p01 = 1.0 / 4.0
	p10 = 1.0 / 10.0
	p11 = 4.0 / 10.0

	w00 = Qobj([[4],[1],[9],[7]])
	qprint('w00', w00)
	w00 = w00.unit()
	rho00 = w00 * w00.dag()
	qprint('rho00', rho00)

	w01 = Qobj([[8],[1],[2],[0]])
	qprint('w01', w01)
	w01 = w01.unit()
	rho01 = w01 * w01.dag()
	qprint('rho01', rho01)

	w10 = Qobj([[5],[9],[1],[8]])
	qprint('w10', w10)
	w10 = w10.unit()
	rho10 = w10 * w10.dag()
	qprint('rho10', rho10)

	w11 = Qobj([[0],[1],[2],[1]])
	qprint('w11', w11)
	w11 = w11.unit()
	rho11= w11 * w11.dag()
	qprint('rho11', rho11)

	rhoXE = p00 * tensor(ket2dm(qb00), rho00) + \
			p01 * tensor(ket2dm(qb01), rho01) + \
			p10 * tensor(ket2dm(qb10), rho10) + \
			p11 * tensor(ket2dm(qb11), rho11)
	qprint('rhoXE', rhoXE)

	Hmin = udr_calc_Hmin(rhoXE)
	print("Hmin=%g" % Hmin)



def main():
	#test01()
	#test02()
	#test03()
	#test04()
	#page1_question1()
	#page2_question1()
	page3_question1()
	#page3_question2()

if __name__ == "__main__":
	main()


