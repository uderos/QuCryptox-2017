#!/home/ubeuser/anaconda3/bin/python

###########################################################
# WEEK 2 - LAB
###########################################################

import math
import numpy as np
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

def normalizer(N):
	return 1.0 / math.sqrt(N)

#def state2dm(psi):
#	return psi * psi.dag()

def state2dm(psi1, psi2):
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
	
		

def test01():
  print("\n\n== TEST 01 ==")
  ptOp = partial_trace_operator(2)

def test02():
	print("\n\n== TEST 02 ==")
	rho = 0.25 * (state2dm(qb00, qb00) + state2dm(qb01, qb10) + state2dm(qb10, qb01) + state2dm(qb11, qb11))
	is_separable(rho)
	rho = 0.5 * (state2dm(qb00, qb00) + state2dm(qb11, qb11))
	is_separable(rho)
	rho = 0.25 * (state2dm(qb00, qb00) + state2dm(qb00, qb11) + state2dm(qb11, qb00) + state2dm(qb11, qb11))
	is_separable(rho)

def test03():
	print("\n\n== TEST 03 ==")
	print("\n\n")
	rho = state2dm(qb00, qb00)
	is_entangled(rho)

	print("\n\n")
	rho = 0.5 * (state2dm(qb00, qb00) + state2dm(qb00, qb01) + \
				 state2dm(qb01, qb00) + state2dm(qb01, qb01))
	is_entangled(rho)

	print("\n\n")
	rho = 0.5 * (state2dm(qb00, qb00) + state2dm(qb00, qb11) + \
				 state2dm(qb11, qb00) + state2dm(qb11, qb11))
	is_entangled(rho)

	print("\n\n")
	ss = math.sin(math.pi/7.0)
	cc = math.cos(math.pi/7.0)
	rho = (ss*ss*state2dm(qb00, qb00) - ss*cc*state2dm(qb00, qb11) - \
		   ss*cc*state2dm(qb11, qb00) + cc*cc*state2dm(qb11, qb11))
	is_entangled(rho)

#	print("\n\n")
#	f = (-1.0) * math.sqrt(2)/3.0
#	rho = ((1.0/3.0)*state2dm(qb00, qb00) + f*state2dm(qb00, qb11) + \
#		  f*state2dm(qb11, qb00) + (2.0/3.0)*state2dm(qb11, qb11))
#	is_entangled(rho)

def test04():
	print("\n\n== TEST 04 ==")
	psi = normalizer(2) * (qb00 + qb11)
	rho = ket2dm(psi)
	qprint('rho',rho)

	measKet = normalizer(2) * (qb00 + qb01)
	measBra = measKet.dag()
	qprint('measKet', measKet)
	qprint('measBra', measBra)
	p = measBra.overlap(rho * measKet)
	qprint('p', p)
	#print("#### measBra * measKet = %r" % p.data[0])

	qprint('rho * measKet',rho * measKet)
	qprint('measKet.dag() * rho',measKet.dag() * rho)

	expV = measKet.dag() * rho * measKet
	print("expV=%g" % expV)


def page1_question1():
	print("\n\n== PAGE 1 QUESTION 01 ==")
	rhoAB = (0.5*state2dm(qb00, qb00) - 0.3*state2dm(qb00, qb11) -0.3*state2dm(qb11, qb00) + 0.5*state2dm(qb11, qb11))
	qprint('rhoAB', rhoAB)
	is_separable(rhoAB)

def page1_question2():
	print("\n\n== PAGE 1 QUESTION 02 ==")
	rhoAB = (0.5*state2dm(qb00, qb00) + 0.5*state2dm(qb11, qb11))
	qprint('rhoAB', rhoAB)
	is_separable(rhoAB)

def page1_question3():
	print("\n\n== PAGE 1 QUESTION 03 ==")
	rhoAB = 0.25 * (state2dm(qb00, qb00) + state2dm(qb00, qb01) + \
					state2dm(qb01, qb00) + state2dm(qb01, qb01) + \
					state2dm(qb10, qb10) + state2dm(qb10, qb11) + \
					state2dm(qb11, qb10) + state2dm(qb11, qb11))
	qprint('rhoAB', rhoAB)
	is_separable(rhoAB)

def page2_question1():
	print("\n\n== PAGE 2 QUESTION 01 ==")
	f = math.sqrt(2)/3.0
	rho1 = ((1.0/3.0)*state2dm(qb00, qb00) - f*state2dm(qb00, qb11) -f*state2dm(qb11, qb00) + (2.0/3.0)*state2dm(qb11, qb11))
	rho2 = ((1.0/3.0)*state2dm(qb01, qb01) - f*state2dm(qb01, qb10) -f*state2dm(qb10, qb01) + (2.0/3.0)*state2dm(qb10, qb10))
	are_uhlmann_equivalent(rho1, rho2, 0)
	are_uhlmann_equivalent(rho1, rho2, 1)

def page2_question2():
	print("\n\n== PAGE 2 QUESTION 02 ==")
	f = math.sqrt(2)/3.0
	rho1 = ((1.0/3.0)*state2dm(qb00, qb00) - f*state2dm(qb00, qb11) -f*state2dm(qb11, qb00) + (2.0/3.0)*state2dm(qb11, qb11))
	rho2 = ((1.0/3.0)*state2dm(qb01, qb01) + (2.0/3.0)*state2dm(qb10, qb10))
	#qprint('rho1', rho1)
	#qprint('rho2', rho2)
	are_uhlmann_equivalent(rho1, rho2, 0)
	are_uhlmann_equivalent(rho1, rho2, 1)

def show_shmidt_coeffs(psi):
	rho = state2dm(psi, psi)
	coeffs = shmidt_coeffs = calc_shmidt_coeffs(rho)
	#qprint('psi', psi)
	#qprint('rho', rho)
	print("coeffs={}".format(coeffs))

def page3_question1():
	print("\n\n== PAGE 3 QUESTION 01 ==")
	show_shmidt_coeffs(qb00)
	show_shmidt_coeffs(normalizer(2)*(qb00 + qb01))
	show_shmidt_coeffs(normalizer(2)*(qb00 + qb11))

	ss = math.sin(math.pi/7.0)
	cc = math.cos(math.pi/7.0)
	show_shmidt_coeffs(ss*qb00 + cc*qb11)

#def meas_outcome_d2(rhoAB, meas_A_ket, meas_B_ket):
#	op = tensor(meas_A_ket, meas_B_ket)
#	meas_outcome = op.dag().overlap(rhoAB * op)
#	return meas_outcome 

def meas_outcome_d2(rhoAB, meas_A_ket, meas_B_ket):
	meas_ket = tensor(meas_A_ket, meas_B_ket)
	meas_dm = ket2dm(meas_ket)
	outcome = (rhoAB * meas_dm).tr()
	return outcome

def calc_chsh_pwin(rhoAB, thetaA0, thetaA1, thetaB0, thetaB1):
	measA0_x0 = math.cos(thetaA0)*qb0 + math.sin(thetaA0)*qb1
	measA1_x0 = (-1.0)*math.sin(thetaA0)*qb0 + math.cos(thetaA0)*qb1
	measA0_x1 = math.cos(thetaA1)*qb0 + math.sin(thetaA1)*qb1
	measA1_x1 = (-1.0)*math.sin(thetaA1)*qb0 + math.cos(thetaA1)*qb1
	measB0_y0 = math.cos(thetaB0)*qb0 + math.sin(thetaB0)*qb1
	measB1_y0 = (-1.0)*math.sin(thetaB0)*qb0 + math.cos(thetaB0)*qb1
	measB0_y1 = math.cos(thetaB1)*qb0 + math.sin(thetaB1)*qb1
	measB1_y1 = (-1.0)*math.sin(thetaB1)*qb0 + math.cos(thetaB1)*qb1

	qprint('measA0_x0',measA0_x0)
	qprint('measA1_x0',measA1_x0)
	qprint('measA0_x1',measA0_x1)
	qprint('measA1_x1',measA1_x1)
	qprint('measB0_y0',measB0_y0)
	qprint('measB1_y0',measB1_y0)
	qprint('measB0_y1',measB0_y1)
	qprint('measB1_y1',measB1_y1)

	#ab xy
	p00_00 = meas_outcome_d2(rhoAB, measA0_x0, measB0_y0)
	p11_00 = meas_outcome_d2(rhoAB, measA1_x0, measB1_y0)
	p00_10 = meas_outcome_d2(rhoAB, measA0_x1, measB0_y0)
	p11_10 = meas_outcome_d2(rhoAB, measA1_x1, measB1_y0)
	p00_01 = meas_outcome_d2(rhoAB, measA0_x0, measB0_y1)
	p11_01 = meas_outcome_d2(rhoAB, measA1_x0, measB1_y1)
	p10_11 = meas_outcome_d2(rhoAB, measA1_x1, measB0_y1)
	p01_11 = meas_outcome_d2(rhoAB, measA0_x1, measB1_y1)

	print("p00_00=%r p11_00=%r" % (p00_00, p11_00))
	print("p00_10=%r p11_10=%r" % (p00_10, p11_10))
	print("p00_01=%r p11_01=%r" % (p00_01, p11_01))
	print("p10_11=%r p01_11=%r" % (p10_11, p10_11))

	pwin_00 = p00_00 + p11_00
	pwin_01 = p00_01 + p11_01
	pwin_10 = p00_10 + p11_10
	pwin_11 = p10_11 + p10_11

	pwin = 0.25*(pwin_00 + pwin_01 + pwin_10 + pwin_11).real
	print("Pwin = %r" % pwin)

	return pwin

def page4_question1_A():
	print("\n\n== PAGE 4 QUESTION 01-A ==")
	thetaA0 = 0
	thetaA1 = math.pi / 4.0
	thetaB0 = math.pi / 8.0
	thetaB1 = math.pi / (-8.0)

	psi_good = normalizer(2) * (qb00 + qb11)
	rho_good = ket2dm(psi_good)
	rho_bad  = ket2dm(qb01)
	rhoAB    = 0.75 * rho_good + 0.25 * rho_bad

	qprint('psi_good', psi_good)
	qprint('rho_good', rho_good)
	qprint('rho_bad ', rho_bad )
	qprint('rhoAB   ', rhoAB   )

	pWin = calc_chsh_pwin(rhoAB, thetaA0, thetaA1, thetaB0, thetaB1)
	print("\n\t ## PWIN = %r" % pWin)

def page4_question1_B():
	print("\n\n== PAGE 4 QUESTION 01-B ==")
	thetaA0 = 0
	thetaA1 = math.pi / 4.0
	thetaB0 = math.pi / 8.0
	thetaB1 = math.pi / (-8.0)

	psi_good = normalizer(2) * (qb00 + qb11)
	rho_good = ket2dm(psi_good)
	rho_bad  = ket2dm(qb01)
	rhoAB    = 0.25 * rho_good + 0.75 * rho_bad

	qprint('psi_good', psi_good)
	qprint('rho_good', rho_good)
	qprint('rho_bad ', rho_bad )
	qprint('rhoAB   ', rhoAB   )

	pWin = calc_chsh_pwin(rhoAB, thetaA0, thetaA1, thetaB0, thetaB1)
	print("\n\t ## PWIN = %r" % pWin)

def page4_question1_C():
	print("\n\n== PAGE 4 QUESTION 01-B ==")
	thetaA0 = 0
	thetaA1 = math.pi / 4.0
	thetaB0 = 0
	thetaB1 = math.pi / 4.0

	psi_good = normalizer(2) * (qb00 + qb11)
	rhoAB    = ket2dm(psi_good)
	qprint('rhoAB   ', rhoAB   )

	pWin = calc_chsh_pwin(rhoAB, thetaA0, thetaA1, thetaB0, thetaB1)
	print("\n\t ## PWIN = %r" % pWin)

def main():
	#test01()
	#test02()
	#test03()
	#test04()
	#page1_question1()
	#page1_question2()
	#page1_question3()
	#page2_question1()
	#page2_question2()
	#page3_question1()
	page4_question1_A()
	page4_question1_B()
	page4_question1_C()

if __name__ == "__main__":
	main()


