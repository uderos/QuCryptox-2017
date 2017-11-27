#!/home/ubeuser/anaconda3/bin/python

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

def state2dm(psi):
	return psi * psi.dag()

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
		

def test01():
  print("== TEST 01 ==")
  psi1 = (1/math.sqrt(2)) * (qb00 + qb10)
  dm1 = state2dm(psi1);
  print("state={}\ndm={}".format(psi1, dm1))

def test02():
	print("== TEST 02 ==")
	psi1 = (1/math.sqrt(2))*(qb00 + qb10)
	psi2 = (1/math.sqrt(2))*(qb01 + qb10)
	psi = tensor(psi1, psi2)
	rho = state2dm(psi)
	qprint('psi1', psi1)
	qprint('psi2', psi2)
	qprint('psi', psi)
	qprint('rho', rho)
	tr0 = rho.ptrace(0)
	tr1 = rho.ptrace(1)
	tr2 = rho.ptrace(2)
	qprint('tr0', tr0)
	qprint('tr1', tr1)
	qprint('tr2', tr2)
	tr = rho.ptrace([0, 1, 2])
	qprint('tr', tr)
	
def test03():
	print("== TEST 03 ==")
	w2 = makeWstate(2)
	qprint('w2', w2)
	w3 = makeWstate(3)
	qprint('w3', w3)


def problem_2_1():
	N=4
	ghzN = state2dm(ghz_state(N))
	ghzTrN = ghzN.ptrace(list(range(0,N-1)))
	qprint('ghzN', ghzN)
	qprint('ghzTrN', ghzTrN)
	#rank_ghzN = np.linalg.matrix_rank(ghzN)
	#rank_ghzTrN = np.linalg.matrix_rank(ghzTrN)
	#print("rank_ghzN=%g" % rank_ghzN)
	#print("rank_ghzTrN=%g" % rank_ghzTrN)

	N=4
	wN = state2dm(makeWstate(N))
	wTrN = wN.ptrace(list(range(0,N-1)))
	qprint('wN', N*wN)
	qprint('wTrN', N*wTrN)

	#rankW = np.linalg.matrix_rank(wN)
	#rankTrW = np.linalg.matrix_rank(wTrN)
	#print("rankW=%g" % rankW)
	#print("rankTrW=%g" % rankTrW)

def problem_2_2():
	rho = state2dm(qb0)
	qprint("rho", rho)
	rho2 = rho*rho
	qprint("rho2", rho2)
	tr = rho2.tr()
	qprint("trace(rho2)", tr)
	
def problem_2_4():
	N=4
	psi = ghz_state(N)
	print("ghz{}={}", N, psi)
	rho = state2dm(psi)
	trN = rho.ptrace([1,2,3])
	trN2 = trN * trN
	qprint('rho', rho)
	qprint('trN', trN)
	qprint('trN2', trN2)
	purity = trN2.tr()
	print("purity=%g" % (purity))

def problem_2_3():
	psi = normalizer(4) * (qb00 + qb01 + qb10 + qb11)
	rho = state2dm(psi)
	trA = rho.ptrace(1)
	trA2 = trA * trA
	purity = trA2.tr()
	qprint('psi', psi)
	print("purity=%g" % purity)

	psi = normalizer(3) * (qb00 + qb10 + qb11)
	rho = state2dm(psi)
	trA = rho.ptrace(1)
	trA2 = trA * trA
	purity = trA2.tr()
	qprint('psi', psi)
	print("purity=%g" % purity)

	psi = build_state([qb00, qb11])
	rho = state2dm(psi)
	trA = rho.ptrace(1)
	trA2 = trA * trA
	purity = trA2.tr()
	qprint('psi', psi)
	print("purity=%g" % purity)

	psi = build_state([qb10, qb01])
	rho = state2dm(psi)
	trA = rho.ptrace(1)
	trA2 = trA * trA
	purity = trA2.tr()
	qprint('psi', psi)
	print("purity=%g" % purity)


def problem_2_5():
	for N in range(4,5):
		print("N=%d =============" % N)
		psi = makeWstate(N)
		rho = state2dm(psi)
		trN = rho.ptrace(list(range(N-1)))
		trN2 = trN * trN
		qprint('psi', psi)
		qprint('trN', trN)
		qprint('trN2', trN2)
		purity = trN2.tr()
		print("purity=%g\n"  % purity)
	
def build_state(l):
	N = len(l)
	state = (1.0 / math.sqrt(N)) * sum(l)
	return state

def problem_4_1():
	psi0 = normalizer(2) * (qb000 + qb111)
	psi1 = normalizer(2) * (qb000 - qb111)
	qprint('psi0', psi0)
	qprint('psi1', psi1)
	rho0 = state2dm(psi0)
	rho1 = state2dm(psi1)
	rhoA0 = rho0.ptrace(0)
	rhoA1 = rho1.ptrace(0)
	qprint('rhoA0', rhoA0)
	qprint('rhoA1', rhoA1)

def problem_4_2():
	psi0 = normalizer(2) * (qb000 + qb111)
	psi1 = normalizer(2) * (qb000 - qb111)
	qprint('psi0', psi0)
	qprint('psi1', psi1)
	rho0 = state2dm(psi0)
	rho1 = state2dm(psi1)
	rhoAB0 = rho0.ptrace([0,1])
	rhoAB1 = rho1.ptrace([0,1])
	qprint('rhoAB0', rhoAB0)
	qprint('rhoAB1', rhoAB1)

def problem_5_2():
	psi0 = normalizer(2) * (qb00 + qb11)
	psi1 = normalizer(2) * (qb00 - qb11)
	psi2 = normalizer(2) * (qb01 + qb10)
	psi3 = normalizer(2) * (qb01 - qb10)

	qprint('psi', psi0)

	#U0 = tensor(qeye(2), qeye(2))
	#U1 = tensor(qeye(2), sigmax())
	#U2 = tensor(qeye(2), sigmaz())
	#U3 = tensor(qeye(2), sigmax()*sigmaz())

	U0 = tensor(qeye(2), qeye(2))
	U1 = tensor(qeye(2), sigmaz())
	U2 = tensor(qeye(2), sigmax())
	U3 = tensor(qeye(2), sigmax()*sigmaz())

	#U0 = tensor(qeye(2), qeye(2))
	#U1 = tensor(qeye(2), sigmax())
	#U2 = tensor(qeye(2), sigmax()*sigmaz())
	#U3 = tensor(qeye(2), sigmaz())

	phi0 = U0 * psi0
	phi1 = U1 * psi1
	phi2 = U2 * psi2
	phi3 = U3 * psi3

	qprint('phi0', phi0)
	qprint('phi1', phi1)
	qprint('phi2', phi2)
	qprint('phi3', phi3)

	#print("delta0={}".format(psi0-phi0))
	#print("delta1={}".format(psi1-phi1))
	#print("delta2={}".format(psi2-phi2))
	#print("delta3={}".format(psi3-phi3))

def main():
	#test01()
	#test02()
	#test03()
	problem_2_1()
	#problem_2_2()
	#problem_2_3()
	#problem_2_4()
	#problem_2_5()
	#problem_4_1()
	#problem_4_2()
	#problem_5_2()

if __name__ == "__main__":
	main()


