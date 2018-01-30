#!/home/ubeuser/anaconda3/bin/python

#####################################################################
# WEEK 9 - LAB
#####################################################################

import math
import numpy as np
from random import randint
#import shipy as sp
#import matplotlib as plt
from qutip import *
from udr import *

DBG_PRINT_FLAG = False

def dbg_print(msg):
	if DBG_PRINT_FLAG:
		print(msg)

		
def exercise_1():
	print("== EXERCISE 1")
	phi = ((1.0/math.sqrt(3.0)) * qb0) + (math.sqrt(2.0/3.0) * qb1)
	rho = state2dm(phi)
	qprint('phi', phi)
	qprint('rho', rho)

	tr_rho = rho.tr()
	print("[1] Tr(rho) = {}".format(tr_rho))

	p0 = (qb0.dag() * rho * qb0).tr()
	print("[2] p(0) = {}".format(p0))

	p1 = (qb1.dag() * rho * qb1).tr()
	print("[3] p(1) = {}".format(p1))

	N_rho = 0.5 * rho
	print("[4] Tr(N_rho) = {}".format(N_rho.tr()))

	p0 = (qb0.dag() * N_rho * qb0).tr()
	print("[5] p(0) = {} ({})".format(p0, 1.0/6.0))

	p1 = (qb1.dag() * N_rho * qb1).tr()
	print("[6] p(1) = {} ({})".format(p1, 1.0/3.0))



def exercise_2_123():
	print("== EXERCISE 2-123")
	rho = ((1.0/3.0) * state2dm(qb0)) + ((2.0/3.0) * state2dm(qb1))
	qprint('rho', rho)
	N_rho = ((1.0/3.0) * state2dm(qb0)) + ((2.0/3.0) * state2dm(qb1))

	print("[1] Th(rho)={} Tr(N_rho)={}".format(rho.tr(), N_rho.tr()))

	p0R = (qb0.dag() * rho * qb0).tr()
	p1R = (qb1.dag() * rho * qb1).tr()
	p0N = (qb0.dag() * N_rho * qb0).tr()
	p1N = (qb1.dag() * N_rho * qb1).tr()
	print("[2] poR={} p1R={} p0N={} p1N={}".format(p0R,p1R,p0N,p1N))


def exercise_2_4():
	print("== EXERCISE 2-4")
	phi = None
	#qprint('phi-init', phi)
	for i in range(2):
		for j in range(2):
			for k in range(2):
				for l in range(2):
					#print("i={} j={} k={} l={}".format(i,j,k,l))
					bob = basis(2,i) * basis(2,j).dag()
					#qprint('bob', bob)
					eve = basis(2,k) * basis(2,l).dag()
					#qprint('eve', eve)
					temp = tensor(bob, eve)
					#qprint('temp', temp)
					if phi is None:
						phi = temp
					else:
						phi = phi + temp
					#qprint('phi', phi)
	qprint('PHI', phi)

	psi = None
	#qprint('psi-init', psi)
	for i in range(2):
		for j in range(2):
			for k in range(2):
				for l in range(2):
					#print("i={} j={} k={} l={}".format(i,j,k,l))
					bob = basis(2,j) * basis(2,i).dag()
					#qprint('bob', bob)
					eve = basis(2,k) * basis(2,l).dag()
					#qprint('eve', eve)
					temp = tensor(bob, eve)
					#qprint('temp', temp)
					if psi is None:
						psi = temp
					else:
						psi = psi + temp
					#qprint('psi', psi)
	qprint('PSI', psi)

	l_phi = phi.eigenenergies()
	print("l_phi={}".format(l_phi))

	es_phi = phi.eigenstates()
	print("es_phi={}".format(es_phi))

	


def test01():
	print("== TEST 01 ==")

def test02():
	print("== TEST 02 ==")
	
def test03():
	print("== TEST 03 ==")


def main():
	#test01()
	#test02()
	#test03()
	#exercise_1()
	#exercise_2_123()
	exercise_2_4()

if __name__ == "__main__":
	main()


