#!/home/ubeuser/anaconda3/bin/python

#####################################################################
# WEEK 6 - LAB
#####################################################################

import math
import numpy as np
#import shipy as sp
#import matplotlib as plt
from qutip import *
from udr import *



def exercise_1():
	print("== EXERCISE 1")
	rho0 = state2dm(qb0)
	rho1 = state2dm(qbHp)
	M00 = state2dm(qb0)
	Mpp = state2dm(qbHp)
	qprint('rho0', rho0)
	qprint('rho1', rho1)
	qprint('M00', M00)
	qprint('Mpp', Mpp)

	# FIRST ROUND
	rho = rho0 # This is what Alice sends to Bob
	print("First round - Alice sends {}".format(rho))
	p0_r1 = (M00 * rho).tr()
	print("Probability of revealing b=0 and Bob accepts={}".format(p0_r1))
	p1_r1 = (Mpp * rho).tr()
	print("Probability of revealing b=1 and Bob accepts={}".format(p1_r1))

	# SECOND ROUND
	rho = rho1 # This is what Alice sends to Bob
	print("Second round - Alice sends {}".format(rho))
	p0_r2 = (M00 * rho).tr()
	print("Probability of revealing b=0 and Bob accepts={}".format(p0_r2))
	p1_r2 = (Mpp * rho).tr()
	print("Probability of revealing b=1 and Bob accepts={}".format(p1_r2))
	
	p_cheat = max(abs(p1_r1), abs(p0_r2))
	print("Alice can cheat with probability up to {}".format(p_cheat))




def exercise_2():
	print("== EXERCISE 2")
	psi0 = qb0
	psi1 = qbHp
	m0 = qb0
	m1 = qbHp
	for theta in range(0,91,5):
		s = math.sin(theta*math.pi/180.0)
		c = math.cos(theta*math.pi/180.0)
		psi = (c*psi0) + (s*psi1)
		psi = psi.unit() # UBEDEBUG
		p0 = abs(psi.overlap(m0))**2
		p1 = abs(psi.overlap(m1))**2
		#qprint('psi', psi.dag())
		print("tdeg=%g p0=%.4f p1=%.4f |psi|=%.4f" % (theta, p0, p1, psi.norm()))


def exercise_3():
	print("== EXERCISE 3 ... damn !!")

def exercise_4():
	print("== EXERCISE 4")


	
	


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
	exercise_1()
	exercise_2()
	exercise_3()
	exercise_4()

if __name__ == "__main__":
	main()


