#!/home/ubeuser/anaconda3/bin/python

#####################################################################
# WEEK 9 - HOMEWORK
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

class pr_box:
	def __init__(self, r):
		self.r = r
		self.x = None
		self.y = None
		self.a = None
		self.b = None

	def set_x(self, x):
		if not ((x == 0) or (x == 1)):
			raise RuntimeError("Invalid input: {}".format(x))
		if not self.x is None:
			raise RuntimeError("Input X already set")
		self.x = x
		self.__calc_outputs()

	def set_y(self, y):
		if not ((y == 0) or (y == 1)):
			raise RuntimeError("Invalid input: {}".format(y))
		if not self.y is None:
			raise RuntimeError("Input Y already set")
		self.y = y
		self.__calc_outputs()

	def get_a(self):
		if not self.__has_outputs():
			raise RuntimeError("PR box output not available")
		return self.a

	def get_b(self):
		if not self.__has_outputs():
			raise RuntimeError("PR box output not available")
		return self.b

	def __calc_outputs(self):
		if self.__has_inputs() and not self.__has_outputs():
			self.a = self.r
			self.b = self.r ^ (self.x * self.y)
			dbg_print("PRB x={} y={} r={} a={} b={}".format(
				self.x,self.y,self.r,self.a,self.b))
		
	def __has_inputs(self):
		return not ((self.x is None) or (self.y is None))

	def __has_outputs(self):
		return not ((self.a is None) or (self.a is None))


def chsh_score(x, y, a, b):
	score = 0
	if (a ^ b) == (x * y):
		score = 1
	dbg_print("CHSH(x={} y={} a={} b={} score={}".format(x,y,a,b,score))
	return score
		
def exercise_1():
	print("== EXERCISE 1")


def exercise_2_1():
	print("== EXERCISE 2.1")
	score = 0
	counter = 0.0
	for x in range(0,2):
		for y in range(0,2):
			for r in range(0,2):
				counter += 1.0
				prb = pr_box(r)
				prb.set_x(x)
				prb.set_y(1)
				box_a = prb.get_a()
				box_b = prb.get_b()
				a = x + r
				b = y * r
				score += chsh_score(x, y, a, b)
				print("cnt={} x={} y={} r={} temp_score={}".format(counter,x,y,r,score))
	result = score / counter
	print("counter={} score={} RESULT={}".format(counter, score, result))

def exercise_2_2():
	print("== exercise 2.2")
	score = 0
	counter = 0.0
	for x in range(0,2):
		for y in range(0,2):
			for r in range(0,2):
				counter += 1.0
				prb = pr_box(r)
				prb.set_x(1)
				prb.set_y(y)
				box_a = prb.get_a()
				box_b = prb.get_b()
				a = x ^ r
				b = y * r
				score += chsh_score(x, y, a, b)
				print("cnt={} x={} y={} r={} temp_score={}".format(counter,x,y,r,score))
	result = score / counter
	print("counter={} score={} result={}".format(counter, score, result))

def exercise_2_3():
	print("== exercise 2.3")
	score = 0
	counter = 0.0
	for x in range(0,2):
		for y in range(0,2):
			for r in range(0,2):
				counter += 1.0
				prb = pr_box(r)
				prb.set_x(x)
				prb.set_y(y)
				box_a = prb.get_a()
				box_b = prb.get_b()
				a = r
				b = (x*y) ^ r
				score += chsh_score(x, y, a, b)
				print("cnt={} x={} y={} r={} temp_score={}".format(counter,x,y,r,score))
	result = score / counter
	print("counter={} score={} result={}".format(counter, score, result))

			

	


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
	exercise_2_1()
	exercise_2_2()
	exercise_2_3()
	exercise_4()

if __name__ == "__main__":
	main()


