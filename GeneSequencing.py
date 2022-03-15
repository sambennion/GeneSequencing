#!/usr/bin/python3

# from tkinter.tix import MAX
import time
from operator import le
from re import A

from numpy import left_shift
from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT6':
	from PyQt6.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time
import random

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1



class GeneSequencing:

	def __init__( self ):
		pass

# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
# how many base pairs to use in computing the alignment

	#Time complexity: O(n) where n is the length of the larger of seq1 and seq2
	#Space Complexity: O(kn) as this uses an array of back_pointers size kn
	def get_sequences_banded(self, back_pointers, seq1, seq2):
		m = len(seq2)
		n = len(seq1)
		k = 2 * MAXINDELS + 1
		j = k-4+n-m

		

		align1 = ""
		align2 = ""


		while m > 3:
			if back_pointers[m][j] == "up":
				align2 = seq2[m-1] + align2
				align1 = "-" + align1
				m = m - 1
				j = j + 1
			elif back_pointers[m][j] == "lef":
				align1 = seq1[n-1] + align1
				align2 = "-" + align2
				j = j - 1
				n = n - 1
			else:
				align1 = seq1[n-1] + align1
				align2 = seq2[m-1] + align2
				m = m - 1
				n = n - 1
		#This will just run 3 times because m <=3
		begin_align1, begin_align2 = self.get_sequences(back_pointers, seq1[:n], seq2[:m])
		#combine first 3 alignment with the rest of the alignment
		align1 = begin_align1 + align1
		align2 = begin_align2 + align2
		return align1, align2


	def get_sequences(self, back_pointers, seq1, seq2):
		n = len(seq1)
		m = len(seq2)
		align1 = ""
		align2 = ""
		while n != 0 and m != 0:
			if back_pointers[m][n] == "up":
				align1 = "-" + align1
				align2 = seq2[m-1] + align2
				m = m-1
			elif back_pointers[m][n] == "lef":
				align1 = seq1[n-1] + align1
				align2 = "-" + align2
				n = n - 1
			else:
				align1 = seq1[n-1] + align1
				align2 = seq2[m-1] + align2
				m = m - 1
				n = n - 1
		return align1, align2


	#Time Complexity: O(nm)
	#Space Complexity: O(nm)
	def unrestricted(self, seq1, seq2):
		
		n = len(seq1)
		m = len(seq2)
		#O((n+1) * O(m+1)) or simply O(nm) space complexity
		table = [[0 for i in range(n+1)] for j in range(m+1)] 
		#Also O(nm) space complexity. I could've made each node in the table an object that contains a score and back_pointer, but it doesn't simplify it anymore.
		back_pointers  = [[0 for i in range(n+1)] for j in range(m+1)] 
		# print(table)
		for i in range(m+1):
			for j in range(n+1):
				if i == 0 or j == 0: #top or left-most column
					table[i][j] = (INDEL * i) + (INDEL * j)
				else:
					matched = seq1[j-1] == seq2[i-1]
					#take the min score from the surrounding
					#back pointers are being set at each phase if the condition is matched, but only the last one sticks.
					if matched:
						min = table[i-1][j-1] + MATCH #diagnal cell - Last priority in case of tie
					else:
						min = table[i-1][j-1] + SUB
					back_pointers[i][j] = 'di'
					if min >= table[i-1][j] + INDEL: #top cell - Second priority in case of tie
						min = table[i-1][j] + INDEL
						back_pointers[i][j] = 'up'

					if min >= table[i][j-1] + INDEL: #left cell - Top priority in case of tie
						min = table[i][j-1] + INDEL
						back_pointers[i][j] = 'lef'
					 
					table[i][j] = min
		align1, align2 = self.get_sequences(back_pointers, seq1, seq2)
		# print(table)
		# time.sleep(100)
		return table[m][n], align1, align2 #final score

	def restricted(self, seq1, seq2):
		# print(seq1)
		# print(seq2)
		n = len(seq1)
		m = len(seq2)
		k = 2 * MAXINDELS + 1

		if abs(n - m) > MAXINDELS:

			return math.inf, "No Alignment Possible", "No Alignment Possible"
		offset = 0  #This will help find the start of string1
		table = [[0 for i in range(k)] for j in range(m+1)] 
		back_pointers  = [[0 for i in range(k)] for j in range(m+1)] 


		for i in range(m + 1):

			if i > 3:#this will keep track of the offset for the band and increment it
				offset = offset + 1
			for j in range(k):
				#Set cell to inf if outside band for places where band is less than k.
				if (i == 0 and j > 3) or (i == 1 and j > 4) or (i == 2 and j > 5) or (i == m and j > 4) or (i == m-1 and j > 5) or (j + offset > n): #honestly, I may have put too many conditions here, but it works.
						# print("i = " + str(i) + " j = " + str(j) + " j + offset> n = " + str(j + offset> n))
						table[i][j] = math.inf

				elif i == 0 or j + offset == 0: #top or left-most column
					table[i][j] = (INDEL * i) + (INDEL * j)
				else:
					left_shift_up = 0
					if offset > 0:
						left_shift_up = 1 #After there's an offset, table comparisons need to be left shifted
					# print(offset + j)
					# print(seq1[j+offset-1])
					matched = seq1[j+offset-1] == seq2[i-1] #True if chars matchs

					if matched:
						min = table[i-1][j + left_shift_up - 1] + MATCH #diagnal cell - Last priority in case of tie
					else:
						min = table[i-1][j-1 + left_shift_up] + SUB
					back_pointers[i][j] = 'di'
					if j != k-1 and min >= table[i-1][j + left_shift_up] + INDEL: #top cell - Second priority in case of tie
						min = table[i-1][j + left_shift_up] + INDEL
						back_pointers[i][j] = 'up'

					if j != 0 and min >= table[i][j-1] + INDEL: #left cell - Top priority in case of tie
						min = table[i][j-1] + INDEL
						back_pointers[i][j] = 'lef'
					table[i][j] = min
					# print(table)
		
		align1, align2 = self.get_sequences_banded(back_pointers, seq1, seq2)
		# print(table)
		# time.sleep(100)
		# if table[m][k-4] == -474:
		# 	print(len(table))
		# 	time.sleep(100)

		return table[m][k-4+n-m], align1, align2 #k-4 seems to be the default if they are the same length. The n-m will change it based on how close together they are.


	

	def align( self, seq1, seq2, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length
		# print("Debug: " + seq1 + ", " + seq2)
		seq1 = seq1[0:align_length]
		seq2 = seq2[0:align_length]
		n = len(seq1)
		m = len(seq2)
		# print(n)
		# print(m)
		# opt = [[0 for i in range(n+1)] for j in range(m+1)]
		score = 0
		if not banded:
			score, alignment1, alignment2 = self.unrestricted(seq1, seq2)
		else:
			score, alignment1, alignment2 = self.restricted(seq1, seq2)

		alignment1 = alignment1[0:100]
		alignment2 = alignment2[0:100]
		# print("Score = " + str(score))
		
	


		
###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
		# score = random.random()*100
		# alignment1 = 'abc-easy  DEBUG:({} chars,align_len={}{})'.format(
		# 	len(seq1), align_length, ',BANDED' if banded else '')
		# alignment2 = 'as-123--  DEBUG:({} chars,align_len={}{})'.format(
		# 	len(seq2), align_length, ',BANDED' if banded else '')
###################################################################################################

		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}
	
