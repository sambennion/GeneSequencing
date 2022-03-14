#!/usr/bin/python3

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

class Cell(object):
		cost = 0
		direction = None
		def __init__(self):
			
			Cell.cost = 0
			Cell.direction = None

class GeneSequencing:

	def __init__( self ):
		pass

# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
# how many base pairs to use in computing the alignment


	def get_sequences(self, table, back_pointers, seq1, seq2):
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
	def unrestricted(self, seq1, seq2, align_length):
		
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
		align1, align2 = self.get_sequences(table, back_pointers, seq1, seq2)
		return table[m][n], align1, align2 #final score

	def restricted(self, seq1, seq2, align_length):
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
		align1, align2 = self.get_sequences(table, back_pointers, seq1, seq2)
		return table[m][n], align1, align2 #final score
				

	

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
			score, alignment1, alignment2 = self.unrestricted(seq1, seq2, align_length)
		else:
			score, alignment1, alignment2 = self.restricted(seq1, seq2, align_length)


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
	
