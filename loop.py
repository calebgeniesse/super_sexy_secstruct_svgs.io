import util
from util import flatten

class Loop:
	def __init__( self, sequence, numbers, stem1, apical=True, tail=False, stem2=None ):
		"""
		Initialize a loop with a particular sequence. Loops are associated with stems
		(we will not create SS diagrams for bare loops); either one or two.
		"""
		self.sequence = sequence
		self.numbers = numbers
		self.stem1 = stem1
		self.apical = apical
		self.tail = tail
		self.stem2 = stem2
		# The coordinate frame of a loop is intrinsically more complicated, since
		# a loop is tracing out a path and even if stems can only be dull things,
		# we still have notions of beginning and end here.
		#
		# This will also be good practice for manipulating more complicated 
		# coordinate frames later.
		self.coordinate_frame = [[0,0],[0,0]]
	
