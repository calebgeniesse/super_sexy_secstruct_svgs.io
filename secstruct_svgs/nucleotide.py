class Nucleotide:
	def __init__( self, one_letter_code, seqpos ):
		self.name = one_letter_code
		self.seqpos = seqpos
		self.x = 0
		self.y = 0
		self.bp_partner = None
		