class Nucleotide:

	def update_absolute_coords( self, unused=None ):
		if self.ref_nt is None:
			self.x = self.dx
			self.y = self.dy
		else:
			self.x = self.ref_nt.x + self.dx
			self.y = self.ref_nt.y + self.dy

	def update_relative_coords( self, unused=None ):
		if self.ref_nt is None:
			self.x = self.dx
			self.y = self.dy
		else:
			self.dx = self.x - self.ref_nt.x
			self.dy = self.y - self.ref_nt.y

	def __init__( self, one_letter_code, seqpos ):
		self.name = one_letter_code
		self.seqpos = seqpos
		self.x = 0
		self.y = 0
		self.bp_partner = None
		# Note: we're not doing 'internal coordinates' but
		# truly just relative Cartesian coordinates.
		# Still helps for propagation.
		# Root is always first nt; no rerooting (yet?)
		# ref nt generally i-1, but can be just some BPed nt
		self.ref_nt = None
		self.dx = 0
		self.dy = 0
