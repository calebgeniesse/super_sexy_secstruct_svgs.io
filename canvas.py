import svgwrite
from util import color, loop_interpolate, interpolate
from coordinate_frame import CoordinateFrame

class Canvas:
	def __init__( self, dwg, width_per_stem=80, width_per_bp=30, height_per_bp=20, orig=None):
	#	if orig == None:
	#		self.ctor( dwg, width_per_stem, width_per_bp, height_per_bp )
	#	else:
	#		self.copyctor( orig )
	#
	#def ctor( self, dwg, width_per_stem, width_per_bp, height_per_bp ):
		self.dwg = dwg
		self.stem_offset_width = width_per_stem
		self.bp_offset_width   = width_per_bp
		self.bp_offset_height  = height_per_bp
		self.font_height_offset = 3
		self.stems = []
		self.loops = []
		self.nucleotides = {}

	#def copyctor( self, orig ):
	#	self.dwg = orig.dwg
	#	self.stem_offset_width = orig.width_per_stem
	#	self.bp_offset_width   = orig.width_per_bp
	#	self.bp_offset_height  = orig.height_per_bp
	#	self.font_height_offset = 3
	#	self.stems = [ Stem( orig.stem ) for 
	#	self.loops = []
	#	self.nucleotides = {}

	
	def draw_text( self, text, loc, color ):
		self.dwg.add(self.dwg.text(text, insert=loc, fill=color ) )
	
	def draw_line( self, loc1, loc2 ):
		""" 
		Only supports black for now.
		"""
		self.dwg.add(self.dwg.line(loc1, loc2, stroke=svgwrite.rgb(10, 10, 16, '%')))

	def draw_nt( self, nt ):
		self.dwg.add(self.dwg.text(nt.name, insert=(nt.x, nt.y), fill=color(nt.name) ) )

		# numbers if multiple of 5
		# TODO: adjust size
		if nt.seqpos % 5 == 0: self.draw_text( nt.seqpos, (nt.x-15, nt.y+5), color(nt.name) )	
	
	def draw_bp( self, bp1, bp2 ):
		"""
		Needs some kind of line centering relative to the height of the characters.
		"""
		self.draw_nt( bp1 )
		self.draw_nt( bp2 )
		# Assumed horizontal. Don't have a good idea of how to figure out this line otherwise.
		# Probably draw about 80% of (x1,y1)-(x2,y2)
		self.draw_line( (bp1.x + 10, bp1.y - self.font_height_offset), 
						(bp1.x + 25, bp1.y - self.font_height_offset) )
		
	def add_stem( self, stem ):
		"""
		We override the default coordinate frame with one that uses a little more information
		about the current set of stems. Of course, this may in turn be updated if we find out
		information about coaxiality etc.
		"""

		# Previously we would re-construct a coordinate frame here. That seems dumb.
		stem.coordinate_frame.position.x = len(self.stems) * self.stem_offset_width
		self.stems.append( stem )
		for seqpos in stem.nucleotides.keys(): self.nucleotides[ seqpos ] = stem.nucleotides[ seqpos ]
		
		# Set initial positions

		# draw basepairs
		x_offset = stem.coordinate_frame.position.x
		y = stem.coordinate_frame.position.y
		for bp_idx, bp in enumerate(stem.base_pairs): 
			# need to compute some offet? ...
			print "orientation is ", stem.coordinate_frame.orientation
			y_offset = bp_idx * self.bp_offset_height * stem.coordinate_frame.orientation
			bp.nt1.x = x_offset
			bp.nt1.y = y + y_offset
			bp.nt2.x = x_offset + self.bp_offset_width
			bp.nt2.y = y + y_offset


	def add_loop( self, loop ):
		if loop.apical:
			"""
			y offsets don't make a ton of sense here. 
			"""
			# Just trace path from coordinate frame start, plus y offset for
			# height, to other side of stem.
			loop.coordinate_frame.position.x  = loop.stem1.coordinate_frame.position.x
			loop.coordinate_frame.position.y  = loop.stem1.coordinate_frame.position.y
			# Offset for stem length
			loop.coordinate_frame.position.y += len(loop.stem1.base_pairs) * self.bp_offset_height - self.bp_offset_height
			# Extra for alignment (may change with font size?)
			loop.coordinate_frame.position.y += 2
			
			# Set end to start, but add bp_offset_width
			loop.coordinate_frame.position2.x  = loop.coordinate_frame.position.x
			loop.coordinate_frame.position2.x += self.bp_offset_width
			loop.coordinate_frame.position2.y  = loop.coordinate_frame.position.y

			x1, y1 = loop.coordinate_frame.position.x,  loop.coordinate_frame.position.y
			x2, y2 = loop.coordinate_frame.position2.x, loop.coordinate_frame.position2.y
			for loop_idx, loop_nt in enumerate(loop.nucleotides.keys()):
				fraction_done_with_loop = (float(loop_idx)+0.5) / float(len( loop.numbers ))
				loop.nucleotides[loop_nt].x, loop.nucleotides[loop_nt].y = loop_interpolate( x1,y1,x2,y2, 0.75, fraction_done_with_loop )
			
		elif loop.tail:
			"""
			We can kind of make some choices here. When stem orientation is variable
			this will have to become a little more complicated, but at the moment we
			will be able to set the end of the coordinate frame based on length.
			"""
			# Is the stem the beginning or end of the tail?
			if min(loop.numbers)-1 in flatten(stem.numbers): # beginning
				# Right now that means: extend to 'right'
				# TODO: nt width magic number.
				loop.coordinate_frame.position.x  = loop.stem1.coordinate_frame.position.x 
				loop.coordinate_frame.position.y  = loop.stem1.coordinate_frame.position.y
				loop.coordinate_frame.position2.x = loop.coordinate_frame.position.x  + 7 * len(loop.numbers)
				loop.coordinate_frame.position2.y = loop.coordinate_frame.position2.y 
			if max(loop.numbers)+1 in flatten(stem.numbers): # beginning
				# Right now that means: extend to 'right'
				# TODO: nt width magic number.
				loop.coordinate_frame.position2.x = loop.stem1.coordinate_frame.position.x
				loop.coordinate_frame.position2.y = loop.stem1.coordinate_frame.position.y
				loop.coordinate_frame.position.x  = loop.coordinate_frame.position2.x  - 7 * len(loop.numbers)
				loop.coordinate_frame.position.y  = loop.coordinate_frame.position2.y 

			# TODO: for these, init positions.
		else:
			"""
			As long as we are just drawing this dumb left-to-right structure
			with no stacking, we should assume we are going from the right side of 
			stem1 to the left side of stem2. Note that we will really have to think
			in terms of the positions of individual NTs eventually.
			"""
			loop.coordinate_frame.position.x  = loop.stem1.coordinate_frame.position.x
			loop.coordinate_frame.position.x += self.bp_offset_width
			loop.coordinate_frame.position.y  = loop.stem1.coordinate_frame.position.y
			loop.coordinate_frame.position.y -= 10

			loop.coordinate_frame.position2.x  = loop.stem2.coordinate_frame.position.x
			loop.coordinate_frame.position2.y  = loop.stem2.coordinate_frame.position.y
			loop.coordinate_frame.position2.y -= 10

			x1, y1 = loop.coordinate_frame.position.x,  loop.coordinate_frame.position.y
			x2, y2 = loop.coordinate_frame.position2.x, loop.coordinate_frame.position2.y
		
			# AMW TEMP: we are right now moving towards a system where a loop already knows its position(s)
			# per nucleotide but we aren't there yet. So for now, note that loop_nt is suddenly a key to a dict
			# of Nucleotides -- not a number.
			#for loop_idx, loop_nt in enumerate(junction_loop.numbers):
			for loop_idx, loop_nt in enumerate(loop.nucleotides.keys()):
				# This fraction goes 0, 0.25, 0.5, 0.75
				# We need something more aggressive, esp. due to character alignment.
				fraction_done_with_loop = (float(loop_idx)+0.5) / float(len( loop.numbers ))
				loop.nucleotides[loop_nt].x = interpolate(x1, x2, fraction_done_with_loop )
				loop.nucleotides[loop_nt].y = interpolate(y1, y2, fraction_done_with_loop )
			
		self.loops.append(loop)
		for seqpos in loop.nucleotides.keys(): self.nucleotides[ seqpos ] = loop.nucleotides[ seqpos ]


	def set_stems_coaxial( self, idx1, idx2 ):
		"""
		Manipulate the coordinate frame of stem idx2 to be coaxial to idx1.
		This is the first time a concept of orientation has come into play, i.e.
		'at x,y AND GOING UP' vs. down.

		Note that there are multiple ways to set two stems coaxial -- many orientations etc.
		Especially with nonlocal coaxialities, as it were, these matter quite a lot. Not
		bothering yet.
		"""

		cf1 = self.stems[idx1].coordinate_frame
		cf2 = self.stems[idx2].coordinate_frame
		
		# idx2 is opposite of idx1 in orientation. Not necessarily going to be true in
		# all kinds of coaxiality, but in this first, simplest kind...
		cf2.orientation = -1 * cf1.orientation

		cf2.position.y = cf1.position.y + cf2.orientation * ( self.bp_offset_height + 1 )
		cf2.position.x = cf1.position.x # coaxiality


	def draw_stem( self, stem ):
		# draw basepairs
		x_offset = stem.coordinate_frame.position.x
		y = stem.coordinate_frame.position.y
		for bp_idx, bp in enumerate(stem.base_pairs): 
			# need to compute some offet? ...
			print "orientation is ", stem.coordinate_frame.orientation
			y_offset = bp_idx * self.bp_offset_height * stem.coordinate_frame.orientation
			
			self.draw_bp( bp.nt1, bp.nt2 )
			
	def draw_apical_loop( self, apical_loop ):
		for loop_idx, loop_nt in enumerate(apical_loop.nucleotides.keys()):
			self.draw_nt(apical_loop.nucleotides[loop_nt])
			
	def draw_junction_loop( self, junction_loop ):
		for loop_idx, loop_nt in enumerate(junction_loop.nucleotides.keys()):
			self.draw_nt(junction_loop.nucleotides[loop_nt])
			
	def render( self ):
		for stem in self.stems: self.draw_stem( stem )
		for apical in [ loop for loop in self.loops if loop.apical ]:
			self.draw_apical_loop( apical )
		for junction in [ loop for loop in self.loops if not loop.apical ]:
			self.draw_junction_loop( junction )
		self.dwg.save()
