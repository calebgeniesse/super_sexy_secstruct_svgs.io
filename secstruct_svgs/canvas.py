import svgwrite
from util import color, loop_interpolate, interpolate, flatten
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
		self.list_of_coaxialities = []

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
	
	def draw_line( self, loc1, loc2, stroke=svgwrite.rgb(10, 10, 16, '%') ):
		self.dwg.add(self.dwg.line(loc1, loc2, stroke=stroke ) )

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
		
		center1 = [bp1.x,bp1.y]
		center2 = [bp2.x,bp2.y]
		
		# Aim is 80%
		frac = 0.8
		f1 = ( 1.0 + frac ) / 2
		f2 = ( 1.0 - frac ) / 2

		#beg = [ f1 * center1[0] + f2 * center2[0], f1 * center1[1] + f2 * center2[1] ]
		#end = [ f2 * center1[0] + f1 * center2[0], f2 * center1[1] + f1 * center2[1] ]
		beg = [ center1[0] + f1*(center2[0]-center1[0]), center1[1] + f1*(center2[1]-center1[1]) ]
		end = [ center1[0] + f2*(center2[0]-center1[0]), center1[1] + f2*(center2[1]-center1[1]) ]

		#self.draw_line( (bp1.x + 10, bp1.y - self.font_height_offset), 
		#				(bp1.x + 25, bp1.y - self.font_height_offset) )
		self.draw_line( beg, end )
				
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
		
		# We need an unstructured stem.numbers for a second.
		stem_numbers = flatten( stem.numbers )

		# Set initial positions
		x_offset = stem.coordinate_frame.position.x
		y = stem.coordinate_frame.position.y
		stem.nucleotides[ min(stem_numbers) ].x = x_offset
		stem.nucleotides[ min(stem_numbers) ].y = y
		# If there is an i-1, this isn't root
		if min(stem_numbers) - 1 in self.nucleotides.keys():
			stem.nucleotides[ min(stem_numbers) ].ref_nt = self.nucleotides[ min(stem_numbers) - 1 ]

		for bp_idx, bp in enumerate(stem.base_pairs): 
			if bp.nt1.seqpos == min(stem_numbers) - 1:
				bp.nt1.dy = bp.nt1.y - bp.nt1.ref_nt.y
			else:	
				bp.nt1.dy = self.bp_offset_height * stem.coordinate_frame.orientation
			
			if bp.nt2.seqpos == min(stem_numbers) - 1:
				bp.nt2.dx = bp.nt2.x - bp.nt2.ref_nt.x
			else:	
				bp.nt2.dx = self.bp_offset_width
			
			bp.nt1.update_absolute_coords()
			bp.nt2.update_absolute_coords()

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
			loop.coordinate_frame.position.y += loop.stem1.coordinate_frame.orientation * (len(loop.stem1.base_pairs) * self.bp_offset_height - self.bp_offset_height)
			# Extra for alignment (may change with font size?)
			loop.coordinate_frame.position.y += loop.stem1.coordinate_frame.orientation *2
			# Extra for flipped loops( no explanation)
			if loop.stem1.coordinate_frame.orientation == -1:
				loop.coordinate_frame.position.y -= 2 * self.bp_offset_height

			
			# Set end to start, but add bp_offset_width DEPENDING on stem orientation
			loop.coordinate_frame.position2.x  = loop.coordinate_frame.position.x
			loop.coordinate_frame.position2.x += self.bp_offset_width * loop.stem1.coordinate_frame.orientation
			loop.coordinate_frame.position2.y  = loop.coordinate_frame.position.y

			# Gentler arc for greater loop separations 
			# 0.75 for bp_offset_width, i.e. apicals
			# scale down to 0.3 for stem_offset_width 
			salient_difference = max( abs(loop.coordinate_frame.position.x - loop.coordinate_frame.position2.x), abs(loop.coordinate_frame.position.y - loop.coordinate_frame.position2.y) )
			total_loop_fraction = 0.75 - 0.45 * (salient_difference-self.bp_offset_width) / (self.stem_offset_width-self.bp_offset_width) 

			x1, y1 = loop.coordinate_frame.position.x,  loop.coordinate_frame.position.y
			x2, y2 = loop.coordinate_frame.position2.x, loop.coordinate_frame.position2.y
			
			# We are GUARANTEED that min(numbers)-1 exists.
			loop.nucleotides[ min(loop.numbers) ].ref_nt = self.nucleotides[ min(loop.numbers) - 1 ]
			
			for seqpos, nt in loop.nucleotides.iteritems():
				if seqpos == min(loop.numbers): continue
				nt.ref_nt = loop.nucleotides[ seqpos - 1 ]
				
			for seqpos, loop_nt in loop.nucleotides.iteritems():
				loop_idx = seqpos - min(loop.nucleotides.keys())
				fraction_done_with_loop = (float(loop_idx)+0.5) / float(len( loop.numbers ))
				# Note: loop_interpolate really ought to think about stem orientation
				# Maybe add 180 for flipped (-1) stems
				# Note that here we're for sure apical
				add_angle = 180 if loop.stem1.coordinate_frame.orientation == -1 else 0

				print "#### Fraction done with loop is: ", fraction_done_with_loop
				loop_nt.x, loop_nt.y = loop_interpolate( x1,y1,x2,y2, total_loop_fraction, fraction_done_with_loop, add_angle=add_angle )
				print "#### Placing loop nucleotide ", loop_nt.name, loop_idx, " at ", loop.loop_nt.x, loop_nt.y
				
				loop_nt.update_relative_coords()
						
						

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

			# Above-described algorithm sucks. Let's do something better:
			# instead of using stem coordinate frames, use the bonded NT.
			#print loop.stem1.nucleotides.keys()
			#print loop.stem2.nucleotides.keys()
			#print loop.numbers
			# Goes from small to large whichever stem it is... hrm.
			if min(loop.numbers)-1 in loop.stem1.nucleotides.keys():
				loop.coordinate_frame.position.x = loop.stem1.nucleotides[min(loop.numbers)-1].x
				loop.coordinate_frame.position.y = loop.stem1.nucleotides[min(loop.numbers)-1].y
				loop.coordinate_frame.position2.x = loop.stem2.nucleotides[max(loop.numbers)+1].x
				loop.coordinate_frame.position2.y = loop.stem2.nucleotides[max(loop.numbers)+1].y
			else: # max(loop.numbers)+1 in loop.stem1.nucleotides.keys()
				loop.coordinate_frame.position.x = loop.stem2.nucleotides[min(loop.numbers)-1].x
				loop.coordinate_frame.position.y = loop.stem2.nucleotides[min(loop.numbers)-1].y
				loop.coordinate_frame.position2.x = loop.stem1.nucleotides[max(loop.numbers)+1].x
				loop.coordinate_frame.position2.y = loop.stem1.nucleotides[max(loop.numbers)+1].y


			if loop.stem1.coordinate_frame.orientation == 1 and loop.stem2.coordinate_frame.orientation == 1:
				loop.coordinate_frame.position.y -= self.bp_offset_height * 6
				loop.coordinate_frame.position2.y -= self.bp_offset_height * 6
			elif loop.stem1.coordinate_frame.orientation == -1 and loop.stem2.coordinate_frame.orientation == -1:
				loop.coordinate_frame.position.y += self.bp_offset_height * 2
				loop.coordinate_frame.position2.y += self.bp_offset_height * 2
			elif loop.stem1.coordinate_frame.orientation == 1 and loop.stem2.coordinate_frame.orientation == -1:
				loop.coordinate_frame.position.x += self.bp_offset_width 
				loop.coordinate_frame.position2.x += self.bp_offset_width 
			else: # loop.stem1.coordinate_frame.orientation == -1 and loop.stem2.coordinate_frame.orientation == 1:
				loop.coordinate_frame.position.x -= self.bp_offset_width
				loop.coordinate_frame.position2.x -= self.bp_offset_width 

			x1, y1 = loop.coordinate_frame.position.x,  loop.coordinate_frame.position.y
			x2, y2 = loop.coordinate_frame.position2.x, loop.coordinate_frame.position2.y
			
			# Gentler arc for greater loop separations 
			# 0.75 for bp_offset_width, i.e. apicals
			# scale down to 0.3 for stem_offset_width 
			salient_difference = max( abs(x1-x2), abs(y1-y2) )
			total_loop_fraction = 0.75 - 0.45 * (salient_difference-self.bp_offset_width) / (self.stem_offset_width-self.bp_offset_width) 
		
			# AMW TEMP: we are right now moving towards a system where a loop already knows its position(s)
			# per nucleotide but we aren't there yet. So for now, note that loop_nt is suddenly a key to a dict
			# of Nucleotides -- not a number.
			for loop_key, loop_nt in loop.nucleotides.iteritems():
				#loop_idx, loop_nt in enumerate(loop.nucleotides.keys()):

				# silly update.........
				loop_idx = loop_key - min(loop.nucleotides.keys())

				fraction_done_with_loop = (float(loop_idx)+0.5) / float(len( loop.numbers ))
				# What if we could loop_interpolate 180 degrees for a straight thing? Maybe this is possible!
				#loop.nucleotides[loop_nt].x = interpolate(x1, x2, fraction_done_with_loop )
				#loop.nucleotides[loop_nt].y = interpolate(y1, y2, fraction_done_with_loop )
				# + to - means +90
				# - to + means +270 
				add_angle = 0
				loop_sense = 1
				if loop.stem1.coordinate_frame.orientation == -1 and loop.stem2.coordinate_frame.orientation == 1:
					add_angle = 90
				elif loop.stem1.coordinate_frame.orientation == 1 and loop.stem2.coordinate_frame.orientation == -1:
					add_angle = 270
				elif loop.stem1.coordinate_frame.orientation == 1 and loop.stem2.coordinate_frame.orientation == 1:
					add_angle = 0
					loop_sense = -1 # clockwise, not ccw
				else: # -1, -1
					add_angle = 180

				print "#### Fraction done with loop is: ", fraction_done_with_loop
				loop_nt.x, loop_nt.y = loop_interpolate( x1,y1,x2,y2, total_loop_fraction, fraction_done_with_loop, add_angle=add_angle, loop_sense=loop_sense )
				print "#### Placing loop nucleotide ", loop_nt.name, loop_idx, " at ", loop_nt.x, loop_nt.y

		self.loops.append(loop)
		for seqpos, loop_nt in loop.nucleotides.iteritems(): self.nucleotides[ seqpos ] = loop_nt


	def set_stems_coaxial( self, idx1, idx2 ):
		"""
		Manipulate the coordinate frame of stem idx2 to be coaxial to idx1.
		This is the first time a concept of orientation has come into play, i.e.
		'at x,y AND GOING UP' vs. down.

		Note that there are multiple ways to set two stems coaxial -- many orientations etc.
		Especially with nonlocal coaxialities, as it were, these matter quite a lot. Not
		bothering yet.
		"""

		self.list_of_coaxialities.append( [idx1, idx2] )

		# Uh oh! We had committed to this "single coordinate frame per stem" idea. But now we really want
		# to interact with nucleotides directly!

		cf1 = self.stems[idx1].coordinate_frame
		cf2 = self.stems[idx2].coordinate_frame
		
		# idx2 is opposite of idx1 in orientation. Not necessarily going to be true in
		# all kinds of coaxiality, but in this first, simplest kind...
		cf2.orientation = -1 * cf1.orientation

		cf2.position.y = cf1.position.y + cf2.orientation * ( self.bp_offset_height + 1 )
		cf2.position.x = cf1.position.x + self.bp_offset_width # coaxiality, keeping in mind the flipped 'sense'

		# We can learn a lot from add_stem...
		x_offset = cf2.position.x
		y = cf2.position.y
		for bp_idx, bp in enumerate(self.stems[idx2].base_pairs): 
			# need to compute some offet? ...
			print "orientation is ", self.stems[idx2].coordinate_frame.orientation
			y_offset = bp_idx * self.bp_offset_height * self.stems[idx2].coordinate_frame.orientation
			bp.nt1.x = x_offset
			bp.nt1.y = y + y_offset
			bp.nt2.x = x_offset + self.bp_offset_width * cf2.orientation # Flip for flipped.
			bp.nt2.y = y + y_offset

		# Oh! Update any dependent loops!
		for loop in self.loops:
			if self.stems[idx2] == loop.stem1 and loop.stem2 == None:
				# apical
				# at some point reduce duplication wtih add loop!!!!
				loop.coordinate_frame.position.x  = loop.stem1.coordinate_frame.position.x
				loop.coordinate_frame.position.y  = loop.stem1.coordinate_frame.position.y
				# Offset for stem length
				loop.coordinate_frame.position.y += len(loop.stem1.base_pairs) * self.bp_offset_height - self.bp_offset_height
				# Extra for alignment (may change with font size?)
				loop.coordinate_frame.position.y += 2
			
				# Set end to start, but add bp_offset_width DEPENDING on stem orientation
				loop.coordinate_frame.position2.x  = loop.coordinate_frame.position.x
				loop.coordinate_frame.position2.x += self.bp_offset_width * loop.stem1.coordinate_frame.orientation
				loop.coordinate_frame.position2.y  = loop.coordinate_frame.position.y

				#Technically, handle tail here or have a better condition
			else: # junction
				pass

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
		for loop_key, loop_nt in apical_loop.nucleotides.iteritems():
			self.draw_nt(loop_nt)
			
	def draw_junction_loop( self, junction_loop ):
		for loop_key, loop_nt in junction_loop.nucleotides.iteritems():
			self.draw_nt(loop_nt)
			
	def draw_sequence_line( self, pos1, pos2 ):
		"""
		Draw a line between consecutive nucleotides. Thin and grey. 
		Figure out where the subjective centers of the nts are ( depends on font but likely about 
		5-6 px up and right of (x,y) ) and draw about 80% of that vector.
		"""

		center1 = [ self.nucleotides[ pos1 ].x + self.font_height_offset, self.nucleotides[ pos1 ].y - self.font_height_offset ]
		center2 = [ self.nucleotides[ pos2 ].x + self.font_height_offset, self.nucleotides[ pos2 ].y - self.font_height_offset ]

		#beg = [ 0.9 * center1[0] + 0.1 * center2[0], 0.9 * center1[1] + 0.1 * center2[1] ]
		#end = [ 0.1 * center1[0] + 0.9 * center2[0], 0.1 * center1[1] + 0.9 * center2[1] ]

		# Aim is 60%
		frac = 0.5
		f1 = ( 1.0 + frac ) / 2
		f2 = ( 1.0 - frac ) / 2

		#beg = [ f1 * center1[0] + f2 * center2[0], f1 * center1[1] + f2 * center2[1] ]
		#end = [ f2 * center1[0] + f1 * center2[0], f2 * center1[1] + f1 * center2[1] ]
		beg = [ center1[0] + f1*(center2[0]-center1[0]), center1[1] + f1*(center2[1]-center1[1]) ]
		end = [ center1[0] + f2*(center2[0]-center1[0]), center1[1] + f2*(center2[1]-center1[1]) ]

		self.draw_line( beg, end, 'gray' )

	def render( self ):
		for stem in self.stems: self.draw_stem( stem )
		for apical in [ loop for loop in self.loops if loop.apical ]:
			self.draw_apical_loop( apical )
		for junction in [ loop for loop in self.loops if not loop.apical ]:
			self.draw_junction_loop( junction )

		for seqpos in self.nucleotides.keys():
			if seqpos + 1 in self.nucleotides.keys():
				self.draw_sequence_line( seqpos, seqpos + 1 )

		self.dwg.save()
