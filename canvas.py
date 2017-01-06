import svgwrite
from util import color, loop_interpolate

class Canvas:
	def __init__( self, dwg, width_per_stem=80, width_per_bp=30, height_per_bp=20):
		self.dwg = dwg
		self.stem_offset_width = width_per_stem
		self.bp_offset_width   = width_per_bp
		self.bp_offset_height  = height_per_bp
		self.font_height_offset = 3
		self.stems = []
		self.loops = []
	
	def draw_text( self, text, loc, color ):
		self.dwg.add(self.dwg.text(text, insert=loc, fill=color ) )
	
	def draw_line( self, loc1, loc2 ):
		""" 
		Only supports black for now.
		"""
		self.dwg.add(self.dwg.line(loc1, loc2, stroke=svgwrite.rgb(10, 10, 16, '%')))

	def draw_nt( self, text, loc ):
		self.dwg.add(self.dwg.text(text, insert=loc, fill=color(text) ) )
	
	def draw_bp( self, bp1, bp2, loc ):
		"""
		Needs some kind of line centering relative to the height of the characters.
		"""
		self.draw_nt( bp1, loc )
		self.draw_nt( bp2, (loc[0]+self.bp_offset_width, loc[1]) )
		self.draw_line( (loc[0] + 10, loc[1] - self.font_height_offset), 
						(loc[0] + 25, loc[1] - self.font_height_offset) )

	def add_stem( self, stem ):
		"""
		We override the default coordinate frame with one that uses a little more information
		about the current set of stems. Of course, this may in turn be updated if we find out
		information about coaxiality etc.
		"""
		stem.coordinate_frame = CoordinateFrame( [ len(self.stems) * self.stem_offset_width, 0 ], 0 )
		self.stems.append( stem )

	def add_loop( self, loop ):
		if loop.apical:
			"""
			y offsets don't make a ton of sense here. 
			"""
			# Just trace path from coordinate frame start, plus y offset for
			# height, to other side of stem.
			loop.coordinate_frame.position.x  = loop.stem1.coordinate_frame.position.x
			loop.coordinate_frame.position.y  = loop.stem1.coordinate_frame.position.x
			# Offset for stem length
			loop.coordinate_frame.position.y += len(loop.stem1.base_pairs) * self.bp_offset_height - self.bp_offset_height
			# Extra for alignment (may change with font size?)
			loop.coordinate_frame.position.y += 2
			
			# Set end to start, but add bp_offset_width
			loop.coordinate_frame.position2.x  = loop.coordinate_frame.position.x
			loop.coordinate_frame.position2.x += self.bp_offset_width
			loop.coordinate_frame.position2.y  = loop.coordinate_frame.position.y
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

		self.loops.append(loop)


	def draw_stem( self, stem ):
		# draw basepairs
		x_offset = stem.coordinate_frame.position.x
		y = stem.coordinate_frame.position.y
		for bp_idx, [bp1, bp2] in enumerate(stem.base_pairs): 
			# need to compute some offet? ...
			y_offset = bp_idx * self.bp_offset_height
			self.draw_bp( bp1, bp2, (x_offset, y + y_offset) )

			# numbers if multiple of 5
			# TODO: adjust size
			# Also TODO: have these occur independent of stem vs loop
			# but right now, convenient
			num1, num2 = stem.numbers[bp_idx]
			if num1 % 5 == 0: self.draw_text( num1, (x_offset-15, y+y_offset+5), color(bp1) )
			if num2 % 5 == 0: self.draw_text( num2, (x_offset+self.bp_offset_width+10, y+y_offset+5), color(bp2) )
	
	def set_stems_coaxial( self, idx1, idx2 ):
		"""
		Manipulate the coordinate frame of stem idx2 to be coaxial to idx1.
		This is the first time a concept of orientation has come into play, i.e.
		'at x,y AND GOING UP' vs. down.

		Note that there are multiple ways to set two stems coaxial -- many orientations etc.
		Especially with nonlocal coaxialities, as it were, these matter quite a lot. Not
		bothering yet.
		"""

		# idx2 is opposite of idx1 in orientation. Not necessarily going to be true in
		# all kinds of coaxiality, but in this first, simplest kind...
		self.stems[idx2].coordinate_frame.orientation = -1 * self.stems[idx2].coordinate_frame.orientation



	def draw_apical_loop( self, apical_loop ):
		"""
		Right now this looks identical to draw_junction_loop.
		But we could imagine different interpolation/shape behaviors
		are appropriate for each.
		"""
		#[[x1,y1],orientation,[x2,y2]] = apical_loop.coordinate_frame
		x1, y1 = apical_loop.coordinate_frame.position.x,  apical_loop.coordinate_frame.position.y
		x2, y2 = apical_loop.coordinate_frame.position2.x, apical_loop.coordinate_frame.position2.y
	
		for loop_idx, loop_nt in enumerate(apical_loop.numbers):
			fraction_done_with_loop = (float(loop_idx)+0.5) / float(len( apical_loop.numbers ))
			# Here's a good example of an apical-junction deviation.
			[x, y] = loop_interpolate( x1,y1,x2,y2, 0.75, fraction_done_with_loop )
			self.draw_nt(apical_loop.sequence[loop_idx], (x, y))
			
			# When stems have 'orientation' draw these 'away from' stem
			num = apical_loop.numbers[loop_idx]
			if num % 5 == 0: self.draw_text( num, (x, y+10), color(apical_loop.sequence[loop_idx]) )
		
	def draw_junction_loop( self, junction_loop ):
		#[[x1,y1],[x2,y2]] = junction_loop.coordinate_frame
		x1, y1 = junction_loop.coordinate_frame.position.x,  junction_loop.coordinate_frame.position.y
		x2, y2 = junction_loop.coordinate_frame.position2.x, junction_loop.coordinate_frame.position2.y
		
		for loop_idx, loop_nt in enumerate(junction_loop.numbers):
			# This fraction goes 0, 0.25, 0.5, 0.75
			# We need something more aggressive, esp. due to character alignment.
			fraction_done_with_loop = (float(loop_idx)+0.5) / float(len( junction_loop.numbers ))
			x = interpolate(x1, x2, fraction_done_with_loop )
			y = interpolate(y1, y2, fraction_done_with_loop )
			self.draw_nt(junction_loop.sequence[loop_idx], (x,-10))
			
			# When stems have 'orientation' draw these away in the direction 
			# perpendicular to the vector between the stems I guess.
			num = junction_loop.numbers[loop_idx]
			if num % 5 == 0: self.draw_text( num, (x, y-10), color(junction_loop.sequence[loop_idx]) )

	def render( self ):
		for stem in self.stems: self.draw_stem( stem )
		for apical in [ loop for loop in self.loops if loop.apical ]:
			self.draw_apical_loop( apical )
		for junction in [ loop for loop in self.loops if not loop.apical ]:
			self.draw_junction_loop( junction )
		self.dwg.save()
