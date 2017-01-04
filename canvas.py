import svgwrite
from util import color

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
		stem.coordinate_frame = [ len(self.stems) * self.stem_offset_width, 0 ]
		self.stems.append( stem )

	def add_loop( self, loop ):
		if loop.apical:
			# Just trace path from coordinate frame start, plus y offset for
			# height, to other side of stem.
			loop.coordinate_frame[0][0]  = loop.stem1.coordinate_frame[0]
			loop.coordinate_frame[0][1]  = loop.stem1.coordinate_frame[1]
			# Offset for stem length
			loop.coordinate_frame[0][1] += len(loop.stem1.base_pairs) * self.bp_offset_height 
			# Extra for alignment (may change with font size?)
			loop.coordinate_frame[0][1] += 2
			
			# Set end to start, but add bp_offset_width
			loop.coordinate_frame[1][0]  = loop.coordinate_frame[0][0]
			loop.coordinate_frame[1][0] += self.bp_offset_width
			loop.coordinate_frame[1][1]  = loop.coordinate_frame[0][1]
		else:
			pass



	def draw_stem( self, stem ):
		# draw basepairs
		x_offset = stem.coordinate_frame[0]
		y = stem.coordinate_frame[1]
		for bp_idx, [bp1, bp2] in enumerate(stem.base_pairs): 
			# need to compute some offet? ...
			y_offset = bp_idx * self.bp_offset_height
			self.draw_bp( bp1, bp2, (x_offset, y + y_offset) )

			# numbers if multiple of 5
			# TODO: adjust size
			num1, num2 = stem.numbers[bp_idx]
			if num1 % 5 == 0: self.draw_text( num1, (x_offset-10, y_offset+5), color(bp1) )
			if num2 % 5 == 0: self.draw_text( num2, (x_offset+self.bp_offset_width+10, y_offset+5), color(bp2) )
		
	def interpolate( self, v1, v2, frac ):
		return v1 + float(v2-v1)*frac

	def draw_apical_loop_reasonably_with_respect_to_stem( self, apical_loop ):
		"""
		y offsets don't make a ton of sense here. 
		"""
		x1 = apical_loop.coordinate_frame[0][0]
		x2 = apical_loop.coordinate_frame[1][0]
		y1 = apical_loop.coordinate_frame[0][1]
		y2 = apical_loop.coordinate_frame[1][1]
		for loop_idx, loop_nt in enumerate(apical_loop.numbers):
			fraction_done_with_loop = (float(loop_idx)+0.5) / float(len( apical_loop.numbers ))
			x = self.interpolate( x1, x2, fraction_done_with_loop ) 
			y = self.interpolate( y1, y2, fraction_done_with_loop ) 
			self.draw_nt(apical_loop.sequence[loop_idx], (x, y))
		
	def draw_junction_loop_reasonably_with_respect_to_stems( self, junction_loop ):
		"""
		As long as we are just drawing this dumb left-to-right structure
		with no stacking, we should assume we are going from the right side of 
		stem1 to the left side of stems2
		"""
		x1 = junction_loop.stem1.coordinate_frame[0] + self.bp_offset_width
		x2 = junction_loop.stem2.coordinate_frame[0]
		for loop_idx, loop_nt in enumerate(junction_loop.numbers):
			# This fraction goes 0, 0.25, 0.5, 0.75
			# We need something more aggressive, esp. due to character alignment.
			fraction_done_with_loop = (float(loop_idx)+0.5) / float(len( junction_loop.numbers ))
			x = self.interpolate(x1, x2, fraction_done_with_loop )
			self.draw_nt(junction_loop.sequence[loop_idx], (x,-10))

