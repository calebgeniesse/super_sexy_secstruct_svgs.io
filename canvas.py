import svgwrite
from util import color

class Canvas:
	def __init__( self, dwg, width_per_stem=80, width_per_bp=30, height_per_bp=20):
		self.dwg = dwg
		self.stem_offset_width = width_per_stem
		self.bp_offset_width   = width_per_bp
		self.bp_offset_height  = height_per_bp
		self.stems = []
	
	def add_text( self, text, loc, color ):
		self.dwg.add(self.dwg.text(text, insert=loc, fill=color ) )
	
	def add_line( self, loc1, loc2 ):
		""" 
		Only supports black for now.
		"""
		self.dwg.add(self.dwg.line(loc1, loc2, stroke=svgwrite.rgb(10, 10, 16, '%')))

	def add_nt( self, text, loc ):
		self.dwg.add(self.dwg.text(text, insert=loc, fill=color(text) ) )
	
	def add_bp( self, bp1, bp2, loc ):
		self.add_nt( bp1, loc )
		self.add_nt( bp2, (loc[0]+self.bp_offset_width, loc[1]) )
		self.add_line( (loc[0] + 10, loc[1]), (loc[0]+ 20, loc[1]) )
	
	def add_stem( self, stem ):
		# draw basepairs
		stem_idx = len(self.stems)
		print  stem_idx, stem.sequence
		x_offset = stem_idx * self.stem_offset_width
		
		# update coordinate frame of stem
		# Let's say the x units of frame are units of stem number for now. brilliant, right?
		# wait no -- we can calculate using canvas members like stem_offset_width and
		# bp_offset_width. Brilliant!
		#stem.coordinate_frame = ( len(self.stems) * self.stem_offset_width + (len(self.stems)-1) * self.bp_offset_width, 0 )
		stem.coordinate_frame = ( len(self.stems) * self.stem_offset_width, 0 )

		for bp_idx, [bp1, bp2] in enumerate(stem.base_pairs): 
			# need to compute some offet? ...
			y_offset = bp_idx * self.bp_offset_height
			self.add_bp( bp1, bp2, (x_offset, y_offset) )

			# numbers if multiple of 5
			# TODO: adjust size
			num1, num2 = stem.numbers[bp_idx]
			if num1 % 5 == 0: self.add_text( num1, (x_offset-10, y_offset+5), color(bp1) )
			if num2 % 5 == 0: self.add_text( num2, (x_offset+self.bp_offset_width+10, y_offset+5), color(bp2) )

		self.stems.append( stem )
		return self.dwg
		
	def draw_apical_loop_reasonably_with_respect_to_stem( self, apical_loop ):
		for loop_idx, loop_nt in enumerate(apical_loop.numbers):
			x_offset = apical_loop.stem1.coordinate_frame[0] + loop_idx * 8
			# This would come out of stem length
			y_offset = len(apical_loop.stem1.base_pairs) * self.bp_offset_height + 10
			self.dwg.add(self.dwg.text(apical_loop.sequence[loop_idx], insert=(x_offset, y_offset), fill=color(apical_loop.sequence[loop_idx])))#'blue'))
		return self.dwg
		
	def draw_junction_loop_reasonably_with_respect_to_stems( self, junction_loop ):
		#y1 = 0 - 10
		#y2 = 0 - 10
		x1 = junction_loop.stem1.coordinate_frame[0]
		x2 = junction_loop.stem2.coordinate_frame[0]
		for loop_idx, loop_nt in enumerate(junction_loop.numbers):
			fraction_done_with_loop = float(loop_idx) / float(len( junction_loop.numbers ))
			print fraction_done_with_loop
			print (x1+(float(x2-x1)* fraction_done_with_loop ) )
			self.dwg.add(self.dwg.text(junction_loop.sequence[loop_idx], insert=((x1+(float(x2-x1)* fraction_done_with_loop ) ), -10 ), fill=color(junction_loop.sequence[loop_idx])))#'blue'))
		return self.dwg

