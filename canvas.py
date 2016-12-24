import svgwrite
from util import color

class Canvas:
	def __init__( self, dwg, width_per_stem=80, width_per_bp=30, height_per_bp=20):
		self.dwg = dwg
		self.stem_offset_width = width_per_stem
		self.bp_offset_width   = width_per_bp
		self.bp_offset_height  = height_per_bp
		self.stems = []

	def add_stem( self, stem ):
		# draw basepairs
		stem_idx = len(self.stems)
		print  stem_idx, stem.sequence
		x_offset = stem_idx * self.stem_offset_width
		for bp_idx, [bp1, bp2] in enumerate(stem.base_pairs): 
			# need to compute some offet? ...
			y_offset = bp_idx * self.bp_offset_height
			self.dwg.add(self.dwg.line((x_offset + 10, y_offset), (x_offset + 20, y_offset), 
					stroke=svgwrite.rgb(10, 10, 16, '%')))
			self.dwg.add(self.dwg.text(bp1, insert=(x_offset+0, y_offset), fill=color(bp1)))#'blue'))
			self.dwg.add(self.dwg.text(bp2, insert=(x_offset+self.bp_offset_width, y_offset), fill=color(bp2)))#'yellow'))
	
			# numbers if multiple of 5
			# TODO: adjust size
			num1, num2 = stem.numbers[bp_idx]
			if num1 % 5 == 0: self.dwg.add(self.dwg.text(num1, insert=(x_offset+(-10), y_offset+5), fill=color(bp1)))#'blue'))
			if num2 % 5 == 0: self.dwg.add(self.dwg.text(num2, insert=(x_offset+30+(10), y_offset+5), fill=color(bp2)))#'blue'))

		self.stems.append( stem )
		return self.dwg
		
	def draw_apical_loop_reasonably_with_respect_to_stem( self, apical_loop, stem_idx, stem, seq ):
		"""
		In a better world, this stem would have attached coordinates.
		Instead I'm just passing its stem_idx and computing the offset.
		"""
		for loop_idx, loop_nt in enumerate(apical_loop):
			x_offset = stem_idx * self.stem_offset_width + loop_idx * 8
			# This would come out of stem length
			y_offset = len(stem.base_pairs) * self.bp_offset_height + 10
			self.dwg.add(self.dwg.text(seq[loop_nt-1], insert=(x_offset, y_offset), fill=color(seq[loop_nt-1])))#'blue'))
		return self.dwg
		
	def draw_junction_loop_reasonably_with_respect_to_stems( self, junction_loop, stem1_idx, stem1, stem2_idx, stem2, seq ):
		#y1 = 0 - 10
		#y2 = 0 - 10
		x1 = stem1_idx * self.stem_offset_width + self.bp_offset_width
		x2 = stem2_idx * self.stem_offset_width
		for loop_idx, loop_nt in enumerate(junction_loop):
			print stem1_idx, stem2_idx, self.stem_offset_width, self.bp_offset_width
			fraction_done_with_loop = float(loop_idx) / float(len( junction_loop ))
			print fraction_done_with_loop
			print (x1+(float(x2-x1)* fraction_done_with_loop ) )
			self.dwg.add(self.dwg.text(seq[loop_nt-1], insert=((x1+(float(x2-x1)* fraction_done_with_loop ) ), -10 ), fill=color(seq[loop_nt-1])))#'blue'))
		return self.dwg

