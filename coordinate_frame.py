

class CoordinateFrame:
	def __init__( self, pos, orientation, pos2=None ):
		"""
		Note pos2 is default None (stems don't have a terminus... yet?)
		"""
		self.position.x = pos[0]
		self.position.y = pos[1]
		self.orientation = orientation
		if pos2 is not None:
			self.position2.x = pos2[0]
			self.position2.y = pos2[1]
		
class Position:
	def __init__( self, x, y ):
		self.position.x = x
		self.position.y = y