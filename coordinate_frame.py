

class CoordinateFrame:
	def __init__( self, pos, orientation, pos2=None ):
		"""
		Note pos2 is default None (stems don't have a terminus... yet?)
		"""
		self.position = Position( pos[0], pos[1] )
		self.orientation = orientation
		if pos2 is not None:
			self.position2 = Position( pos2[0], pos2[1] )
		
class Position:
	def __init__( self, x, y ):
		self.position.x = x
		self.position.y = y