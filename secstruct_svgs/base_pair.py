#from nucleotide import Nucleotide
class BasePair(object):
    def __init__( self, nt1, nt2 ):
        self.nt1 = nt1
        self.nt2 = nt2
        self.nt1.bp_partner = nt2
        self.nt2.bp_partner = nt1