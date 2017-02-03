#import util
from secstruct_svgs.util import flatten
from secstruct_svgs.coordinate_frame import CoordinateFrame
from secstruct_svgs.base_pair import BasePair
from secstruct_svgs.nucleotide import Nucleotide

# Thoughts:
# 1. A Stem or Loop shouldn't know details of how it is rendered.
# 2. An object containing multiple Stems and Loops might know how to orient them relative
# to each other
# 3. A Stem and Loop can carry information directly pertinent to rendering. For example, 
# a Stem doesn't just have a sequence (with a length); it has a number of base pairs
# with types.

WATSON_CRICK = 1

class Stem(object):
    def __init__( self, sequence, stem_pattern=None ):
        """
        Initialize a stem with a given sequence. The length of the sequence 
        is used to size a variety of variables, such as its (relative/internal)
        coordinates.
        """
        
        sequence = "".join( [c for c in sequence if c != ',' and c != ' ' ] )
        
        if stem_pattern is not None:
            sequence = "".join( [ sequence[i-1] for i in sorted(flatten(stem_pattern) ) ] )

        self.sequence = sequence

        # AMW: Check later that this works for nontrivial examples
        # For example: will break on PK and intermolecular stems.
        #self.numbers  = sorted(flatten(stem_pattern) )
        # AMW: Actually, this may end up redoing the same organization as input
        # but I like being explicit
        numbers = sorted(flatten(stem_pattern) )
        self.numbers = [ ( numbers[i], numbers[len(numbers)-i-1] )  
                         for i in xrange(len(numbers)/2 ) ]

        # Trim common spacers from sequence
        # Fold sequence into BPs
        self.base_pairs = [ BasePair( 
            Nucleotide(sequence[i], numbers[i]), 
            Nucleotide(sequence[len(sequence)-i-1], numbers[len(sequence)-i-1]) )  
                            for i in xrange(len(sequence)/2 ) ]

        # Build up nucleotides.
        self.nucleotides = {}
        for bp in self.base_pairs:
            self.nucleotides[ bp.nt1.seqpos ] = bp.nt1
            self.nucleotides[ bp.nt2.seqpos ] = bp.nt2

        # Right now, set internal coordinates. First nt of every stem will seem 
        # like a ref right now, but the canvas can connect things up as needed.
        # second half of seq ref is first half.
        for seqpos, nt in self.nucleotides.iteritems():
            if seqpos == min( numbers ): continue # min(numbers) stays None, 0, 0
            
            if seqpos in [ bpnum[0] for bpnum in self.numbers ]:
                #self.numbers.index(seqpos) < len(self.numbers)/2:
                nt.ref_nt = self.nucleotides[ seqpos-1 ]
            else: #if self.numbers.index(seqpos) >= len(self.numbers)/2
                nt.ref_nt = nt.bp_partner
            nt.absolute_coords_need_updating = False
            nt.relative_coords_need_updating = True
            nt.update_relative_coords()
                            
        # All watson-crick for now.
        self.base_pair_types = [ WATSON_CRICK for _ in self.base_pairs ]

        # First rotation notion. Pair with 1 (DOWN) or -1 (UP)
        self.coordinate_frame = CoordinateFrame( [ 0, 0 ], 1 )
        print "in ctor, orientation ", self.coordinate_frame.orientation

    def update_absolute_coords(self):
        for nt in self.nucleotides.itervalues():
            if nt.ref_nt is None:
                continue
            print("Stem: updating abs coords for", nt.seqpos)
            print("NT coords pre: ", nt.x, nt.y, nt.ref_nt.seqpos)
            nt.update_absolute_coords()
            print("NT coords post:", nt.x, nt.y, nt.ref_nt.seqpos)

    def update_relative_coords(self):
        for nt in self.nucleotides.itervalues():
            nt.update_relative_coords()
