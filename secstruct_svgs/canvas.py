from __future__ import print_function
import svgwrite
import math
from secstruct_svgs.util import distance, color, loop_interpolate, flatten#, interpolate
#from coordinate_frame import CoordinateFrame

class Canvas(object):
    """
    A Canvas is where all this stuff is added and spread out.
    """
    

    def __init__( self, dwg, width_per_stem=80, width_per_bp=30, height_per_bp=20 ):#, orig=None):
    #    if orig == None:
    #        self.ctor( dwg, width_per_stem, width_per_bp, height_per_bp )
    #    else:
    #        self.copyctor( orig )
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
    #    self.dwg = orig.dwg
    #    self.stem_offset_width = orig.width_per_stem
    #    self.bp_offset_width   = orig.width_per_bp
    #    self.bp_offset_height  = orig.height_per_bp
    #    self.font_height_offset = 3
    #    self.stems = [ Stem( orig.stem ) for 
    #    self.loops = []
    #    self.nucleotides = {}
        
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

        # Set initial position of stem-root
        x_offset = stem.coordinate_frame.position.x
        y = stem.coordinate_frame.position.y
        stem.nucleotides[ min(stem_numbers) ].x = x_offset
        stem.nucleotides[ min(stem_numbers) ].y = y
        # If there is an i-1, this isn't root
        if min(stem_numbers) - 1 in self.nucleotides.keys():
            stem.nucleotides[ min(stem_numbers) ].ref_nt = self.nucleotides[ min(stem_numbers) - 1 ]
            # Is this what we really want?
            bp.nt1.relative_coords_need_updating = True
            bp.nt1.absolute_coords_need_updating = False
            bp.nt1.update_relative_coords()
            
            bp.nt2.relative_coords_need_updating = True
            bp.nt2.absolute_coords_need_updating = False
            bp.nt2.update_relative_coords()
            
        for bp in stem.base_pairs: 
            if bp.nt1.seqpos != min(stem_numbers):# - 1:
                bp.nt1.dy = self.bp_offset_height * stem.coordinate_frame.orientation
                #if bp.nt1.seqpos - 1 in stem_numbers:
                #    bp.nt1.dy = bp.nt1.y - bp.nt1.ref_nt.y

            if bp.nt2.seqpos != min(stem_numbers):# - 1:
                bp.nt2.dx = self.bp_offset_width
            
            bp.nt1.absolute_coords_need_updating = True
            bp.nt1.relative_coords_need_updating = False
            bp.nt1.update_absolute_coords()

            bp.nt2.absolute_coords_need_updating = True
            bp.nt2.relative_coords_need_updating = False
            bp.nt2.update_absolute_coords()

        self.report_stem_coords()

    def report_stem_coords(self):
        print("Coordinates of stems:")
        for idx, stem in enumerate(self.stems):
            print("Stem", idx)
            for bp in stem.base_pairs:
                print("(%f, %f) -- (%f, %f)" % (bp.nt1.x, bp.nt1.y, bp.nt2.x, bp.nt2.y))

    def add_apical_loop(self, loop):
        """
        y offsets don't make a ton of sense here.
        """

        # Just trace path from coordinate frame start, plus y offset for
        # height, to other side of stem.
        loop.coordinate_frame.position.x  = loop.stem1.nucleotides[min(loop.numbers)-1].x
        loop.coordinate_frame.position.y  = loop.stem1.nucleotides[min(loop.numbers)-1].y
        # Offset for stem length
        #loop.coordinate_frame.position.y += loop.stem1.coordinate_frame.orientation * (len(loop.stem1.base_pairs) * self.bp_offset_height - self.bp_offset_height)
        # Extra for alignment (may change with font size?)
        loop.coordinate_frame.position.y += loop.stem1.coordinate_frame.orientation * 2
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

        for seqpos, nt in loop.nucleotides.items():
            if seqpos == min(loop.numbers): continue
            nt.ref_nt = loop.nucleotides[ seqpos - 1 ]

        for seqpos, loop_nt in loop.nucleotides.items():
            loop_idx = seqpos - min(loop.nucleotides.keys())
            fraction_done_with_loop = (float(loop_idx)+0.5) / float(len( loop.numbers ))
            # Note: loop_interpolate really ought to think about stem orientation
            # Maybe add 180 for flipped (-1) stems
            # Note that here we're for sure apical
            add_angle = 180 if loop.stem1.coordinate_frame.orientation == -1 else 0

            print("#### <A> Loop %s is %f done" % (loop.numbers, fraction_done_with_loop))
            loop_nt.x, loop_nt.y = loop_interpolate( x1,y1,x2,y2, total_loop_fraction, fraction_done_with_loop, add_angle=add_angle )
            print("#### <A> Placing", loop.numbers[loop_idx], loop_nt.name, " at ", loop_nt.x, loop_nt.y)
            loop_nt.absolute_coords_need_updating = False
            loop_nt.relative_coords_need_updating = True
            loop_nt.update_relative_coords()

    def add_junction_loop( self, loop ):
        """
        As long as we are just drawing this dumb left-to-right structure
        with no stacking, we should assume we are going from the right side of
        stem1 to the left side of stem2. Note that we will really have to think
        in terms of the positions of individual NTs eventually.

        Okay, now that we have matured a little and are looking at relative coordinates
        and more complicated topologies, junction loops can be really hard and nonlocal.
        We need to condense most of the problems, as it were, into one chainbreak that we then
        try to close. (Otherwise, we just have nutso coordinates and we have to fix 'all
        bond lengths' instead of one gap.)
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

        # AMW TEMP: we are right now moving towards a system where a loop already knows its position(s)
        # per nucleotide but we aren't there yet. So for now, note that loop_nt is suddenly a key to a dict
        # of Nucleotides -- not a number.
        
        # Don't even need the BL, just need dx dy
        #loop_bond_length = min(20, distance(loop.coordinate_frame.position,loop.coordinate_frame.position2)/len(loop.numbers))
        frac_dx = (x2 - x1)/len(loop.numbers)
        frac_dy = (y2 - y1)/len(loop.numbers)
        if frac_dx**2 + frac_dy**2 > 400:
            frac_dx *= 20 / math.sqrt(frac_dx**2 + frac_dy**2)
            frac_dy *= 20 / math.sqrt(frac_dx**2 + frac_dy**2)

        # Junction loops for coaxial stacks may suck for now... but maybe we can fix with MC now that
        # kinematics are more reliable.
        for loop_key, loop_nt in loop.nucleotides.items():
            # 'prior' nt guaranteed to exist.
            loop_nt.ref_nt = self.nucleotides[loop_key - 1]
            loop_nt.dx, loop_nt.dy = frac_dx, frac_dy
            loop_nt.relative_coords_need_updating = False
            loop_nt.absolute_coords_need_updating = True
            loop_nt.update_absolute_coords()
            print("#### <J> Placing", loop_key, loop_nt.name, " at ", loop_nt.x, loop_nt.y)

    def add_tail_loop(self, loop):
        """
        We can kind of make some choices here. When stem orientation is variable
        this will have to become a little more complicated, but at the moment we
        will be able to set the end of the coordinate frame based on length.
        """
        # Is the stem the beginning or end of the tail?
        if min(loop.numbers)-1 in flatten(loop.stem1.numbers): # beginning
            # Right now that means: extend to 'right'
            # TODO: nt width magic number.
            loop.coordinate_frame.position.x  = loop.stem1.nucleotides[min(loop.numbers)-1].x 
            loop.coordinate_frame.position.y  = loop.stem1.nucleotides[min(loop.numbers)-1].y
            loop.coordinate_frame.position2.x = loop.coordinate_frame.position.x  + 7 * len(loop.numbers)
            loop.coordinate_frame.position2.y = loop.coordinate_frame.position.y 
        if max(loop.numbers)+1 in flatten(loop.stem1.numbers): # beginning
            # Right now that means: extend to 'right'
            # TODO: nt width magic number.
            loop.coordinate_frame.position2.x = loop.stem1.nucleotides[max(loop.numbers)+1].x
            loop.coordinate_frame.position2.y = loop.stem1.nucleotides[max(loop.numbers)+1].y
            loop.coordinate_frame.position.x  = loop.coordinate_frame.position2.x  - 7 * len(loop.numbers)
            loop.coordinate_frame.position.y  = loop.coordinate_frame.position.y 

    def add_loop( self, loop ):
        self.loops.append(loop)
        for seqpos, loop_nt in loop.nucleotides.items(): self.nucleotides[ seqpos ] = loop_nt
        if loop.apical:
            self.add_apical_loop( loop )
        elif loop.tail:
            self.add_tail_loop(loop)
        else:
            self.add_junction_loop( loop )

    def set_stems_coaxial( self, idx1, idx2 ):
        """
        JUST add to list of coaxialities. Previous use of this function had to do with coordinate frames
        and later NT positions... this is now the job of initialize_stem_relative_positions().
        """

        self.list_of_coaxialities.append( [idx1, idx2] )

    def initialize_stem_relative_positions( self ):
        # Uh oh! We had committed to this "single coordinate frame per stem" idea. But now we really want
        # to interact with nucleotides directly!

        # All coaxial pairs: first nt of idx2 now has a root of final from idx1
        for idx1, idx2 in self.list_of_coaxialities:
            st2 = self.stems[idx2]
            st1 = self.stems[idx1]
            st2.nucleotides[min(flatten(st2.numbers))].ref_nt = st1.nucleotides[max(flatten(st1.numbers))]
            # dy is 1 plus self.bp_offset_height
            st2.nucleotides[min(flatten(st2.numbers))].dy = st2.coordinate_frame.orientation * (self.bp_offset_height + 1)
            st2.nucleotides[min(flatten(st2.numbers))].relative_coords_need_updating = False
            st2.nucleotides[min(flatten(st2.numbers))].absolute_coords_need_updating = True

            cf1 = st1.coordinate_frame
            cf2 = st2.coordinate_frame
        
            # idx2 is opposite of idx1 in orientation. Not necessarily going to be true in
            # all kinds of coaxiality, but in this first, simplest kind...
            cf2.orientation = -1 * cf1.orientation

            # may be irrelevant now... are we using this?
            cf2.position.y = cf1.position.y + cf2.orientation * ( self.bp_offset_height + 1 )
            cf2.position.x = cf1.position.x + self.bp_offset_width # coaxiality, keeping in mind the flipped 'sense'

            # since orientation is now opposite we have to flip all our dx and dy signs.
            #x_offset = cf2.position.x
            #y = cf2.position.y
            for bp in st2.base_pairs: 
                bp.nt1.dx = -bp.nt1.dx
                bp.nt1.dy = -bp.nt1.dy
                bp.nt1.absolute_coords_need_updating = True
                bp.nt1.relative_coords_need_updating = False
                
                bp.nt2.dx = -bp.nt2.dx
                bp.nt2.dy = -bp.nt2.dy
                bp.nt2.absolute_coords_need_updating = True
                bp.nt2.relative_coords_need_updating = False
                # need to compute some offet? ...
            
            st1.update_absolute_coords()
            st2.update_absolute_coords()
        
        # Now, set all additional stems based on black magic basically
        for idx, stem in enumerate(self.stems):
            if idx in flatten(self.list_of_coaxialities): 
                continue
            nt_base_idx = min(flatten(stem.numbers))
            if idx == 0:
                # These are actually the defaults, but let's be defensive...
                stem.nucleotides[nt_base_idx].ref_nt = None
                stem.nucleotides[nt_base_idx].dx = 0
                stem.nucleotides[nt_base_idx].dy = 0
                
                stem.nucleotides[nt_base_idx].absolute_coords_need_updating = True
                stem.nucleotides[nt_base_idx].relative_coords_need_updating = False

                continue
            #print("Stem %d has nt base %d and its ref_nt will be set to number %d" % (idx, nt_base_idx, min(flatten(self.stems[idx-1].numbers))))
            stem.nucleotides[nt_base_idx].ref_nt = self.stems[idx-1].nucleotides[min(flatten(self.stems[idx-1].numbers))]
            stem.nucleotides[nt_base_idx].dx = self.stem_offset_width
            stem.nucleotides[nt_base_idx].dy = 0
            
            stem.nucleotides[nt_base_idx].absolute_coords_need_updating = True
            stem.nucleotides[nt_base_idx].relative_coords_need_updating = False
           
            for bp in stem.base_pairs:
                #print("BP: ", bp.nt1.dx, bp.nt1.dy, bp.nt2.dx, bp.nt2.dy)
                bp.nt1.absolute_coords_need_updating = True
                bp.nt1.relative_coords_need_updating = False
                bp.nt2.absolute_coords_need_updating = True
                bp.nt2.relative_coords_need_updating = False

            #stem.nucleotides[min(flatten(stem.numbers))].x = idx * self.stem_offset_width
            #stem.nucleotides[min(flatten(stem.numbers))].y = 0

            stem.update_absolute_coords()

            #print("Just now gonna report after this stem", idx)
            #self.report_stem_coords()

        # Oh! Update any dependent loops!
        # You know what, I bet the coordinates aren't even really set yet, and the issue is that they're 
        # gonna be set based on the coordinate frames, not the nt absolute coordinates.
        for loop in self.loops:
            print("Updating loop containing", loop.nucleotides)
            # Just update them all
            for _, nt in loop.nucleotides.items():
                nt.absolute_coordinates_need_updating = True
                nt.update_absolute_coords()

def svg_render(canvas):
    def draw_text(canvas, text, loc, color):
        canvas.dwg.add(canvas.dwg.text(text, insert=loc, fill=color ) )

    def draw_line(canvas, loc1, loc2, stroke=svgwrite.rgb(10, 10, 16, '%')):
        canvas.dwg.add(canvas.dwg.line(loc1, loc2, stroke=stroke))

    def draw_nt(canvas, nt):
        canvas.dwg.add(canvas.dwg.text(nt.name, insert=(nt.x, nt.y), fill=color(nt.name) ) )

        # numbers if multiple of 5
        # TODO: adjust size
        if nt.seqpos % 5 == 0: draw_text(canvas, nt.seqpos, (nt.x-15, nt.y+5), color(nt.name))

    def draw_bp(canvas, bp1, bp2):
        """
        Needs some kind of line centering relative to the height of the characters.
        """
        draw_nt(canvas, bp1)
        draw_nt(canvas, bp2)
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
        #                (bp1.x + 25, bp1.y - self.font_height_offset) )
        draw_line(canvas, beg, end)

    def draw_stem(canvas, stem):
        # draw basepairs
        for bp in stem.base_pairs: 
            draw_bp(canvas, bp.nt1, bp.nt2)

    def draw_apical_loop(canvas, apical_loop):
        for loop_nt in apical_loop.nucleotides.values():
            draw_nt(canvas, loop_nt)

    def draw_junction_loop(canvas, junction_loop):
        for loop_nt in junction_loop.nucleotides.values():
            draw_nt(canvas, loop_nt)

    def draw_sequence_line(canvas, nt1, nt2):
        """
        Draw a line between consecutive nucleotides. Thin and grey. 
        Figure out where the subjective centers of the nts are ( depends on font but likely about 
        5-6 px up and right of (x,y) ) and draw about 80% of that vector.
        """

        center1 = [nt1.x + canvas.font_height_offset, nt1.y - canvas.font_height_offset]
        center2 = [nt2.x + canvas.font_height_offset, nt2.y - canvas.font_height_offset]

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

        draw_line(canvas, beg, end, 'gray')

    for stem in canvas.stems: draw_stem(canvas, stem)
    for apical in [loop for loop in canvas.loops if loop.apical]:
        draw_apical_loop(canvas, apical)
    for junction in [loop for loop in canvas.loops if not loop.apical]:
        draw_junction_loop(canvas, junction)
    for seqpos in canvas.nucleotides.keys():
        if seqpos + 1 in canvas.nucleotides.keys():
            draw_sequence_line(canvas, canvas.nucleotides[seqpos], canvas.nucleotides[seqpos + 1])
    
    canvas.dwg.save()
