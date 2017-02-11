"""
I don't believe this file ought to require a docstring. TODO.
"""

from __future__ import print_function

import sys
import os.path
import math
import copy

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

#import secstruct_svgs
from secstruct_svgs.stem import Stem
from secstruct_svgs.loop import Loop
from secstruct_svgs.canvas import Canvas, svg_render
from secstruct_svgs.util import get_stems, flatten, consecutive_segments, seq_for, score, color #, distance, distance_squared, angle, closer_than, flat_harmonic_penalty, harmonic_penalty
import numpy as np

import svgwrite

# test
#from util import loop_interpolate
#x1, y1 = [-1,0]
#x2, y2 = [ 1,0]
#print loop_interpolate( x1,y1,x2,y2, 0.75, 1.0/9.0 )
#exit()

# General theory:
#
# We are interested in the relative arrangement of a set of stems.
#
# This is the most important problem to solve -- the placement of loop nucleotides
# is both less geometrically constrained by "traditional publication-quality aesthetics"
# and easier to resolve (in the sense of a loop-modeling problem, really).
#
# This is why we're first trying to *develop an API* against which to solve this problem.
#
# Stems have a position and orientation relative to the 'root stem' (at the moment, just
# stem 0; at the moment, the position and orientation *cannot* vary). So we can compute
# 'stem internal coordinates' -- trivial -- then transform them.

def get_loops(sequence, true_stems):
    """
    Given a sequence and the numbers that have been committed to stems,
    this function extracts a set of loops -- first in terms of numbers
    ('consecutive segments'), then in terms of Loop objects with Stem
    attachments.
    TODO: move to 'sequence utils' file.
    """
    loop_nt = [i + 1 for i in range(len(sequence)) if i+1 not in flatten([stem.numbers for stem in true_stems])]
    #print "loop_nt", loop_nt
    loops = consecutive_segments( loop_nt )
    print(loops)
    real_loops = []
    for loop in loops:
        adjacent_stems = [ stem for stem in true_stems if min(loop)-1 in flatten(stem.numbers) or max(loop)+1 in flatten(stem.numbers) ]
        apical_stems = [ stem for stem in true_stems if min(loop)-1 in flatten(stem.numbers) and max(loop)+1 in flatten(stem.numbers) ]
        if len(adjacent_stems) == 0:
            print("loop", loop, "has no adjacent stems!")
            continue
        if len(apical_stems) == 1:
            print("apical loop", loop)
            print("conn to", apical_stems[0].numbers)
            real_loops.append( Loop( seq_for( loop, sequence ), loop, apical_stems[0] ) )
        elif len(adjacent_stems) == 1: # tail
            print("tail loop", loop)
            print("conn to", adjacent_stems[0].numbers)
            real_loops.append( Loop( seq_for( loop, sequence ), loop, adjacent_stems[0], False, True ) )
        else:
            print("junct loop", loop)
            print("conn to", adjacent_stems[0].numbers)
            print("conn to", adjacent_stems[1].numbers)
            real_loops.append( Loop( seq_for( loop, sequence ), loop, adjacent_stems[0], False, False, adjacent_stems[1] ) )
    return real_loops


def apicals( stems, loops ):
    """
    Finds all the apical loops in the given set of loops. All apical loops stretch
    from i+1 to j-1 for some base pair in a stem.
    """

    # AMW: for draw_apical_loop_reasonably_with_respect_to_stem, we actually
    # want + need the stem_idx for now.
    apicals = []
    for loop in loops:
        print("Apical?", loop)
        for stem_idx, stem in enumerate(stems):
            print("\t", stem.numbers)
            if (min(loop)-1, max(loop)+1) in stem.numbers:
                apicals.append([(stem_idx, stem), loop])
                break

    return apicals

def junctions( stems, loops ):
    """
    Junctions aren't just 'non-apical'; they're also 'non-tail'
    So instead of just taking the non apical loops, we instead have to ensure that our
    candidate goes from max(stem1)+1 to min(stem2)-1
    """
    junctions = []
    for loop in loops:
        for stem1_idx, stem1 in enumerate(stems):
            for stem2_idx, stem2 in enumerate(stems):
                if stem1_idx == stem2_idx: continue
                if stem2_idx != stem1_idx + 1: continue
                if max(flatten(stem1.numbers))+1 == min(loop) and min(flatten(stem2.numbers))-1 == max(loop):
                    junctions.append([(stem1_idx, stem1), (stem2_idx, stem2), loop])

    print(junctions)
    return junctions


def perturb_loops( canvas ):
    """
    Random perturbation to each loop nt, propagating
    to adjacent nts
    """

    # Thought: ignore non-junction, as apicals are
    # local and likely nearly ideal to start with.

    new_canvas = copy.deepcopy(canvas)#Canvas(canvas)

    for loop in new_canvas.loops:

        if loop.apical: continue

        # Achieve OK stem kinematics by perturbing one nt in a stem and marking the whole thing as moved.
        moved_this_turn = { seqpos: False for seqpos in loop.nucleotides.keys() }

        for seqpos, nt in loop.nucleotides.items():
            if moved_this_turn[ seqpos ]: continue

            moved_this_turn[ seqpos ] = True
            dx = np.random.normal(0,3)
            dy = np.random.normal(0,3)

            #print "Moving nt %d from ( %f, %f )... " % ( seqpos, new_canvas.nucleotides[seqpos].x, new_canvas.nucleotides[seqpos].y ),
            nt.x += dx
            nt.y += dy
            nt.relative_coordinates_need_updating = True
            nt.update_relative_coords()
            #print "... to ( %f, %f )!" % ( new_canvas.nucleotides[seqpos].x, new_canvas.nucleotides[seqpos].y )

            # Each of the other nts in the loop do a proportion of the same motion
            # depending on how far away
            for seqpos2, nt2 in loop.nucleotides.items():
                if seqpos2 == seqpos: continue
                prefactor =  0.8 / abs( seqpos - seqpos2 )
                nt2.x += prefactor * dx
                nt2.y += prefactor * dy
                nt2.relative_coordinates_need_updating = True
                nt2.update_relative_coords()
            

    return new_canvas


def perturb( canvas ):
    """
    Random perturbation to each NT. Gaussian?
    Should ultimately have realistic kinematics... stems should move together, or at least BPs should.
    AMW: Right now we have partial kinematics, but really we should just update coords to propagate changes!
    """


    # Achieve OK stem kinematics by perturbing one nt in a stem and marking the whole thing as moved.
    moved_this_turn = { seqpos: False for seqpos in canvas.nucleotides.keys() }

    new_canvas = copy.deepcopy(canvas)#Canvas(canvas)
    
    # Perturb stems directly
    for idx, stem in enumerate(new_canvas.stems):

        if np.random.uniform() < 0.1:
            # Try a rotation
            for bp_idx, bp in enumerate(stem.base_pairs):
                complement_idx = len(stem.base_pairs) - bp_idx - 1
                comp_bp = stem.base_pairs[complement_idx]
                bp.nt1.x = comp_bp.nt2.x
                bp.nt1.y = comp_bp.nt2.y
                bp.nt2.x = comp_bp.nt1.x
                bp.nt2.y = comp_bp.nt1.y
                for nt in (bp.nt1, bp.nt2, comp_bp.nt1, comp_bp.nt2):
                    nt.relative_coordinates_need_updating = True
                    nt.update_relative_coords()

            continue

        # Move 'root'
        del_x = np.random.normal(0,3)
        del_y = np.random.normal(0,3)
        nt = stem.nucleotides[min(flatten(stem.numbers))]
        nt.x += del_x
        nt.y += del_y
        for seqpos, nt in stem.nucleotides.items():
            moved_this_turn[ seqpos ] = True
            if seqpos != min(flatten(stem.numbers)):
                # this one was just updated
                nt.absolute_coordinates_need_updating = True
                nt.update_absolute_coords()

        # update dependent loops
        for loop in new_canvas.loops:
            if loop.stem1 is stem:
                # we fold from stem
                for seqpos, nt in loop.nucleotides.items():
                    nt.absolute_coordinates_need_updating = True
                    nt.update_absolute_coords()
      
    for loop in new_canvas.loops:
        if loop.apical: continue
        for seqpos, nt in loop.nucleotides.items():
            if moved_this_turn[ seqpos ]: continue

            moved_this_turn[ seqpos ] = True
            dx = np.random.normal(0,3)
            dy = np.random.normal(0,3)

            # This was a perturbation to a loop. Partial-propagate.
            # still skip apicals
            nt.x += dx
            nt.y += dy
            for seqpos2, nt2 in loop.nucleotides.items():
                if seqpos2 == seqpos: continue
                moved_this_turn[ seqpos2 ] = True
                prefactor =  0.8 / abs( seqpos - seqpos2 )
                nt2.x += prefactor * dx
                nt2.y += prefactor * dy
        
    return new_canvas

def metropolis( new, old, temp ): 
    return np.random.uniform() < math.exp( ( old - new ) / temp )

def prepack(canvas):
    """
    Utterly broken, but this gives you the idea.
    """
    for coaxiality in canvas.list_of_coaxialities:
        old_score = score(canvas)
        reversed_canvas = copy.deepcopy(canvas)
        # there is a better way to do this
        reversed_canvas.list_of_coaxialities[reversed_canvas.list_of_coaxialities.index(coaxiality)] = [coaxiality[1], coaxiality[0]]
        reversed_canvas.initialize_stem_relative_positions()
        for stem in reversed_canvas.stems:
            stem.update_absolute_coords()
        new_score = score(reversed_canvas)
        print("Attempting coaxiality reversal, old score", old_score, "new score", new_score)
        if new_score < old_score:
            canvas = reversed_canvas
    return canvas

def mc_loops( canvas ):
    cycles = 1000
    for x in range(cycles):
        old_score = score(canvas)
        new_canvas = perturb_loops(canvas)
        new_score = score(new_canvas)
        if new_score < old_score or metropolis(new_score, old_score, temp=1):
            # accept
            print("Loop phase cycle %d: Accepted perturbation from %f to %f." % (x, old_score, new_score))
            canvas = new_canvas

        else:
            print("Loop phase cycle %d: Rejected perturbation from %f to %f." % (x, old_score, new_score))
    return canvas

def mc( canvas ):
    """
    This is an in-progress framework for doing MCMC simulations on a canvas.
    The idea is that you score NT configurations, update, etc.
    """
    cycles = 1000
    for x in range(cycles):
        old_score = score(canvas)
        new_canvas = perturb(canvas)
        new_score = score(new_canvas)
        if new_score < old_score or metropolis( new_score, old_score, temp=1 ):
            # accept
            print("All-phase cycle %d: Accepted perturbation from %f to %f." % (x, old_score, new_score))
            canvas = new_canvas

        else:
            print("All-phase cycle %d: Rejected perturbation from %f to %f." % (x, old_score, new_score))
    return canvas



def tk_render(canvas, tk_canvas):
    def draw_text(tk_canvas, text, loc, color):
        tk_canvas.create_text(text, position=loc, fill=color)

    def draw_line(tk_canvas, loc1, loc2, stroke=svgwrite.rgb(10, 10, 16, '%')):
        tk_canvas.create_line(loc1[0], loc1[1], loc2[0], loc2[1])#, stroke=stroke)

    def draw_nt(tk_canvas, nt):
        draw_text(tk_canvas, nt.name, (nt.x, nt.y), color(nt.name))

        # numbers if multiple of 5
        # TODO: adjust size
        if nt.seqpos % 5 == 0: draw_text(tk_canvas, nt.seqpos, (nt.x-15, nt.y+5), color(nt.name))

    def draw_bp(tk_canvas, bp1, bp2):
        """
        Needs some kind of line centering relative to the height of the characters.
        """
        draw_nt(tk_canvas, bp1)
        draw_nt(tk_canvas, bp2)
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
        draw_line(tk_canvas, beg, end)

    def draw_stem(tk_canvas, stem):
        # draw basepairs
        for bp in stem.base_pairs: 
            draw_bp(tk_canvas, bp.nt1, bp.nt2)

    def draw_apical_loop(tk_canvas, apical_loop):
        for loop_nt in apical_loop.nucleotides.values():
            draw_nt(tk_canvas, loop_nt)

    def draw_junction_loop(tk_canvas, junction_loop):
        for loop_nt in junction_loop.nucleotides.values():
            draw_nt(tk_canvas, loop_nt)

    def draw_sequence_line(tk_canvas, nt1, nt2):
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

        draw_line(tk_canvas, beg, end, 'gray')

    for stem in canvas.stems: draw_stem(canvas, stem)
    for apical in [loop for loop in canvas.loops if loop.apical]:
        draw_apical_loop(tk_canvas, apical)
    for junction in [loop for loop in canvas.loops if not loop.apical]:
        draw_junction_loop(tk_canvas, junction)
    for seqpos in canvas.nucleotides.keys():
        if seqpos + 1 in canvas.nucleotides.keys():
            draw_sequence_line(tk_canvas, canvas.nucleotides[seqpos], canvas.nucleotides[seqpos + 1])
    
    canvas.dwg.save()


def main():


    fn = 's_s_ss.svg'

    seq = 'ccccuuucccccaaaaggggguuuccccccaaaagggggguuucccccccaaaaggggggguuuugggg'
    ss  = '((((...(((((....)))))...((((((....))))))...(((((((....)))))))....))))'

    # Different loop lengths
    #seq = 'ccaaccgcaagguuggaucccauguuaaaaacgcaug'
    #ss  = '((((((....))))))....((((.........))))'

    # Different stem lengths
    #seq = 'ccaaccgcaagguuggaucccauguucgcaug'
    #ss  = '((((((....))))))....((((....))))'
    #seq = 'ccccgcaaggggaucccauguucgcaug'
    #ss  = '((((....))))....((((....))))'
    # Eventually: support chainbreaks

    # This produces old-style numerical stems. I want to turn them into
    # sequence stems. Use the sequence + numerical stem ctor.
    stems = [Stem(seq, stem) for stem in get_stems(ss, seq)]
    # Find loops by process of elimination
    loops = get_loops(seq, stems)

    dwg = svgwrite.Drawing(fn)
    canvas = Canvas(dwg)
    
    for stem in stems: canvas.add_stem(stem)
    
    # Inform the canvas that stem 1 and stem 2 are coaxial.
    # This doesn't make that great semantic sense, but the canvas
    # is the only stem-collection type object so it is aware of this
    # sort of thing.

    # Also, how is the API user really going to know?
    # You'd imagine they'd actually want to set it based on
    # labels containing residues or something...
    canvas.set_stems_coaxial(0, 1)

    canvas.initialize_stem_relative_positions()


    # We have to do this before adding any loops because (initial) positions
    # are calculated after addition time.
    for loop in loops: canvas.add_loop(loop)

    def void_mc(canvas): canvas = mc(canvas)
    def void_mc_loops(canvas): canvas = mc_loops(canvas)

    import tkinter
    from tkinter import Button
    from tkinter import Canvas as tkCanvas
    window = tkinter.Tk()
    # Big moves that don't require subtlety. For example: what's the best-scoring stem coaxiality?
    # (given a particular list of coaxialities)
    do_prepack = Button(window, text="Do prepack", command=lambda: prepack(canvas))
    do_prepack.pack()
    do_render = Button(window, text="Render to file", command=lambda: svg_render(canvas))
    do_render.pack()
    do_mc_loops = Button(window, text="Monte Carlo loops", command=lambda: void_mc_loops(canvas))
    do_mc_loops.pack()
    do_mc = Button(window, text="Monte Carlo overall", command=lambda: void_mc(canvas))
    do_mc.pack()
    tk_canvas = tkCanvas(window, width=500, height=500)
    tk_render(canvas, tk_canvas)
    tk_canvas.create_line(0, 0, 500, 500)
    tk_canvas.pack()

    window.mainloop()

    # Maybe a relaxation phase where we only move loops -- and move
    # them in some sort of sensible synchrony -- is in order.
    #canvas = mc_loops(canvas)
    #canvas = mc(canvas)
    
    #svg_render(canvas)

if __name__=="__main__":
    main()
    