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
from secstruct_svgs.canvas import Canvas
from secstruct_svgs.util import get_stems, flatten, consecutive_segments, seq_for, score #, distance, distance_squared, angle, closer_than, flat_harmonic_penalty, harmonic_penalty
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
    loop_nt = [i + 1 for i in xrange(len(sequence)) if i+1 not in flatten([stem.numbers for stem in true_stems])]
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

        for seqpos, nt in loop.nucleotides.iteritems():
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
            for seqpos2, nt2 in loop.nucleotides.iteritems():
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
        for seqpos, nt in stem.nucleotides.iteritems():
            moved_this_turn[ seqpos ] = True
            if seqpos != min(flatten(stem.numbers)):
                # this one was just updated
                nt.absolute_coordinates_need_updating = True
                nt.update_absolute_coords()

        # update dependent loops
        for loop in new_canvas.loops:
            if loop.stem1 is stem:
                # we fold from stem
                for seqpos, nt in loop.nucleotides.iteritems():
                    nt.absolute_coordinates_need_updating = True
                    nt.update_absolute_coords()
      
    for loop in new_canvas.loops:
        if loop.apical: continue
        for seqpos, nt in loop.nucleotides.iteritems():
            if moved_this_turn[ seqpos ]: continue

            moved_this_turn[ seqpos ] = True
            dx = np.random.normal(0,3)
            dy = np.random.normal(0,3)

            # This was a perturbation to a loop. Partial-propagate.
            # still skip apicals
            nt.x += dx
            nt.y += dy
            for seqpos2, nt2 in loop.nucleotides.iteritems():
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
    for x in xrange(cycles):
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
    for x in xrange(cycles):
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

if __name__=="__main__":

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

    # Big moves that don't require subtlety. For example: what's the best-scoring stem coaxiality?
    # (given a particular list of coaxialities)
    canvas = prepack(canvas)
    
    # Maybe a relaxation phase where we only move loops -- and move
    # them in some sort of sensible synchrony -- is in order.
    #canvas = mc_loops(canvas)
    #canvas = mc(canvas)
    
    canvas.render()
