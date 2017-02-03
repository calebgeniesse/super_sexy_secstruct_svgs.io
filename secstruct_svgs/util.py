import math

def mod(vec): return math.sqrt( vec[0] ** 2 + vec[1] ** 2 )

def mod_squared(vec): return vec[0] ** 2 + vec[1] ** 2

def dot(v1, v2): return v1[0] * v2[0] + v1[1] * v2[1]

def vector(nt1, nt2): return [nt2.x - nt1.x, nt2.y - nt1.y]

def distance(nt1, nt2): return mod(vector(nt1, nt2))

def distance_squared(nt1, nt2): return mod_squared(vector(nt1, nt2))

def vectors_close(v1, v2): 
    return v1[0] - v2[0] < 0.0000001 and v1[1] - v2[1] < 0.0000001


def angle( nt1, nt2, nt3 ):
    """ Math looks right on this one but I am getting the supplement. So, fixing..."""
    if vector(nt1, nt2) == [0,0]:
        print("nt1", nt1.seqpos, " at ", nt1.x, nt1.y, " is at the same position as nt2", nt2.seqpos)
    if vector(nt2, nt3) == [0,0]:
        print("nt2", nt2.seqpos, " at ", nt2.x, nt2.y, " is at the same position as nt3", nt3.seqpos)
    #print(vector(nt1, nt2), vector(nt2, nt3))
    if vectors_close(vector(nt1, nt2), vector(nt2, nt3)):
        # These vectors are identical and that is messing with the ability to call two things parallel?
        return 180.0
    return 180.0 - math.degrees(math.acos(dot(vector(nt1, nt2), vector(nt2, nt3)) / (mod(vector(nt1, nt2)) * mod(vector(nt2, nt3)))))

def harmonic_penalty( dist, ideal, spring_constant ):
    return spring_constant * ( dist - ideal ) ** 2

def flat_harmonic_penalty( dist, min_d, max_d, spring_constant ):
    if dist < min_d:
        return spring_constant * ( dist - min_d ) ** 2
    elif dist > max_d:
        return spring_constant * ( dist - max_d ) ** 2
    else:
        return 0

def closer_than( val1, val2, diff ): return  val2 - val1 < diff or val1 - val2 < diff

def seq_for( loop, seq ):
    return "".join([ seq[idx-1] for idx in sorted(loop) ])

def color(nt):
    if nt == 'a': return 'yellow'
    if nt == 'c': return 'green'
    if nt == 'g': return 'red'
    if nt == 'u': return 'blue'
    
    print "WARNING: unrecognized nucleotide."

def interpolate( v1, v2, frac ): return v1 + float(v2-v1)*frac

# loops are consecutive sequences of non-stem residues
def flatten( obj ):
    arr = []
    try:
        for i in obj:
            arr.extend( flatten( i ) )
        return arr
    except:
        return obj

def consecutive_segments( lst ):
    lists = []
    newlist = []
    for idx, item in enumerate(lst):
        #print lst, lists, newlist
        if idx == 0:
            newlist.append(item)
            continue
        elif item - lst[idx-1] > 1:
            lists.append(newlist)
            newlist = []
            newlist.append(item)
            continue
        else: newlist.append(item)
    lists.append(newlist)
        
    if lists[-1] == []: return lists[:-1]
    return lists

def get_fasta_entities(sequence_for_fasta):
    fasta_entities = []
    i = 0
    spacers = [ ',', ' ' ]
    while i < len(sequence_for_fasta):
        if sequence_for_fasta[i] in spacers: continue
        if i == len(sequence_for_fasta) - 1 or sequence_for_fasta[i+1] != '[':
            # simple case
            fasta_entities.append(sequence_for_fasta[i])
        else:
            entity = sequence_for_fasta[i]
            i += 1
            while sequence_for_fasta[i] != ']':
                entity += sequence_for_fasta[i]
                i += 1
            entity += ']'
            fasta_entities.append(entity)

        i += 1
    #print fasta_entities
    return fasta_entities

def parse_ss( fasta_entities, line ):
    left_brackets = []
    pair_map = {}
    all_pairs = []

    complement = {'a':['u', 'X[OMU]', 'X[PSU]', 'X[5MU]'], 'u':['a','g', 'X[OMG]'], 'c':['g'], 'g':['c','X[5MC]', 'X[OMC]', 'u', 'X[OMU]', 'X[PSU]', 'X[5MU]']}
    left_bracket_char = '('
    right_bracket_char = ')'
    spacers = [ ',', ' ' ]
    count = 0
    for c in line:
        if c in spacers: continue
        count += 1
        if c == left_bracket_char:  left_brackets.append( count )
        if c == right_bracket_char:
            if len( left_brackets ) == 0:
                print "ValidationError: Number of right brackets does not match left brackets"
                exit()
            res1 = left_brackets[-1]
            res2 = count
            del left_brackets[-1]
            pair_map[ res1 ] = res2
            pair_map[ res2 ] = res1
            #all_pairs.append( [res1,res2] )
            if fasta_entities[ res2-1 ] in complement.keys() and len( fasta_entities ) > 0 and not fasta_entities[res1-1] in complement[ fasta_entities[res2-1] ]:
                print "ValidationError: Not complementary at positions %s%d and %s%d!"  % (fasta_entities[res1-1],res1,fasta_entities[res2-1],res2)
                exit()

    if len (left_brackets) > 0:
        print "ValidationError: Number of right brackets does not match left brackets"
        exit()
    return pair_map#, all_pairs

def get_stems( line, sequence_for_fasta ):
    """
    Not concerned with handling pseudoknots for now so let's just grab stems from our
    SS and sequence and get out of here. We are assuming that all stems are, in this
    sense, 'coequal.'

    A reasonable approach eventually would be to mandate that () is reserved for non-PK
    stems.
    """

    chainbreak_pos = []
    
    # THIS, not the sequence itself, should be equal in size to the secstruct etc.
    fasta_entities = get_fasta_entities(sequence_for_fasta)
    
    pair_map = parse_ss( fasta_entities, line )
    numres = len(line)

    # Parse out stems
    already_in_stem = { i: 0 for i in xrange( numres ) }
    
    stems = []
    for i in xrange( 1, numres + 1 ):
        if not pair_map.has_key( i ) or already_in_stem[ i ]: continue

        k = i 
        stem_res = []

        stem_res.append( [k, pair_map[k]] )
        already_in_stem[ k ] = 1
        already_in_stem[ pair_map[k] ] = 1

        while pair_map.has_key( k + 1 ) and pair_map[ k+1 ] == pair_map[ k ] - 1 and not already_in_stem[k+1] and not k in chainbreak_pos and not pair_map[k+1] in chainbreak_pos:
            k += 1
            stem_res.append( [k, pair_map[k]] )
            already_in_stem[ k ] = 1
            already_in_stem[ pair_map[k] ] = 1

        stems.append( stem_res )
    return stems

def loop_interpolate( x1,y1,x2,y2, fraction_of_circle, fraction_done_with_loop, add_angle=0, loop_sense=1 ):
    """
    Suppose (x1,y1) and (x2,y2) are two points separated by fraction_of_circle of
    arc. Then, return points fraction_done_with_loop along that arc between those
    two points.
    """
    #    C_________D
    #    /\       /)|
    #   /  \     /) |   find coordinates of D
    #  /    \   )   |   given its fraction of total angle
    # /______\_/____|      away from B or whatever.
    #A    M   B     E
    #        x1,y1   xD,y(AB*mAE/mAB)

    # Nonsensically, loops are going in the wrong direction.
    # Equivalent to careful debugging of below is just taking the complement.

    fraction_done_with_loop = 1 - fraction_done_with_loop

    mMCB = (360-fraction_of_circle*360)/2
    mAB = math.sqrt( (x2-x1)**2 + (y2-y1)**2 )
    mMC = 0.5 * mAB * math.tan(math.radians(mMCB))
    mBC = 0.5 * mAB / math.sin(math.radians(mMCB))
    r = mBC

    # Thus all points satisfy x^2 + y^2 = mBC^2
    mBCD = fraction_done_with_loop * fraction_of_circle*360

    # MC perp MB, so C can be found: mMC in perp dir from av A, B
    unit = None
    if y1 == y2:
        unit = [0,1]
    elif x1 == x2:
        unit = [1,0]
    else:
        print y2-y1
        print x2-x1
        unit = [1,float(-1)/float((y2-y1)/float(x2-x1))]
    xC,yC = [ (unit[0] * mMC) + (x1+x2)/2, (unit[1] * mMC) + (y1+y2)/2]
    #print xC, yC


    # xC + r * cos a
    # yC + r * sin a
    # a is angle from horizontal right
    # AMW TODO: for slightly nudged stems, this won't be -90
    # Add add_angle, which is 180 for flipped stems
    a = loop_sense*(-90 + add_angle + mMCB + mBCD)
    #print "%0.2f %0.2f" % (xC + r * math.cos(math.radians(a)), yC + r * math.sin(math.radians(a)))
    return [ xC + r * math.cos(math.radians(a)), yC + r * math.sin(math.radians(a))]


# constants used below. May end up functions of font.
NT_DISTANCE = 20
NT_MIN_DISTANCE = 10
NT_MAX_DISTANCE = 25
NT_MIN_ANGLE    = 90
NT_MAX_ANGLE    = 180 # maybe should be 300... hm.
# Imagine having different spring constants for stem vs loop vs PK
#spring_constant = .1 suitable for harmonic, but not ideal
spring_constant = 1
angular_spring_constant = 1.0/15.0

def length_score( canvas ):
    """
    1. Nucleotides should be a particular distance from what's adjacent to them in sequence.
    form shouldn't be harmonic. Many distances are 'fine' -- say, 10-20 -- so we need a 'flat harmonic'
    """
    score = 0
    for seqpos, nt in canvas.nucleotides.iteritems():
        if seqpos + 1 not in canvas.nucleotides.keys(): continue

        d = distance( nt, canvas.nucleotides[seqpos+1] )
        #print "Distance between %d and %d is %f" % (seqpos, seqpos+1, d)
        #score += harmonic_penalty( d, NT_DISTANCE, spring_constant )
        score += flat_harmonic_penalty( d, NT_MIN_DISTANCE, NT_MAX_DISTANCE, spring_constant )

    return score


def angle_score( canvas ):
    """
    1. Nucleotides should be a particular distance from what's adjacent to them in sequence.
    form shouldn't be harmonic. Many distances are 'fine' -- say, 10-20 -- so we need a 'flat harmonic'
    """
    score = 0
    for seqpos, nt in canvas.nucleotides.iteritems():
        if seqpos + 1 not in canvas.nucleotides.keys(): continue
        if seqpos + 2 not in canvas.nucleotides.keys(): continue
        a = angle( nt, canvas.nucleotides[seqpos+1], canvas.nucleotides[seqpos+2] )
        #print "Angle between %d and %d and %d is %f" % (seqpos, seqpos+1, seqpos+2, a)
        # Don't worry about angular harmonic func at the moment.
        score += flat_harmonic_penalty( a, NT_MIN_ANGLE, NT_MAX_ANGLE, angular_spring_constant )
    return score

def score( canvas ):
    """
    Assign a score to a configuration of nucleotides.
    TODO: update to use relative coords or something?
    """

    score = 0

    score += length_score(canvas)
    score += angle_score(canvas)

    # 1. Chainbreaks for all junction loops.

    # 2. Nucleotides should be at least 15 from ALL other nts. Don't double-count.
    for seqpos1, nt1 in canvas.nucleotides.iteritems():
        for seqpos2, nt2 in canvas.nucleotides.iteritems():
            if seqpos1 >= seqpos2: continue
            #if distance( nt1, nt2 ) < 15: score += 100
            # Don't do a square root you don't have to
            if distance_squared( nt1, nt2 ) < 225: score += 100


    # 2b. Depending on base pair width, it's possible that the 'best' orientation for two stacks is
    # for them to be overlapping...
    #  nt --- nt
    #  |   nt |-- nt
    #  nt -|- nt  |
    #  |   nt |-- nt
    # so, penalize this ( 'between two bp'ed nts' penalized like hitting a nt )
    for seqpos, nt in canvas.nucleotides.iteritems():
        for stem in canvas.stems:
            for bp in stem.base_pairs:
                if nt == bp.nt1 or nt == bp.nt2:
                    continue
                # if nt in question is 'in between' nt1 and nt2. how to judge?
                if nt.x > bp.nt1.x and nt.x < bp.nt2.x \
                    and closer_than( nt.y, bp.nt1.y, canvas.bp_offset_height*0.8 ):
                    score += 100

            # Look at pairs of adjacent nts in stem
            for seqpos1, nt1 in stem.nucleotides.iteritems():
                if seqpos1 + 1 in stem.nucleotides.keys():
                    nt2 = stem.nucleotides[ seqpos1 + 1 ]
                    if nt.y > nt1.y and nt.y < nt2.y \
                        and closer_than( nt.x, nt1.x, canvas.bp_offset_width*0.8 ):
                        score += 100

    # 3. Stems should be straight and vertical.
    ### Resolved via kinematics?

    # 4. Coaxial stacks, where known, should be colinear.
    for idx1, idx2 in canvas.list_of_coaxialities:
        # idx1 and idx2 are coaxial stems.
        # Any two BPs from these stems should do, but -- because I have no idea
        # how we intend to handle 'alternate-order' stacks (where the first nt in a bp
        # is stacked on the second nt in another), compare either direction.
        bp1 = canvas.stems[idx1].base_pairs[0]
        bp2 = canvas.stems[idx2].base_pairs[0]
        penalty1 = harmonic_penalty( bp1.nt1.x - bp2.nt1.x, 0, 1 ) + harmonic_penalty( bp1.nt2.x - bp2.nt2.x, 0, 1 )
        penalty2 = harmonic_penalty( bp1.nt1.x - bp2.nt2.x, 0, 1 ) + harmonic_penalty( bp1.nt2.x - bp2.nt1.x, 0, 1 )
        score += min( penalty1, penalty2 )


    # 5. Apical loops should be symmetrical and near a circular path.
    for loop in canvas.loops:
        for seqpos in loop.nucleotides.keys():
            pass

    return score