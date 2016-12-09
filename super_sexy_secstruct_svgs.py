import svgwrite


def draw_stems(bps, dwg):
    # draw basepairs
    for bp_idx, (bp1, bp2) in enumerate(bps): 
        # need to compute some offet? ...
        y_offset = bp_idx * 2.
        dwg.add(dwg.line((0, 0), (0, y_offset), 
                stroke=svgwrite.rgb(10, 10, 16, '%')))
        dwg.add(dwg.text(bp1, insert=(0.2, y_offset), fill='blue'))
        dwg.add(dwg.text(bp2, insert=(0.8, y_offset), fill='yellow'))
    return dwg


def draw_secstruct(bps, file_name='default.out'):
    # init drawing
    dwg = svgwrite.Drawing(file_name)

    # draw struct here
    dwg = draw_stems(bps, dwg)

    # save drawing
    dwg.save()
    return dwg



if __name__=="__main__":
    fn = 's_s_ss.svg'
    bps = [('C', 'G'), ('A', 'U')]
    dwg = draw_secstruct(bps, file_name=fn)
