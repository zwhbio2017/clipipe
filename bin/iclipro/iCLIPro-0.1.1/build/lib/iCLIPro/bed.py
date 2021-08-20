"""Simple interface for handling bedGraphs.
"""

import bisect
import os

def nc_open(fname, omode):
    import gzip
    if fname.endswith('.gz'):
        base_fname = os.path.basename(fname)
        return gzip.GzipFile(base_fname, fileobj=open(fname, omode), mode=omode)
    return open(fname, omode)

def save(fname_out, bg, track_name, db_name):
    fout = nc_open(fname_out, "wt")
    for target_strand in ['+', '-']:
        fout.write('track type=bedGraph name="%s strand%s" description="%s strand%s" db=%s\n' % (track_name, target_strand, track_name, target_strand, db_name))
        for (chrome, strand), by_pos in sorted(bg.iteritems()):
            if strand <> target_strand:
                continue
            for pos, v in sorted(by_pos.iteritems()):
                fout.write("%s\t%s\t%s\t%s%s\n" % (chrome, pos, pos+1, strand, v))
    fout.close()

def load(fname_in):
    retd = {}
    for r in nc_open(fname_in, "rt"):
        if r.startswith('track'):
            continue
        rs = r.rstrip('\n\r').split('\t')
        chrome, pos_sta, pos_end, strand_val = rs
        pos_sta = int(pos_sta)
        pos_end = int(pos_end)
        if strand_val[0] in ['-', '+']:
            strand = strand_val[0]
            val = int(strand_val[1:])
        else:
            strand = '+'
            val = int(strand_val)
        retd.setdefault((chrome, strand), {}).setdefault((pos_sta, pos_end), val)
    return retd

def cluster(bg_in, cl_size):
    if cl_size == 0:
        return bg_in

    bg_out = {}
    for k, by_pos in bg_in.iteritems():
        tmp_k_out = bg_out.setdefault(k, {})

        # assign neighbors to their closest peaks
        pos2cl = {}
        for _, top_pos in sorted([(counts, pos) for pos, counts in by_pos.iteritems()]):
            for pos in xrange(top_pos-cl_size, top_pos+cl_size+1):
                pos2cl.setdefault(pos, top_pos)

        # assign counts of neighbors to their closest peaks
        for pos, counts in by_pos.iteritems():
            top_pos = pos2cl[pos]
            tmp_k_out[top_pos] = tmp_k_out.get(top_pos, 0) + counts
            
    return bg_out

def stats(bg):
    tot_sites = 0
    tot_counts = 0
    for _, by_pos in bg.iteritems():
        tot_sites += len(by_pos)
        tot_counts += sum([v for v in by_pos.values()])
    return tot_sites, tot_counts

def compare(b_center, b_other, ws):
    retd = {}
    for (chrome, strand), c_by_pos in b_center.iteritems():
        o_by_pos = b_other.get((chrome, strand), {})
        if strand == '+':
            i = 1
        else:
            i = -1

        #c_poss = sorted(c_by_pos.keys())
        o_poss = sorted(o_by_pos.keys())

        for p_c, count_c in c_by_pos.iteritems():
            scan_up = p_c - ws
            scan_down = p_c + ws
            i_up = bisect.bisect_left(o_poss, scan_up)
            i_down = bisect.bisect(o_poss, scan_down)
            for p_o in o_poss[i_up:i_down]:
                d = (p_o - p_c)*i
                retd[d] = retd.get(d, 0) + o_by_pos[p_o]

    return retd

