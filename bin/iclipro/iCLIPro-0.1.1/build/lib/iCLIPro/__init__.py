"""iCLIPro is a package that tests iCLIP read overlaps.
It can aid in the indentification and correction of
potential systematic errors in iCLIP data.

For documentation, see http://www.biolab.si/iCLIPro

Documentation source can be found in folder docs.

Examples are provided in folder examples.
"""

from _version import __version__

import iCLIPro.bed
import iCLIPro.cmd_utils

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

matplotlib.rcParams['font.size'] = 10 #'small'
matplotlib.rcParams['font.sans-serif'] = 'Arial'
matplotlib.rcParams['axes.labelsize'] = 10 #'small'
matplotlib.rcParams['axes.titlesize'] = 10
matplotlib.rcParams['xtick.labelsize'] = 8 #'small'
matplotlib.rcParams['ytick.labelsize'] = 8 #'small'
matplotlib.rcParams['legend.fontsize'] = 9

import os
import pysam
import numpy

def amb_match(s1, s2):
    #assert len(s1) == len(s2)
    cn = sum([(c1 == 'N' or c2 == 'N' or c1 == c2) for c1, c2 in zip(s1, s2)])
    return len(s1) - cn <= 2

def collapse_bcs(reads, strand):
    retl = []
    # group reads by starting (xlink_pos)
    by_pos = {}
    for xlink_pos, middle_pos, end_pos, bc in reads:
        by_pos.setdefault(xlink_pos, {}).setdefault(bc, []).append((middle_pos, end_pos))

    # reassign reads with ambiguous rnd barcodes to most frequent unambiguous matching rnd barcode
    # determine all unambiguous rnd barcodes
    for xlink_pos, by_bc in by_pos.iteritems():
        accepted_bcs = []
        ambig_bcs = []
        for bc, hits in by_bc.iteritems():
            Ns = bc.count("N")
            if Ns == 0:
                accepted_bcs.append((len(hits), bc))
            else:
                ambig_bcs.append((Ns, bc))
        # sort unambiguous by decreasing frequency (number of reads it is present)
        accepted_bcs = [bc for _, bc in sorted(accepted_bcs, reverse=True)]

        # add those ambiguous (in increasing order of the number of Ns in rnd barcode)
        ambig_bcs.sort()
        for _, abc in ambig_bcs:
            matches = False
            for bc in accepted_bcs:
                if amb_match(abc, bc):
                    matches = True
                    # move read records from abc under bc
                    by_bc[bc].extend(by_bc[abc])
                    del by_bc[abc]
                    break
            if not matches:
                accepted_bcs.append(abc)

        for bc, middle_end_poss in by_bc.iteritems():
            middle_end_poss.sort()
            i = len(middle_end_poss)
            if i % 2 == 0:
                i = i / 2
                if strand == '-':
                    i = i - 1
            else:
                i = i / 2
            middle_pos, end_pos = middle_end_poss[i]
            retl.append((xlink_pos, middle_pos, end_pos))

    return retl

def _read_bam(in_bam_fname, mapq_th, split_flank, bin_size=1000):
    ### read file
    try:
        samfile = pysam.Samfile(in_bam_fname, "rb")
    except:
        print "Error: could not open BAM file: %s" % in_bam_fname
        return None
    # can we rely on direct indexing to get chrome name
    assert [samfile.getrname(i) == rname for i, rname in enumerate(samfile.references)]

    # statistics counters
    all_recs = 0
    notmapped_recs = 0
    mapped_recs = 0
    lowmapq_recs = 0 
    # ignore reads with gaps at beginning or end of reads,
    # because hard to tell if exon-intron or exon-exon
    split_flank_ignored_recs = 0
    used_recs = 0

    # assign reads to bins
    valid_bc_nucs = set('ATCGN')
    bins = {}
    all_bcs = {} # list of all random barcodes and used to reduce memory use
    all_bcs_cn = {} # list of all random barcodes and used to reduce memory use
    split_flank_distance = split_flank - 1
    for r in samfile:
        all_recs += 1
#        if all_recs > 50000: break #DEBUG
        if r.is_unmapped:
            notmapped_recs += 1
            continue

        if r.mapq < mapq_th:
            lowmapq_recs += 1
            continue

        mapped_recs += 1

        ## identify gaps
        # list of positions to which read mapped to
        poss = sorted(r.positions)

        # ignore reads that have gaps whithin the first/last (split_flank) nucleotides
        # longest gap between consecutive nucleotides
        gaps = [y-x for x, y in zip(poss, poss[1:])]

        # two consecutive nucleotides are 1 apart
        if max(gaps[:split_flank_distance]+[0] + gaps[-split_flank_distance:]) > 1:
            split_flank_ignored_recs += 1
            continue

        ## use this read record because it passes all filtering steps
        used_recs += 1

        ## random barcode must be part of read name
        if ':rbc' in r.qname:
            bc = r.qname.rsplit(':rbc', 1)[1].split(':',)[0]
        elif ':' in r.qname:
            bc = r.qname.rsplit(':', 1)[1]
            if set(bc) - valid_bc_nucs:
                # barcode contains invalid characters
                bc = ''
        else:
            bc = ''
        bc = all_bcs.setdefault(bc, bc)
        all_bcs_cn[bc] = all_bcs_cn.get(bc, 0) + 1

        ## position of cross-link is one nucleotide before start of read
        if r.is_reverse:
            strand = '-'
            xlink_pos = poss[-1] + 1
            end_pos = poss[0]
        else:
            strand = '+'
            xlink_pos = poss[0] - 1
            end_pos = poss[-1]
        chrome = samfile.references[r.tid]

        ## middle position is position of middle nucleotide
        # (should not take just middle of mapped region, because we can have spliced reads)
        middle_pos = poss[len(poss)/2]
        read_len = len(r.seq)

        ## place reads in bins
        bin_index = xlink_pos / bin_size
        bins.setdefault((chrome, strand), {}).setdefault(bin_index, []).append((xlink_pos, middle_pos, end_pos, bc, read_len))
    samfile.close()
    return all_recs, notmapped_recs, mapped_recs, lowmapq_recs, split_flank_ignored_recs, used_recs, all_bcs_cn, bins

def nice_lens_str(len_list):
    len_list = sorted(len_list)

    nice_lst = []
    pv = len_list[0]
    cs = pv
    for v in len_list[1:]:
        if v - pv > 1:
            nice_lst.append((cs, pv))
            cs = v
        pv = v
    nice_lst.append((cs, pv))
    nice_str = []
    for start_val, end_val in nice_lst:
        if start_val == end_val:
            nice_str.append(str(start_val))
        else:
            nice_str.append("%s..%s" % (start_val, end_val))
    return ",".join(nice_str)

def overlap_test(in_bam_fname, output_folder, bin_size, bin_min_reads, overlap_ratio_flank,
                 readlen2groups, group_comparisons,
                 mapq_th=10, split_flank=5, group_comparisons_flank=60):

    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)

    output_folder_beds = "%s/beds" % (output_folder)
    if not os.path.isdir(output_folder_beds):
        os.makedirs(output_folder_beds)

    bam_basename = os.path.basename(in_bam_fname).rsplit('.', 1)[0]
    output_folder_pref = output_folder.split('/')[-1]
    out_fname_pref = "%s/%s_%s" % (output_folder, output_folder_pref, bam_basename)
    out_fname_beds_pref = "%s/%s_%s" % (output_folder_beds, output_folder_pref, bam_basename)

    all_recs, notmapped_recs, mapped_recs, lowmapq_recs, split_flank_ignored_recs, used_recs, all_bcs, bins = _read_bam(in_bam_fname, mapq_th, split_flank, bin_size=bin_size)

    ### identify bins with enough reads
    tot_bins = 0
    tot_reads = 0
    selected_reads = 0
    selected_reads_len_distrib = {}
    selected_bins = 0
    bins_to_remove = []
    for (chrome, strand), by_bin_index in bins.iteritems():
        tot_bins += len(by_bin_index)

        for bin_index, reads in by_bin_index.iteritems():
            reads_in_bin = len(reads)

            tot_reads += reads_in_bin
            if reads_in_bin >= bin_min_reads:
                selected_bins += 1
                selected_reads += reads_in_bin
                for _, _, _, _, read_len in reads:
                    selected_reads_len_distrib[read_len] = selected_reads_len_distrib.get(read_len, 0) + 1
            else:
                bins_to_remove.append((chrome, strand, bin_index))
    max_read_len = max(selected_reads_len_distrib.keys())
    # remove bins that are not needed anymore
    for chrome, strand, bin_index in bins_to_remove:
        tmpd = bins[(chrome, strand)]
        del tmpd[bin_index]

    ### generate bedGraphs for each group
    g2lens = {}
    for read_len, read_grps in readlen2groups.iteritems():
        for read_grp in read_grps:
            g2lens.setdefault(read_grp, []).append(read_len)

    grp_bed_start = {}
    grp_bed_middle = {}
    grp_bed_end = {}
    for read_grp, valid_read_lens in g2lens.iteritems():
        bed_start = grp_bed_start.setdefault(read_grp, {})
        bed_middle = grp_bed_middle.setdefault(read_grp, {})
        bed_end = grp_bed_end.setdefault(read_grp, {})

        for (chrome, strand), by_bin_index in bins.iteritems():
            tmpd_bed_start = bed_start.setdefault((chrome, strand), {})
            tmpd_bed_middle = bed_middle.setdefault((chrome, strand), {})
            tmpd_bed_end = bed_end.setdefault((chrome, strand), {})

            for bin_index, reads in by_bin_index.iteritems():
                valid_reads = [(xlink_pos, middle_pos, end_pos, bc) for xlink_pos, middle_pos, end_pos, bc, read_len in reads if read_len in valid_read_lens]

                for xlink_pos, middle_pos, end_pos in collapse_bcs(valid_reads, strand):
                    tmpd_bed_start[xlink_pos] = tmpd_bed_start.get(xlink_pos, 0) + 1
                    tmpd_bed_middle[middle_pos] = tmpd_bed_middle.get(middle_pos, 0) + 1
                    tmpd_bed_end[end_pos] = tmpd_bed_end.get(end_pos, 0) + 1
    # bins not needed after this point
    del bins

    ### save bedGraphs
    for read_grp, bed_start in sorted(grp_bed_start.iteritems()):
        iCLIPro.bed.save("%s_%s_start.bed" % (out_fname_beds_pref, read_grp), bed_start, "%s start %s" % (read_grp, bam_basename), db_name="")
    for read_grp, bed_middle in sorted(grp_bed_middle.iteritems()):
        iCLIPro.bed.save("%s_%s_center.bed" % (out_fname_beds_pref, read_grp), bed_middle, "%s center %s" % (read_grp, bam_basename), db_name="")
    for read_grp, bed_end in sorted(grp_bed_end.iteritems()):
        iCLIPro.bed.save("%s_%s_end.bed" % (out_fname_beds_pref, read_grp), bed_end, "%s end %s" % (read_grp, bam_basename), db_name="")

    ### generate read overlap maps
    fout_graph_pdf_fname = "%s_graphs.pdf" % out_fname_pref
    fout_graph_png_fname = "%s_graphs.png" % out_fname_pref
    fout_graph_tab_fname = "%s_graphs.tab" % out_fname_pref
    fout_heatmap_png_fname = "%s_graphs_heatmap_ref%%s.png" % out_fname_pref
    fout_heatmap_pdf_fname = "%s_graphs_heatmap_ref%%s.pdf" % out_fname_pref
    fout_graph_tab = open(fout_graph_tab_fname, "wt")
    header = ["reads_group", "reference_group", "origin"]
    header = header + range(-group_comparisons_flank, group_comparisons_flank+1)
    fout_graph_tab.write("%s\n" % "\t".join([str(x) for x in header]))
    graph_color = '#0081c5'
    dy = 3
    dx = len(group_comparisons)
    plt.figure(figsize=(dx*1.9, dy*1.7))
    plt.clf()
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.12, top=0.9, wspace=0.55, hspace=0.7)
    si = 1
    heatmaps = {} # group comparisons by ref_grp, plot one heatmap for each
    read_groups_by_ref = {}
    for read_grp, ref_grp in group_comparisons:
        read_groups_by_ref.setdefault(ref_grp, []).append(read_grp)
        si_start = si
        ylim = 0.0
        for origin, grp_bed in [('start', grp_bed_start), ('center', grp_bed_middle), ('end', grp_bed_end)]:
            bed_cur = grp_bed.get(read_grp, {})
            bed_ref = grp_bed.get(ref_grp, {})

            pos_distrib = iCLIPro.bed.compare(bed_ref, bed_cur, group_comparisons_flank)
            if pos_distrib:
                xs, ys = zip(*sorted(pos_distrib.iteritems()))
            else:
                xs, ys = [], []

            ax = plt.subplot(dy, dx, si)
            plt.bar(xs, ys, width=0.7, ec=graph_color, color=graph_color)
            plt.xlabel("offset of %s to %s" % (read_grp, ref_grp))
            ax.get_xaxis().set_label_coords(0.5, -0.22)
            if (si - 1) % dx == 0:
                plt.ylabel("%s\nfrequency" % origin, multialignment='center', verticalalignment='center')
                ax.get_yaxis().set_label_coords(-0.4, 0.5)
            ylim = max(ylim, max(plt.ylim()))

            freqs = [0]*(2*group_comparisons_flank+1)
            for x, y in zip(xs, ys):
                assert -group_comparisons_flank <= x <= group_comparisons_flank
                freqs[group_comparisons_flank + x] = y
            olst = [read_grp, ref_grp, origin] + freqs
            fout_graph_tab.write("%s\n" % "\t".join([str(x) for x in olst]))
            si += dx

            max_freq = float(max(freqs))
            if max_freq == 0.0:
                max_freq = 1.0
            # normalize by row
            freqs = [v/max_freq for v in freqs]
            heatmaps.setdefault(ref_grp, {}).setdefault(origin, []).append(freqs)

        # put yaxis to same scale
        ylim = (0.0, ylim)
        si = si_start

        plt.subplot(dy, dx, si)
        plt.ylim(ylim)
        yticks_vals, _ = plt.yticks()
        yticks_max = yticks_vals[-1]
        yticks = [0, yticks_max*1.0/4.0, yticks_max*2.0/4.0, yticks_max*3.0/4.0, yticks_max]

        xticks_min = -group_comparisons_flank
        xticks_max = group_comparisons_flank
        xticks = [xticks_min, xticks_min*1.0/2.0, 0, xticks_max*1.0/2.0, xticks_max]
     
        si = si_start 
        for origin, grp_bed in [('start', grp_bed_start), ('center', grp_bed_middle), ('end', grp_bed_end)]:
            plt.subplot(dy, dx, si)
            plt.plot([0, 0], ylim, c='k', lw=0.4, alpha=0.7)
            plt.ylim(ylim)
            plt.xlim((-group_comparisons_flank, group_comparisons_flank))
            plt.yticks(yticks, ["%0.0f" % y for y in yticks])
            plt.xticks(xticks, ["%0.0f" % x for x in xticks])
            si += dx

        si = si_start + 1
    fout_graph_tab.close()

    savefig_error = False

    savefig_error_png = False
    try:
        plt.savefig(fout_graph_png_fname)
    except:
        savefig_error = True
        savefig_error_png = True
        print "ERROR: problems saving figure: '%s'" % fout_graph_png_fname
        print "       please check your matplotlib installation"
        if os.path.isfile(fout_graph_png_fname):
            os.remove(fout_graph_png_fname)

    savefig_error_pdf = False
    try:
        plt.savefig(fout_graph_pdf_fname)
    except:
        savefig_error = True
        savefig_error_pdf = True
        print "ERROR: problems saving figure: '%s'" % fout_graph_pdf_fname
        print "       please check your matplotlib installation"
        if os.path.isfile(fout_graph_pdf_fname):
            os.remove(fout_graph_pdf_fname)
    plt.clf()

    #cmap = matplotlib.colors.ListedColormap(['Blue', 'White', 'Yellow', 'Red'], 11)
    #cmap = matplotlib.colors.LinearSegmentedColormap.from_list('BWYR', ['Blue', 'White', 'Orange', 'Red'], gamma=0.3)
    #cmap = matplotlib.colors.LinearSegmentedColormap.from_list('BYOR', ['Blue', 'LightYellow', 'Orange', 'Red'], N=100, gamma=0.85)
    #cmap = matplotlib.colors.LinearSegmentedColormap.from_list('BYOR', ['DarkBlue', 'White', 'LightYellow', 'Orange', 'Red'], N=100, gamma=0.93)
    #cmap = matplotlib.colors.LinearSegmentedColormap.from_list('BYOR', ['DarkBlue', 'LightYellow', 'Orange', 'Red'], N=100, gamma=0.93)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('BYOR', ['DarkBlue', 'White', 'LightYellow', 'Orange', 'Red'], N=100, gamma=0.91)
    #cmap = plt.get_cmap("RdYlBu_r", 7)

    ### generate heatmaps, one file for each ref_grp
    for ref_grp, by_origin in heatmaps.iteritems():
        dy = len(by_origin)
        dx = 1

        read_groups = read_groups_by_ref[ref_grp]
        fig = plt.figure(figsize=(0.14*dx*(group_comparisons_flank*2+1), dy*len(read_groups)*0.14))
        plt.clf()
        plt.subplots_adjust(left=0.1, right=0.9, bottom=0.15, top=0.94, wspace=0.12, hspace=0.22)
        si = 1
        for origin in ['start', 'center', 'end']:
            hm = by_origin[origin]
            # scalling of whole matrix
            max_v = max([max(r) for r in hm])
            min_v = min([min(r) for r in hm])
            range_v = (max_v - min_v)
            if range_v == 0.0:
                range_v = 1.0
            #hm_scaled = [[(x-min_v)/range_v for x in r] for r in hm]
            hm_scaled = [[(x-min(r))/(max(r)-min(r)) for x in r] for r in hm]

            by_origin[origin] = hm_scaled

            hm = numpy.array(hm_scaled[::-1])
            ax = plt.subplot(dy, dx, si, frameon=False)
            ax.spines['left'].set_position(('outward', 10))
            ax.spines['bottom'].set_position(('outward', 10))
            # Hide the right and top spines
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            # Only show ticks on the left and bottom spines
            ax.yaxis.set_ticks_position('left')
            ax.xaxis.set_ticks_position('bottom')

            im = plt.pcolor(hm, linewidths=0.5, shading='faceted', edgecolors='gray', cmap=cmap, vmin=0.0, vmax=1.0)
            #origin='upper', aspect='equal', cmap = plt.get_cmap("RdYlBu_r"), interpolation='nearest', extent=(-group_comparisons_flank, group_comparisons_flank+1, 0, len(group_comparisons)))

            plt.xlim(0, 2*group_comparisons_flank+1)
            xticks = [x for x in range(-group_comparisons_flank, 0, 10)]
            xticks = xticks + [0] + [-v for v in xticks]
            if origin == 'end':
                xlabels = ["%0.0f" % x for x in xticks]
                plt.xlabel("offset of sites\ncompared to reads of group %s (%s)" % (ref_grp, nice_lens_str(g2lens[ref_grp])))
            else:
                xlabels = ['']*len(xticks)
            plt.xticks([group_comparisons_flank+x+0.5 for x in xticks], xlabels)

            plt.ylim(0, len(read_groups))
            plt.yticks([y + 0.5 for y in range(len(read_groups))], ["%s (%s)" % (read_grp, nice_lens_str(g2lens[read_grp])) for read_grp in read_groups[::-1]], fontsize='xx-small')

            plt.title("%s - %s" % (bam_basename, origin))
            si += 1

        fig.subplots_adjust(right=0.85)
        cbar_ax = fig.add_axes([0.9, 0.75, 0.02, 0.12])
        cbar = fig.colorbar(im, cax=cbar_ax, orientation='vertical', ticks=[0.0, 0.25, 0.5, 0.75, 1.0], format="%0.2g", drawedges=False)
        cbar.solids.set_edgecolor("face")

        fout_heatmap_ref_png_fname = fout_heatmap_png_fname % ref_grp
        savefig_error_png = False
        try:
            plt.savefig(fout_heatmap_ref_png_fname, bbox_inches='tight')
        except:
            savefig_error = True
            savefig_error_png = True
            print "ERROR: problems saving figure: '%s'" % fout_heatmap_ref_png_fname
            print "       please check your matplotlib installation"
            if os.path.isfile(fout_heatmap_ref_png_fname):
                os.remove(fout_heatmap_ref_png_fname)

        fout_heatmap_ref_pdf_fname = fout_heatmap_pdf_fname % ref_grp
        savefig_error_pdf = False
        try:
            plt.savefig(fout_heatmap_ref_pdf_fname, bbox_inches='tight')
        except:
            savefig_error = True
            savefig_error_pdf = True
            print "ERROR: problems saving figure: '%s'" % fout_heatmap_ref_pdf_fname
            print "       please check your matplotlib installation"
            if os.path.isfile(fout_heatmap_ref_pdf_fname):
                os.remove(fout_heatmap_ref_pdf_fname)
        plt.clf()

    ### main diagnostic
    ##
    ## we observe the ratio of read probabilities before and after center offset
    overlap_ratio_scores = []
    for ref_grp, by_origin in heatmaps.iteritems():
        read_groups = read_groups_by_ref[ref_grp]
        hm = by_origin['center']
        ref_grp_read_len = min(g2lens[ref_grp])
        ref_half_len = ref_grp_read_len/2

        site_ratios = []
        for read_grp, hm_row in zip(read_groups, hm):
            read_grp_read_len = min(g2lens[read_grp])
            offset_middle_max = int(ref_half_len - read_grp_read_len/2.0) # take largest overlap when read_grp_read_len is odd
            lower_limit_offset = (group_comparisons_flank - 1) - offset_middle_max - overlap_ratio_flank - 1
            # skip center
            upper_limit_offset = (group_comparisons_flank + 1) + offset_middle_max + overlap_ratio_flank + 1

            left_part = hm_row[lower_limit_offset:group_comparisons_flank]
            left_of_center_offset_middle_max = sum(left_part)
            right_part = hm_row[group_comparisons_flank+1:upper_limit_offset+1]
            right_of_center_offset_middle_max = sum(right_part)
            assert len(left_part) == len(right_part)

            if right_of_center_offset_middle_max == 0.0:
                right_of_center_offset_middle_max = 1.0
            ratio = left_of_center_offset_middle_max/right_of_center_offset_middle_max
            site_ratios.append(ratio)

        median_site_ratio = numpy.median(site_ratios)
        mean_site_ratio = numpy.mean(site_ratios)
        SD_site_ration = numpy.std(site_ratios)

        overlap_ratio_scores.append((ref_grp, ref_grp_read_len, median_site_ratio, mean_site_ratio, SD_site_ration))

    ### generate report
    fout_report_fname = "%s_report_iCLIPro.txt" % out_fname_pref
    fout_report = open(fout_report_fname, "wt")

    fout_report.write("# Parameters\n")
    fout_report.write("#   in_bam_fname: %s\n" % in_bam_fname)
    fout_report.write("#   output_folder: %s\n" % output_folder)
    fout_report.write("#   bin_size: %s\n" % bin_size)
    fout_report.write("#   bin_min_reads: %s\n" % bin_min_reads)
    fout_report.write("#   overlap_ratio_flank: %s\n" % overlap_ratio_flank)
    fout_report.write("#   mapq_th: %s\n" % mapq_th)
    fout_report.write("#   split_flank: %s\n" % split_flank)
    fout_report.write("\n")

    fout_report.write("# Read filtering stats\n")
    fout_report.write("#   total number of records: %s\n" % all_recs)
    fout_report.write("#   not mapped: %s (%0.2f%% total)\n" % (notmapped_recs, 100.0*notmapped_recs/all_recs))
    fout_report.write("#   mapped but ignored because MAPQ < %s: %s (%0.2f%% total)\n" % (mapq_th, lowmapq_recs, 100.0*lowmapq_recs/all_recs))
    fout_report.write("#   mapped and with MAPQ >= %s: %s (%0.2f%% total)\n" % (mapq_th, mapped_recs, 100.0*mapped_recs/all_recs))
    fout_report.write("#   additionally ignored because gap in first/last %s nt: %s (%0.2f%% total)\n" % (split_flank, split_flank_ignored_recs, 100.0*split_flank_ignored_recs/all_recs))
    fout_report.write("#   used reads: %s (%0.2f%% total)\n" % (used_recs, 100.0*used_recs/all_recs))
    
    show_only = 30
    all_bcs = sorted([(cn, bc) for bc, cn in all_bcs.iteritems()], reverse=True)
    all_bcs = [(bc, cn) for cn, bc in all_bcs]
    if len(all_bcs) > show_only:
        fout_report.write("#   distinct random barcodes (showing only %s out of all %s): %s\n" % (show_only, len(all_bcs), all_bcs[:show_only]))
    else:
        fout_report.write("#   distinct random barcodes (%s): %s\n" % (len(all_bcs), all_bcs))
    fout_report.write("\n") 

    fout_report.write("# Bin filtering stats\n")
    fout_report.write("#   total number of bins (of size %snt): %s\n" % (bin_size, tot_bins))
    fout_report.write("#   selected bins (with at least %s reads): %s (%0.2f%% total)\n" % (bin_min_reads, selected_bins, 100.0*selected_bins/tot_bins))
    fout_report.write("#   used reads: %s\n" % tot_reads)
    fout_report.write("#   selected reads (in selected bins): %s (%0.2f%% total)\n" % (selected_reads, 100.0*selected_reads/tot_reads))
    fout_report.write("#   selected reads length distribution and grouping:\n")
    fout_report.write("#\tread length\tnumber of reads\t% selected reads\tgroup(s)\n")
    notgrouped_reads = 0
    notgrouped_lengths = []
    for read_len, cn in sorted(selected_reads_len_distrib.iteritems()):
        if read_len not in readlen2groups:
            notgrouped_lengths.append(read_len)
            notgrouped_reads += cn
            grp_names = 'notused'
        else:
            grp_names = readlen2groups[read_len]
        fout_report.write("#\t%s\t%s\t%0.2f%%\t%s\n" % (read_len, cn, 100.0*cn/selected_reads, grp_names))
    if notgrouped_reads:
        fout_report.write("WARNING: %s reads (%0.2f%% selected) of lengths %s are not being used because not assigned to any group.\n" % (notgrouped_reads, 100.0*notgrouped_reads/selected_reads, notgrouped_lengths))
    fout_report.write("\n")

    fout_report.write("# Generated bedGraphs\n")
    fout_report.write("#\tgroup\torigin\ttotal cross-linked sites\ttotal cDNA\tfile name\n")
    for read_grp in sorted(g2lens.iterkeys()):
        bed_start = grp_bed_start.setdefault(read_grp, {})
        start_fname = "%s_%s_start.bed" % (out_fname_beds_pref, read_grp)
        start_num_sites, start_tot_cDNA = iCLIPro.bed.stats(bed_start)

        bed_middle = grp_bed_middle.setdefault(read_grp, {})
        middle_fname = "%s_%s_center.bed" % (out_fname_beds_pref, read_grp)
        middle_num_sites, middle_tot_cDNA = iCLIPro.bed.stats(bed_middle)

        bed_end = grp_bed_end.setdefault(read_grp, {})
        end_fname = "%s_%s_end.bed" % (out_fname_beds_pref, read_grp)
        end_num_sites, end_tot_cDNA = iCLIPro.bed.stats(bed_end)

        fout_report.write("#\t%s\t%s\t%s\t%s\t%s\n" % (read_grp, 'start', start_num_sites, start_tot_cDNA, start_fname))
        fout_report.write("#\t%s\t%s\t%s\t%s\t%s\n" % (read_grp, 'center', middle_num_sites, middle_tot_cDNA, middle_fname))
        fout_report.write("#\t%s\t%s\t%s\t%s\t%s\n" % (read_grp, 'end', end_num_sites, end_tot_cDNA, end_fname))
    fout_report.write("\n")

    fout_report.write("# Generated graphs for given comparisons\n")
    fout_report.write("#   graph data: %s\n" % (fout_graph_tab_fname))
    fout_report.write("#   graph: %s\n" % (fout_graph_pdf_fname))
    fout_report.write("#   graph: %s\n" % (fout_graph_png_fname))
    fout_report.write("\n")

    fout_report.write("# Overlap start site ratios (overlap_ratio_flank=%s)\n" % overlap_ratio_flank)
    fout_report.write("#\treference group\treference group read length\tmedian start site ratio\tmean start site ratio\tSD\n")
    for ref_grp, ref_grp_read_len, median_site_ratio, mean_site_ratio, SD_site_ration in overlap_ratio_scores:
        fout_report.write("#\t%s\t%s\t%0.3f\t%0.3f\t%0.3f\n" % (ref_grp, ref_grp_read_len, median_site_ratio, mean_site_ratio, SD_site_ration))
    fout_report.write("#\n")

    fout_report.write("END\n")
    fout_report.close() 

    for r in open(fout_report_fname, "rt"):
        print r.rstrip('\n\r')
    print

    if savefig_error:
        print "ERROR: problems saving some figures, see above"
        print "       please check your matplotlib installation"
    print "Full report also available at: %s" % fout_report_fname

def split_bam(in_bam_fname, output_folder, readlen2groups, mapq_th=10, split_flank=5):
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)

    bam_basename = os.path.basename(in_bam_fname).rsplit('.', 1)[0]
    out_fname_pref = "%s/%s" % (output_folder, bam_basename)

    all_recs, notmapped_recs, mapped_recs, lowmapq_recs, split_flank_ignored_recs, used_recs, all_bcs, bins = _read_bam(in_bam_fname, mapq_th, split_flank)

    ### select all reads
    selected_reads = 0
    selected_reads_len_distrib = {}
    for (chrome, strand), by_bin_index in bins.iteritems():
        for bin_index, reads in by_bin_index.iteritems():
            reads_in_bin = len(reads)

            selected_reads += reads_in_bin
            for _, _, _, _, read_len in reads:
                selected_reads_len_distrib[read_len] = selected_reads_len_distrib.get(read_len, 0) + 1
    max_read_len = max(selected_reads_len_distrib.keys())


    ### generate bedGraphs for each group
    g2lens = {}
    for read_len, read_grps in readlen2groups.iteritems():
        for read_grp in read_grps:
            g2lens.setdefault(read_grp, []).append(read_len)

    grp_bed_start = {}
    grp_bed_middle = {}
    grp_bed_end = {}
    for read_grp, valid_read_lens in g2lens.iteritems():
        bed_start = grp_bed_start.setdefault(read_grp, {})
        bed_middle = grp_bed_middle.setdefault(read_grp, {})
        bed_end = grp_bed_end.setdefault(read_grp, {})

        for (chrome, strand), by_bin_index in bins.iteritems():
            tmpd_bed_start = bed_start.setdefault((chrome, strand), {})
            tmpd_bed_middle = bed_middle.setdefault((chrome, strand), {})
            tmpd_bed_end = bed_end.setdefault((chrome, strand), {})

            for bin_index, reads in by_bin_index.iteritems():
                valid_reads = [(xlink_pos, middle_pos, end_pos, bc) for xlink_pos, middle_pos, end_pos, bc, read_len in reads if read_len in valid_read_lens]

                for xlink_pos, middle_pos, end_pos in collapse_bcs(valid_reads, strand):
                    tmpd_bed_start[xlink_pos] = tmpd_bed_start.get(xlink_pos, 0) + 1
                    tmpd_bed_middle[middle_pos] = tmpd_bed_middle.get(middle_pos, 0) + 1
                    tmpd_bed_end[end_pos] = tmpd_bed_end.get(end_pos, 0) + 1
    # bins not needed after this point
    del bins

    ### save bedGraphs
    for read_grp, bed_start in sorted(grp_bed_start.iteritems()):
        iCLIPro.bed.save("%s_%s_start.bed" % (out_fname_pref, read_grp), bed_start, "%s start %s" % (read_grp, bam_basename), db_name="")
    for read_grp, bed_middle in sorted(grp_bed_middle.iteritems()):
        iCLIPro.bed.save("%s_%s_center.bed" % (out_fname_pref, read_grp), bed_middle, "%s center %s" % (read_grp, bam_basename), db_name="")
    for read_grp, bed_end in sorted(grp_bed_end.iteritems()):
        iCLIPro.bed.save("%s_%s_end.bed" % (out_fname_pref, read_grp), bed_end, "%s end %s" % (read_grp, bam_basename), db_name="")

    ### generate report
    fout_report_fname = "%s_report_iCLIPro_bam_splitter.txt" % out_fname_pref
    fout_report = open(fout_report_fname, "wt")

    fout_report.write("# Parameters\n")
    fout_report.write("#   in_bam_fname: %s\n" % in_bam_fname)
    fout_report.write("#   output_folder: %s\n" % output_folder)
    fout_report.write("#   mapq_th: %s\n" % mapq_th)
    fout_report.write("#   split_flank: %s\n" % split_flank)
    fout_report.write("\n")

    fout_report.write("# Read filtering stats\n")
    fout_report.write("#   total number of records: %s\n" % all_recs)
    fout_report.write("#   not mapped: %s (%0.2f%% total)\n" % (notmapped_recs, 100.0*notmapped_recs/all_recs))
    fout_report.write("#   mapped but ignored because MAPQ < %s: %s (%0.2f%% total)\n" % (mapq_th, lowmapq_recs, 100.0*lowmapq_recs/all_recs))
    fout_report.write("#   mapped and with MAPQ >= %s: %s (%0.2f%% total)\n" % (mapq_th, mapped_recs, 100.0*mapped_recs/all_recs))
    fout_report.write("#   additionally ignored because gap in first/last %s nt: %s (%0.2f%% total)\n" % (split_flank, split_flank_ignored_recs, 100.0*split_flank_ignored_recs/all_recs))
    fout_report.write("#   used reads: %s (%0.2f%% total)\n" % (used_recs, 100.0*used_recs/all_recs))
    show_only = 30
    all_bcs = sorted([(cn, bc) for bc, cn in all_bcs.iteritems()], reverse=True)
    all_bcs = [(bc, cn) for cn, bc in all_bcs]
    if len(all_bcs) > show_only:
        fout_report.write("#   distinct random barcodes (showing only %s out of all %s): %s\n" % (show_only, len(all_bcs), all_bcs[:show_only]))
    else:
        fout_report.write("#   distinct random barcodes (%s): %s\n" % (len(all_bcs), all_bcs))
    fout_report.write("\n") 

    fout_report.write("# Used reads\n")
    fout_report.write("#   reads: %s\n" % selected_reads)
    fout_report.write("#   reads length distribution and grouping:\n")
    fout_report.write("#\tread length\tnumber of reads\t% selected reads\tgroup(s)\n")
    notgrouped_reads = 0
    notgrouped_lengths = []
    for read_len, cn in sorted(selected_reads_len_distrib.iteritems()):
        if read_len not in readlen2groups:
            notgrouped_lengths.append(read_len)
            notgrouped_reads += cn
            grp_names = 'notused'
        else:
            grp_names = readlen2groups[read_len]
        fout_report.write("#\t%s\t%s\t%0.2f%%\t%s\n" % (read_len, cn, 100.0*cn/selected_reads, grp_names))
    if notgrouped_reads:
        fout_report.write("WARNING: %s reads (%0.2f%% selected) of lengths %s are not being used because not assigned to any group.\n" % (notgrouped_reads, 100.0*notgrouped_reads/selected_reads, notgrouped_lengths))
    fout_report.write("\n")

    fout_report.write("# Generated bedGraphs\n")
    fout_report.write("#\tgroup\torigin\ttotal cross-linked sites\ttotal cDNA\tfile name\n")
    for read_grp in sorted(g2lens.iterkeys()):
        bed_start = grp_bed_start.setdefault(read_grp, {})
        start_fname = "%s_%s_start.bed" % (out_fname_pref, read_grp)
        start_num_sites, start_tot_cDNA = iCLIPro.bed.stats(bed_start)

        bed_middle = grp_bed_middle.setdefault(read_grp, {})
        middle_fname = "%s_%s_center.bed" % (out_fname_pref, read_grp)
        middle_num_sites, middle_tot_cDNA = iCLIPro.bed.stats(bed_middle)

        bed_end = grp_bed_end.setdefault(read_grp, {})
        end_fname = "%s_%s_end.bed" % (out_fname_pref, read_grp)
        end_num_sites, end_tot_cDNA = iCLIPro.bed.stats(bed_end)

        fout_report.write("#\t%s\t%s\t%s\t%s\t%s\n" % (read_grp, 'start', start_num_sites, start_tot_cDNA, start_fname))
        fout_report.write("#\t%s\t%s\t%s\t%s\t%s\n" % (read_grp, 'center', middle_num_sites, middle_tot_cDNA, middle_fname))
        fout_report.write("#\t%s\t%s\t%s\t%s\t%s\n" % (read_grp, 'end', end_num_sites, end_tot_cDNA, end_fname))
    fout_report.write("\n")

    fout_report.write("END\n")
    fout_report.close() 

    for r in open(fout_report_fname, "rt"):
        print r.rstrip('\n\r')
    print

    print "Full report also available at: %s" % fout_report_fname

