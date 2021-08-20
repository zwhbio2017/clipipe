def parse_readlen_groups(s):
    grp_names = set()
    retd = {}
    for o in s.split(','):
        o = o.strip()
        grp_name, grp_len = o.split(':')
        if grp_name in grp_names:
            print "Error parsing read len groups definition: same group name is used twice:", grpname
            return
        ns = grp_len.split('-')
        if len(ns) == 1:
            try:
                n = int(ns[0])
            except:
                print "Error parsing read len groups definition: not an integer:", grp_name, grp_len
                return
            retd.setdefault(n, []).append(grp_name)
        elif len(ns) == 2:
            n1, n2 = ns
            try:
                n1 = int(n1)
                n2 = int(n2)
            except:
                print "Error parsing read len groups definition: not an integer:", grp_name, grp_len
                return
            for i in range(n1, n2+1):
                retd.setdefault(i, []).append(grp_name)
        else:
            print "Error parsing read len groups definition: wrong format for length span:", grp_name, grp_len
            return
    return retd

def parse_group_comparisons(s, readlen2groups):
    # which group names are given in the group read len definition
    retl = []
    valid_grp_names = set()
    for tmps in readlen2groups.itervalues():
        for grp_name in tmps:
            valid_grp_names.add(grp_name)
    for o in s.split(','):
        grp_pair = o.strip().split('-')
        if len(grp_pair) <> 2:
            print "Error parsing group comparisons definition: wrong comparison:", o
            return
        g1, g2 = grp_pair
        g1 = g1.strip()
        g2 = g2.strip()
        if g1 not in valid_grp_names:
            print "Error parsing group comparisons definition: group '%s' used in comparison ('%s') is not defined" % (g1, o)
            return
        if g2 not in valid_grp_names:
            print "Error parsing group comparisons definition: group '%s' used in comparison ('%s') is not defined" % (g2, o)
            return
        retl.append((g1, g2))
    return retl 

def parse_neighbor_clustering(s):
    retl = []
    for o in s.split(','):
        o = o.strip()
        try:
            n = int(o)
            retl.append(n)
        except:
            print "Error parsing list of neighbor clustering max distance, value is not integer: '%s'" % (o)
            return
    return sorted(retl)

