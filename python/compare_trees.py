"""
Comparison of two trees on different leaf sets

Copyright (c) 2018 NJMerge Developers
Erin K. Molloy <molloy.erin.k@gmail.com>
All rights reserved.

License: 3-Clause BSD,
see https://opensource.org/licenses/BSD-3-Clause
"""
import argparse
import dendropy
from dendropy.calculate.treecompare \
    import false_positives_and_negatives
import os
import sys


def compare_trees(tr1, tr2):
    from dendropy.calculate.treecompare \
        import false_positives_and_negatives

    lb1 = set([l.taxon.label for l in tr1.leaf_nodes()])
    lb2 = set([l.taxon.label for l in tr2.leaf_nodes()])

    com = lb1.intersection(lb2)
    if com != lb1 or com != lb2:
        com = list(com)
        tns = dendropy.TaxonNamespace(com)

        tr1.retain_taxa_with_labels(com)
        tr1.migrate_taxon_namespace(tns)

        tr2.retain_taxa_with_labels(com)
        tr2.migrate_taxon_namespace(tns)
    com = list(com)

    tr1.update_bipartitions()
    tr2.update_bipartitions()

    nl = len(com)
    ei1 = len(tr1.internal_edges(exclude_seed_edge=True))
    ei2 = len(tr2.internal_edges(exclude_seed_edge=True))

    [fp, fn] = false_positives_and_negatives(tr1, tr2)
    rf = float(fp + fn) / (ei1 + ei2)

    return(nl, ei1, ei2, fp, fn, rf)


def main(args):
    tax = dendropy.TaxonNamespace()
    tr1 = dendropy.Tree.get(path=args.tree1,
                            schema='newick',
                            rooting='force-unrooted',
                            taxon_namespace=tax)
    tr2 = dendropy.Tree.get(path=args.tree2,
                            schema='newick',
                            rooting='force-unrooted',
                            taxon_namespace=tax)

    tr1.collapse_basal_bifurcation(set_as_unrooted_tree=True)
    tr2.collapse_basal_bifurcation(set_as_unrooted_tree=True)

    [nl, ei1, ei2, fp, fn, rf] = compare_trees(tr1, tr2)
    print("RF distance on %d shared leaves: %d" % (nl, fp + fn))

    # sys.stdout.write('%d %d %d %d %d %f' % (nl, ei1, ei2, fp, fn, rf))
    # sys.stdout.flush()
    # os._exit(0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compares two trees")
    parser.add_argument("-t1", "--tree1", type=str,  required=True,
                        help="File containing newick string for tree 1")
    parser.add_argument("-t2", "--tree2", type=str, required=True,
                        help="File containing newick string for tree 2")
    main(parser.parse_args())
