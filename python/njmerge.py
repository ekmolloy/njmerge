"""
This file is a python prototype of NJMerge from the paper:

Molloy, E.K., Warnow, T. (2018). NJMerge: A generic technique for
    scaling phylogeny estimation methods and its application to
    species trees. In Blanchette, M. and Ouangraoua, A., editors,
    Comparative Genomics. RECOMB-CG 2018. Lecture Notes in Computer
    Science. Springer International Publishing, Cham. To appear.

Copyright (c) 2018 NJMerge Developers
Erin K. Molloy <molloy.erin.k@gmail.com>
All rights reserved.

License: 3-Clause BSD,
see https://opensource.org/licenses/BSD-3-Clause

Parts of the merge_trees_via_nj function is based on the nj_tree
function in DendroPy (c) 2010-2015 Jeet Sukumaran and Mark T. Holder.

NJMerge, like Dendropy, is licensed under the 3-Clause BSD License.
"""
import argparse
from copy import deepcopy
import dendropy
from dendropy.calculate.treecompare import false_positives_and_negatives
import numpy
import os
from string import maketrans
import subprocess
import sys
import time

sys.setrecursionlimit(10000)


def matrix_to_dendropy_pdm(dmat, taxa):
    """Read FastME distance matrix into a dendropy PDM object

    Parameters
    ----------
    dmat : str
        paup* distance matrix file name

    Returns
    -------
    pdm : dendropy phylogenetic distance matrix object

    """
    pdm = dendropy.PhylogeneticDistanceMatrix()
    pdm.taxon_namespace = dendropy.TaxonNamespace()
    pdm._mapped_taxa = set()

    for i, si in enumerate(taxa):
        for j, sj in enumerate(taxa):
            dij = dmat[i, j]

            xi = pdm.taxon_namespace.get_taxon(si)
            if not xi:
                xi = dendropy.Taxon(si)
                pdm.taxon_namespace.add_taxon(xi)
                pdm._mapped_taxa.add(xi)
                pdm._taxon_phylogenetic_distances[xi] = {}

            xj = pdm.taxon_namespace.get_taxon(sj)
            if not xj:
                xj = dendropy.Taxon(sj)
                pdm.taxon_namespace.add_taxon(xj)
                pdm._mapped_taxa.add(xj)
                pdm._taxon_phylogenetic_distances[xj] = {}

            dij = float(dij)
            pdm._taxon_phylogenetic_distances[xi][xj] = dij
    return pdm


def get_leaf_list(subtree):
    """Return list of leaf labels

    Parameters
    ----------
    subtree : dendropy tree or node object

    Returns
    -------
    list of str

    """
    return [l.taxon.label for l in subtree.leaf_nodes()]


def get_leaf_set(subtree):
    """Return set of leaf labels

    Parameters
    ----------
    subtree : dendropy tree or node object

    Returns
    -------
    set of str

    """
    return set([l.taxon.label for l in subtree.leaf_nodes()])


def are_two_trees_incompatible(tree1, tree2):
    """Check if two unrooted trees are equivalent on their shared taxon set

    Parameters
    ----------
    tree1 : dendropy tree object
    tree2 : dendropy tree object

    Returns
    -------
    violates : bool
        True, if trees are NOT compatible
        False, if trees are compatible

    """
    leaves1 = get_leaf_set(tree1)
    leaves2 = get_leaf_set(tree2)
    shared = list(leaves1.intersection(leaves2))

    taxa = dendropy.TaxonNamespace(shared)  # CRITICAL!!!

    # No topological information
    if len(shared) < 4:
        return False

    # Move trees onto shared leaf set
    tree1.retain_taxa_with_labels(shared)
    tree1.migrate_taxon_namespace(taxa)
    tree1.is_rooted = False
    tree1.collapse_basal_bifurcation()
    tree1.update_bipartitions()

    tree2.retain_taxa_with_labels(shared)
    tree2.migrate_taxon_namespace(taxa)
    tree2.is_rooted = False
    tree2.collapse_basal_bifurcation()
    tree2.update_bipartitions()

    # Check for compatibility
    [fp, fn] = false_positives_and_negatives(tree1, tree2)
    if fp > 0 or fn > 0:
        return True
    else:
        return False


def map_splits_to_nodes(tree):
    """Map splits (encoded as integers) to nodes in tree

    NOTE: dendropy had/had the same functionality built-in
          but I never got it to work right...

    Parameters
    ----------
    tree : dendropy tree object

    Returns
    -------
    split_to_edge_map : dictionary
        keys are splits encoded as integers (read below!)
        values are nodes in dendropy tree object

    """
    tns = tree.taxon_namespace
    ntx = len(tns)

    # Update bipartitions and grab bipartition to edge map
    tree.update_bipartitions()
    bipartition_to_edge_map = tree.bipartition_edge_map

    # Bipartitions become splits (integers) and the edges become nodes
    split_to_node_map = {}
    for b in list(bipartition_to_edge_map):
        node = bipartition_to_edge_map[b].head_node

        bitmask1 = b.split_bitmask
        split_to_node_map[bitmask1] = node

        # Not sure if this is necessary, but I re-root the trees a lot.
        # So better safe than sorry until further testing...
        bitstring1 = b.split_as_bitstring()
        bitstring2 = bitstring1.translate(maketrans("10", "01"))
        bitmask2 = int(bitstring2, 2)
        split_to_node_map[bitmask2] = node

    return split_to_node_map


def get_node_from_clade(tree, split_to_node_map, clade):
    """Returns the parent node of the clade

    Parameters
    ----------
    tree : dendropy tree object
    split_to_node_map : dictionary
       keys are splits encoded as integers
        values are nodes in dendropy tree object
    clade : list of str
        taxon labels

    Returns
    -------
    dendropy node object or None

    """
    leaves = get_leaf_set(tree)

    # Check if the clade is the whole tree!
    if leaves == clade:
        return split_to_node_map[0]

    # Check if the clade contains leaves not in the tree itself
    if len(leaves.intersection(clade)) < len(clade):
        return None

    # Encode labels as split (integer) and return the node or none
    split = tree.taxon_namespace.taxa_bitmask(labels=clade)
    if split in split_to_node_map:
        return split_to_node_map[split]
    else:
        return None


def extract_nodes_from_split(tree, node, clade):
    """If you re-root tree at node, then one side of the root
    will correspond to clade. Return the re-rooted tree and the
    parent node of the clade.

    Parameters
    ----------
    tree : dendropy tree object
    node : dendropy node object
        node at the split
    clade : set of str
        taxon labels

    Returns
    -------
    tree : dendropy tree object
        re-rooted at node
    found : dendropy node
        node corresponding to the clade

    """
    # Re-root tree as rooted, i.e., (clade, rest);
    if node is not tree.seed_node:
        tree.reroot_at_edge(node.edge)

    bipartition = tree.seed_node.child_nodes()

    if len(bipartition) != 2:
        raise Exception("Tree is not rooted!\n")

    [node1, node2] = bipartition
    leaves1 = get_leaf_set(node1)
    leaves2 = get_leaf_set(node2)

    if leaves1 == clade:
        found = node1
    elif leaves2 == clade:
        found = node2
    else:
        raise Exception("Cannot handle this bipartition!\n")

    return [tree, found]


def join_nodes_in_both_trees(tree1, nodeAinT1, cladeA,
                             tree2, nodeBinT2, cladeB, test=False):
    """Join clade A and clade B in both trees

    1. Re-root tree 1 as (A, X); and extract A and X
    2. Re-root tree 2 as (B, X'); and extract B and X'
    3. Build new tree 1 as (A, B, X);
    4. Build new tree 2 as (A, B, X');
    5. Return tree 1 and tree 2

    Parameters
    ----------
    tree1 : dendropy tree object
    nodeAinT1 : dendropy node object
    cladeA : list of str
        taxon labels
    tree2 : dendropy node object
    nodeBinT2 : dendropy node object
    cladeB : list of str
        taxon labels

    Returns
    -------
    tree1 : dendropy tree object
    tree2 : dendropy tree object

    """
    cladeA = set(cladeA)
    cladeB = set(cladeB)
    leaves1 = get_leaf_set(tree1)
    leaves2 = get_leaf_set(tree2)

    cladeAisT1 = leaves1 == cladeA
    cladeBisT2 = leaves2 == cladeB

    # Handle adding all of tree1 into tree 2 and vice versa!!
    if cladeAisT1 and cladeBisT2:
        # Join trees 1 and 2 together, i.e.,
        # tree 1 will equal tree 2
        # print("Joining tree1 and tree2...")
        if test:
            return [None, None]
        root = dendropy.Node()
        root.add_child(nodeAinT1)
        root.add_child(nodeBinT2)
        tree1 = dendropy.Tree(seed_node=root)
        tree1.is_rooted = True
        tree2 = None
    elif cladeAisT1:
        # Add all of tree 1 into tree 2, i.e.,
        # tree 1 will be contained in tree 2
        # print("Add all of tree 1 into tree 2")
        if test:
            return [None, None]
        [tree2, nodeBinT2] = extract_nodes_from_split(tree2, nodeBinT2,
                                                      cladeB)
        root = dendropy.Node()
        root.add_child(nodeAinT1)
        root.add_child(tree2.seed_node)
        tree1 = dendropy.Tree(seed_node=root)
        tree1.is_rooted = True
        tree2 = None
    elif cladeBisT2:
        # Add all of tree 2 into tree 1, i.e.,
        # tree will be contained in tree 1
        # print("Add all of tree 2 into tree 1")
        if test:
            return [None, None]
        [tree1, nodeAinT1] = extract_nodes_from_split(tree1, nodeAinT1,
                                                      cladeA)
        root = dendropy.Node()
        root.add_child(tree1.seed_node)
        root.add_child(nodeBinT2)
        tree1 = dendropy.Tree(seed_node=root)
        tree1.is_rooted = True
        tree2 = None
    else:
        # Make the join, i.e.,
        # tree 1 and tree 2 will have a shared leaf
        # print("Making join...")
        [tree1, nodeAinT1] = extract_nodes_from_split(tree1, nodeAinT1,
                                                      cladeA)
        [tree2, nodeBinT2] = extract_nodes_from_split(tree2, nodeBinT2,
                                                      cladeB)

        root1 = dendropy.Node()
        root1.add_child(tree1.seed_node)
        root1.add_child(deepcopy(nodeBinT2))  # TODO: Remove deep copies!
        tree1 = dendropy.Tree(seed_node=root1)
        tree1.is_rooted = True

        root2 = dendropy.Node()
        root2.add_child(tree2.seed_node)
        root2.add_child(deepcopy(nodeAinT1))  # TODO: Remove deep copies!
        tree2 = dendropy.Tree(seed_node=root2)
        tree2.is_rooted = True

    return [tree1, tree2]


def join_nodes_in_one_tree(tree1, nodeAinT1, cladeA, tree2, nodeBinT2,
                           cladeB):
    """Join clade A and clade B in just one tree

    1. Re-root tree 1 as (A,X); and extract (A,X); and A
    2. Re-root tree 2 as (B,Y); and extract (B,Y); and B
    3. Build new tree 2 as (A,(B,Y));
    4. Return tree 1 and tree 2

    Parameters
    ----------
    tree1 : dendropy tree object
    nodeAinT1 : dendropy node object
    cladeA : list of str
        taxon labels below node A
    tree2 : dendropy node object
    nodeBinT2 : dendropy node object
    cladeB : list of str
        taxon labels below node B

    Returns
    -------
    tree1 : dendropy tree object
    tree2 : dendropy tree object

    """
    [tree1, nodeAinT1] = extract_nodes_from_split(tree1, nodeAinT1, cladeA)
    [tree2, nodeBinT2] = extract_nodes_from_split(tree2, nodeBinT2, cladeB)

    root = dendropy.Node()
    root.add_child(deepcopy(nodeAinT1))  # TODO: Remove deep copies!
    root.add_child(tree2.seed_node)
    tree2 = dendropy.Tree(seed_node=root)
    tree2.is_rooted = True

    return [tree1, tree2]


def test_join(trees, leaves, maps, nodeA, nodeB):
    """Test whether joining cladeA and cladeB in one
    or both trees causes the two trees to be incompatible

    Parameters
    ----------
    trees : list of dendropy tree objects
    leaves : list of sets
        leaf set for each constraint tree
    maps : list of dictionaries
        split-to-node map for each constraint tree
    nodeA : dendropy node
    nodeB : dendropy node

    Returns
    -------
    violates : boolean
        True, if trees are NOT compatible
        False, if trees are compatible

    """
    cladeA = get_leaf_set(nodeA)
    cladeB = get_leaf_set(nodeB)
    cladeAB = cladeA.union(cladeB)

    if len(cladeA.intersection(cladeB)) > 0:
        raise Exception("Nodes are not disjoint on their leaf sets!\n")

    nodeAinTs = []
    nodeBinTs = []
    for i, tree in enumerate(trees):
        leaf = leaves[i]

        if cladeA == leaf:
            nodeAinTs.append(nodeA)
        else:
            nodeAinTs.append(get_node_from_clade(tree, maps[i], cladeA))

        if cladeB == leaf:
            nodeBinTs.append(nodeB)
        else:
            nodeBinTs.append(get_node_from_clade(tree, maps[i], cladeB))

    nAinTs = []
    for nodeAinT in nodeAinTs:
        nAinTs.append(nodeAinT is not None)

    if sum(nAinTs) < 1:
        raise Exception("Node A was not found in any tree!\n")

    nBinTs = []
    for nodeBinT in nodeBinTs:
        nBinTs.append(nodeBinT is not None)

    if sum(nBinTs) < 1:
        raise Exception("Node B was not found in any tree!\n")

    violates = False
    for i, tree1 in enumerate(trees[:-1]):
        map1 = maps[i]
        nAinT1 = nAinTs[i]
        nBinT1 = nBinTs[i]
        nodeAinT1 = nodeAinTs[i]
        nodeBinT1 = nodeBinTs[i]
        j = i
        for tree2 in trees[i+1:]:
            j = j + 1
            map2 = maps[j]
            nAinT2 = nAinTs[j]
            nBinT2 = nBinTs[j]
            nodeAinT2 = nodeAinTs[j]
            nodeBinT2 = nodeBinTs[j]
            if (sum([nAinT1, nAinT2]) > 0) and (sum([nBinT1, nBinT2]) > 0):
                if nAinT1 and nAinT2:
                    # nodeA in *both* T1 and T2
                    if nBinT1 and nBinT2:
                        # Case 1: nodeB in *both* T1 and T2
                        # Valid if nodeA and nodeB are siblings in both T1 & T2
                        node1 = get_node_from_clade(tree1, map1, cladeAB)
                        node2 = get_node_from_clade(tree2, map2, cladeAB)
                        if (node1 is None) or (node2 is None):
                            violates = True
                    elif nBinT1:
                        # Case 2: Node B in T1 only
                        # Valid if nodeA and nodeB are siblings in T1
                        node = get_node_from_clade(tree1, map1, cladeAB)
                        if node is None:
                            violates = True
                    elif nBinT2:
                        # Case 3: Node B in T2 only
                        # Valid if nodeA and nodeB are siblings in T2
                        node = get_node_from_clade(tree2, map2, cladeAB)
                        if node is None:
                            violates = True
                    else:
                        raise Exception("Node B not found in either tree!\n")
                elif nAinT1:
                    # nodeA in T1 only
                    if nBinT1 and nBinT2:
                        # Case 4: nodeB in *both* T1 and T2
                        node = get_node_from_clade(tree1, map1, cladeAB)
                        if node is None:
                            violates = True
                    elif nBinT1:
                        # Case 5: Node B in T1 only
                        # Valid if nodeA and nodeB are siblings in T1
                        node = get_node_from_clade(tree1, map1, cladeAB)
                        if node is None:
                            violates = True
                    elif nBinT2:
                        # Case 6: Node B in T2 only
                        # Do join in both trees and test for compatibility
                        # TODO: Remove deep copies!
                        t1 = deepcopy(tree1)
                        t2 = deepcopy(tree2)
                        nA = deepcopy(nodeAinT1)
                        nB = deepcopy(nodeBinT2)
                        [t1, t2] = join_nodes_in_both_trees(t1, nA, cladeA,
                                                            t2, nB, cladeB,
                                                            test=True)
                        if t1 is not None:
                            violates = are_two_trees_incompatible(t1, t2)
                    else:
                        raise Exception("Node B not found in either tree!\n")
                elif nAinT2:
                    # nodeA in T2 only
                    if nBinT1 and nBinT2:
                        # Case 7: nodeB in *both* T1 and T2
                        node = get_node_from_clade(tree2, map2, cladeAB)
                        if node is None:
                            violates = True
                    elif nBinT1:
                        # Case 8 (reverse of Case 6): Node B in T1 only
                        # Do join in both trees and test for compatibility
                        # TODO: Remove deep copies!
                        t1 = deepcopy(tree1)
                        t2 = deepcopy(tree2)
                        nA = deepcopy(nodeAinT2)
                        nB = deepcopy(nodeBinT1)
                        [t1, t2] = join_nodes_in_both_trees(t1, nB, cladeB,
                                                            t2, nA, cladeA,
                                                            test=True)
                        if t1 is not None:
                            violates = are_two_trees_incompatible(t1, t2)
                    elif nBinT2:
                        # Case 9: Node B in T2 only
                        # Only valid if nodeA and nodeB are siblings in T2
                        node = get_node_from_clade(tree2, map2, cladeAB)
                        if node is None:
                            violates = True
                    else:
                        raise Exception("Node B not found in either tree!\n")
                else:
                    raise Exception("Node A not found in either tree!\n")

            if violates:
                return violates

    return violates


def join_nodes(trees, leaves, maps, nodeA, nodeB):
    """Join cladeA and cladeB in one or both trees

    Parameters
    ----------
    trees : list of dendropy tree objects
    leaves : list of sets
        leaves for trees
    maps : list of dictionaries
        split-to-node map for trees
    nodeA : dendropy node
    nodeB : dendropy node

    Returns
    -------
    trees : list of dendropy tree objects
    edits : list of booleans
        whether or not constraint tree was edited
    """
    cladeA = get_leaf_set(nodeA)
    cladeB = get_leaf_set(nodeB)

    if len(cladeA.intersection(cladeB)) > 0:
        raise Exception("Nodes are not disjoint on their leaf sets!\n")

    edits = [False] * len(trees)
    for i, edit in enumerate(edits):
        leaf = leaves[i]

        if cladeA == leaf:
            nodeAinT = nodeA  # TODO: Fix having duplicate constraints
        else:
            nodeAinT = get_node_from_clade(trees[i], maps[i], cladeA)

        if cladeB == leaf:
            nodeBinT = nodeB  # TODO: Fix having duplicate constraints
        else:
            nodeBinT = get_node_from_clade(trees[i], maps[i], cladeB)

        nAinT = nodeAinT is not None
        nBinT = nodeBinT is not None

        if nAinT and nBinT:
            # Node A and node B are both in T, do nothing!
            pass
        elif nAinT:
            # Add node B to T
            edits[i] = True
            root = dendropy.Node()
            root.add_child(deepcopy(nodeB))  # TODO: Remove deep copies!
            if leaves[i] == cladeA:
                nodeAinT.parent_node = None
                root.add_child(nodeAinT)
            else:
                [tree, nodeAinT] = extract_nodes_from_split(trees[i], nodeAinT,
                                                            cladeA)
                root.add_child(tree.seed_node)
            trees[i] = dendropy.Tree(seed_node=root)
            trees[i].is_rooted = True
        elif nBinT:
            # Add node A to T
            edits[i] = True
            root = dendropy.Node()
            root.add_child(deepcopy(nodeA))  # TODO: Remove deep copies!
            if leaves[i] == cladeB:
                nodeBinT.parent_node = None
                root.add_child(nodeBinT)
            else:
                [tree, nodeBinT] = extract_nodes_from_split(trees[i], nodeBinT,
                                                            cladeB)
                root.add_child(tree.seed_node)
            trees[i] = dendropy.Tree(seed_node=root)
            trees[i].is_rooted = True
        else:
            # Neither node A or node B is in T, do nothing!
            pass

    return [trees, edits]


def merge_trees_via_nj(pdm, trees):
    """Return a Neighbor-Joining tree that is compatible
    with the set of constraint trees on disjoint leaf sets

    Parameters
    ----------
    pdm : dendropy phylogenetic distance matrix object
    trees : list of dendropy tree objects
        constraint trees

    Returns
    -------
    tree : dendropy tree object
        constrained NJ tree

    """
    leaves = []
    for tree in trees:
        leaves.append(get_leaf_set(tree))

    # Check trees are on disjoint leaf sets
    for i, li in enumerate(leaves[:-1]):
        for lj in leaves[i+1:]:
            shared = li.intersection(lj)
            if len(shared) != 0:
                raise Exception("Input trees are not on disjoint leaf sets!\n")

    # Check distance matrix and trees have matching leaf sets
    full_leaf_set = set()
    for l in leaves:
        full_leaf_set = full_leaf_set.union(l)
    if full_leaf_set != set([x.label for x in pdm.taxon_namespace]):
        raise Exception("Names in matrix do not match those in trees!\n")

    # Remove some extra nonsense
    for tree in trees:
        # Root trees
        tree.resolve_polytomies(limit=2)
        tree.is_rooted = True

        # Remove branch lengths
        for e in tree.preorder_edge_iter():
            e.length = None

        # Remove bootstrap support
        for n in tree.internal_nodes():
            n.label = None

    # Map splits to nodes
    maps = []
    for tree in trees:
        maps.append(map_splits_to_nodes(tree))

    # Taken from dendropy
    original_dmatrix = pdm._taxon_phylogenetic_distances
    tree_factory = dendropy.Tree
    tree = tree_factory(taxon_namespace=pdm.taxon_namespace)
    tree.is_rooted = False

    # Initialize node pool - taken from dendropy
    node_pool = []
    for t1 in pdm._mapped_taxa:
        nd = tree.node_factory()
        nd.taxon = t1
        nd._nj_distances = {}
        node_pool.append(nd)

    # Initialize factor - taken from dendropy
    n = len(pdm._mapped_taxa)

    # Cache calculations - taken from dendropy
    for nd1 in node_pool:
        nd1._nj_xsub = 0.0
        for nd2 in node_pool:
            if nd1 is nd2:
                continue
            d = original_dmatrix[nd1.taxon][nd2.taxon]
            nd1._nj_distances[nd2] = d
            nd1._nj_xsub += d

    while n > 1:
        print("%d joins to go!" % n)

        # Sort the Q-matrix
        # TODO: Use multi-threading!
        pairs = []
        qvalues = []
        for idx1, nd1 in enumerate(node_pool[:-1]):
            idx2 = idx1 + 1
            for nd2 in node_pool[idx2:]:
                v1 = (n - 2) * nd1._nj_distances[nd2]
                qvalue = v1 - nd1._nj_xsub - nd2._nj_xsub
                pairs.append([idx1, idx2])
                qvalues.append(qvalue)
                idx2 = idx2 + 1

        # Test for constraint violations
        # TODO: Use multi-threading in test_join function!
        nodes_to_join = None
        for idxq in numpy.argsort(qvalues):
            [idx1, idx2] = pairs[idxq]
            nd1 = node_pool[idx1]
            nd2 = node_pool[idx2]
            # Check join does not violate a constraint tree!
            violates = test_join(trees, leaves, maps, nd1, nd2)
            if not violates:
                nodes_to_join = (nd1, nd2)
                break

        if nodes_to_join is None:
            raise Exception("Unable to find valid siblinghood!\n")

        # Nodes to join
        (nd1, nd2) = nodes_to_join

        # Update the constraint trees!
        [trees, edits] = join_nodes(trees, leaves, maps, nd1, nd2)

        if sum(edits) > 0:
            i = 0
            for t, e in zip(trees, edits):
                if e:
                    # Check to see if you can quit early
                    leaves[i] = get_leaf_set(t)
                    if leaves[i] == full_leaf_set:
                        # print("Able to quit early!")
                        return t

                    # Update split-to-node maps
                    maps[i] = map_splits_to_nodes(t)
                i = i + 1

        # Delete duplicate constraint trees, i.e.,
        # constraint trees already on the exact same leaf set
        delete = False
        for i, li in enumerate(leaves[:-1]):
            j = i + 1
            for lj in leaves[j:]:
                # print("i=%d j=%d" % (i,j))
                if li == lj:
                    print("Tree %d equals %d" % (i, j))
                    delete = True
                    del trees[j]
                    del maps[j]
                    del leaves[j]
                    break
                elif li.issubset(lj):
                    print("Tree %d contained in %d" % (i, j))
                    delete = True
                    del trees[i]
                    del maps[i]
                    del leaves[i]
                    break
                elif lj.issubset(li):
                    print("Tree %d contained in %d" % (j, i))
                    delete = True
                    del trees[j]
                    del maps[j]
                    del leaves[j]
                    break
                j = j + 1
        if delete:
            print("%d constraint trees left!" % len(trees))

        # Create the new node - taken from dendropy
        new_node = tree.node_factory()

        # Attach it to the tree - taken from dendropy
        for node_to_join in nodes_to_join:
            new_node.add_child(node_to_join)
            node_pool.remove(node_to_join)

        # Calculate the distances for the new node - taken from dendropy
        new_node._nj_distances = {}
        new_node._nj_xsub = 0.0
        for node in node_pool:
            # actual node-to-node distances
            v1 = 0.0
            for node_to_join in nodes_to_join:
                v1 += node._nj_distances[node_to_join]
            v3 = nodes_to_join[0]._nj_distances[nodes_to_join[1]]
            dist = 0.5 * (v1 - v3)
            new_node._nj_distances[node] = dist
            node._nj_distances[new_node] = dist

            # Adjust/recalculate the values needed for the Q-matrix
            # calculations - taken from dendropy
            new_node._nj_xsub += dist
            node._nj_xsub += dist
            for node_to_join in nodes_to_join:
                node._nj_xsub -= node_to_join._nj_distances[node]

        # Clean up - taken from dendropy
        for node_to_join in nodes_to_join:
            del node_to_join._nj_distances
            del node_to_join._nj_xsub

        # Add the new node to the pool of nodes - taken from dendropy
        node_pool.append(new_node)

        # Adjust count - taken from dendropy
        n -= 1

    # More clean up - taken from dendropy
    tree.seed_node = node_pool[0]
    del tree.seed_node._nj_distances
    del tree.seed_node._nj_xsub
    return tree


def main(args):
    if args.trees is None:
        print("Warning: No constraint trees given, running NJ!")
    else:
        ts = time.time()
        trees = []
        for tree in args.trees:
            trees.append(dendropy.Tree.get(path=tree, schema="newick"))
        te = time.time()
        t = te - ts
        print("Converted read all trees in %f seconds" % t)

    # Read taxa list
    with open(args.taxa_file, "r") as f:
        taxa = [l.replace("\n", "") for l in f]

    # Read distance matrix
    ts = time.time()
    with open(args.matrix_file, "r") as f:
        ntax = int(f.readline())
        dmat = numpy.zeros((ntax, ntax))
        for i, line in enumerate(f):
            row = numpy.array([float(x) for x in line.split()[1:]])
            dmat[i, :] = row
    te = time.time()
    t = te - ts
    print("Read matrix into numpy array in %f seconds" % t)

    # TODO: Make this faster!
    ts = time.time()
    pdm = matrix_to_dendropy_pdm(dmat, taxa)
    te = time.time()
    t = te - ts
    print("Converted numpy matrix into a dictionary in %f seconds" % t)

    # Merge trees
    if args.trees is None:
        ts = time.time()
        merged_tree = pdm.nj_tree()
        merged_tree.is_rooted = True
        te = time.time()
        t = te - ts
        print("Built NJ tree in %f seconds" % t)
    else:
        ts = time.time()
        merged_tree = merge_trees_via_nj(pdm, trees)
        te = time.time()
        t = te - ts
        print("Built NJMerge tree in %f seconds" % t)

    with open(args.output_file, "w") as f:
        f.write(merged_tree.as_string(schema="newick")[5:])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="NJMerge")

    parser.add_argument("-t", "--trees", type=str, nargs="+",
                        help="Contraint tree files", required=False)

    parser.add_argument("-m", "--matrix_file", type=str,
                        help="Matrix file", required=True)

    parser.add_argument("-x", "--taxa_file", type=str,
                        help="Taxa list for matrix rows", required=True)

    parser.add_argument("-o", "--output_file", type=str,
                        help="Output file name", required=True)

    main(parser.parse_args())
