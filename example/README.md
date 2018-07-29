
To run NJMerge on the 100-taxon datasets, use the following commands.
```
cd 100-taxon-dataset
python ../../python/njmerge.py \
    -t subset-1-outof-5.tre \
       subset-2-outof-5.tre \
       subset-3-outof-5.tre \
       subset-4-outof-5.tre \
       subset-5-outof-5.tre \
    -m distance.mat \
    -x distance.mat_taxlist \
    -o njmerge.tre
```

Note that 
+ [`distance.mat`](100-taxon-dataset/distance.mat) is a distance matrix in PHYLIP format on the full taxon set
+ [`distance.mat_taxlist`](100-taxon-dataset/distance.mat_taxlist) is the taxon name corresponding to each row in the distance matrix
+ [`subset-1-outof-5.tre`](100-taxon-dataset/subset-1-outof-5.tre) is the first of five constraint trees; all constraint trees need to be on *disjoint* subsets of taxa, i.e., taxa cannot appear in more than one subset tree
+ njmerge.tre is the output file

Now check that the output of NJMerge obeys the constraint trees, using the following commands:
```
../../python/compare_trees.py -t1 subset-1-outof-5.tre -t2 njmerge.tre
../../python/compare_trees.py -t1 subset-2-outof-5.tre -t2 njmerge.tre
../../python/compare_trees.py -t1 subset-3-outof-5.tre -t2 njmerge.tre
../../python/compare_trees.py -t1 subset-4-outof-5.tre -t2 njmerge.tre
../../python/compare_trees.py -t1 subset-5-outof-5.tre -t2 njmerge.tre
```
