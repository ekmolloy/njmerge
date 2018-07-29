
To run NJMerge on the 100-taxon datasets, use the following commands:
```
cd 100-taxon-dataset
```

```
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

Now you can check that the output of NJMerge obeys the constraint trees, using the following commands:
```
../../python/compare_trees.py -t1 subset-1-outof-5.tre -t2 njmerge.tre
../../python/compare_trees.py -t1 subset-2-outof-5.tre -t2 njmerge.tre
../../python/compare_trees.py -t1 subset-3-outof-5.tre -t2 njmerge.tre
../../python/compare_trees.py -t1 subset-4-outof-5.tre -t2 njmerge.tre
../../python/compare_trees.py -t1 subset-5-outof-5.tre -t2 njmerge.tre
```

You can also check that the output of NJMerge is different thatn the output of NJ:
