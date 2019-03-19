NJMerge
=======
NJMerge is a tool for scaling phylogengy estimation methods to large datasets. NJMerge can be used in a divide-and-conquer framework as follows: 1) divide the species set into disjoint subsets, 2) construct trees on each subset, and 3) combine the subset trees using an associated distance matrix (on the full species set). NJMerge has been sucessfully tested in the context of species tree estimation [(Molloy and Warnow, 2018)](https://doi.org/10.1007/978-3-030-00834-5_15).

REQUIREMENTS
------------
+ Python 2.7 or later
+ [DendroPy](https://www.dendropy.org) 4.3.0

EXAMPLE
-------
To get started using NJMerge, work through [this example](example/README.md).

CITATION
--------
```
@incollection{MolloyWarnow2018,
    author={Erin Molloy and Tandy Warnow},
    editor={Mathieu Blanchette and A{\''{i}}da Ouangraoua},
    title={{NJMerge: A generic technique for scaling phylogeny estimation methods and its application to species trees}},
    booktitle="Comparative Genomics. RECOMB-CG 2018. Lecture Notes in Computer Science},
    volume={11183},
    year={2018},
    publisher={Springer International Publishing},
    address={Cham}
    doi={10.1007/978-3-030-00834-5\_15}
}
```

LICENSE
-------
Please see [LICENSE.txt](LICENSE.txt) for details.
