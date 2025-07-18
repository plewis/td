# op

Calculates either the Billera-Holmes-Vogtmann geodesic distance 
(using the Owens-Provan algorithm) or Kuhner-Felsenstein distance between trees:

It assumes you want tree distances between some reference tree and at 
least one other tree.

If you want geodesic distances (the default) between the first tree in 
trees.tre and all the others,

./td --treefile trees.tre

If you want KF distances between the first tree in trees.tre and all the others,

./td --treefile trees.tre --dist kf

If you want geodesic distances (the default) between the first tree in 
trees.tre and all the others, and you want to see all the details,

./td --treefile trees.tre --dist geodesic --quiet no

Note: in the command above, it is not necessary to specify "--dist geodesic"
because geodesic distances are the default. Also note that --quiet is, by
default, true.

Literature Cited

LJ Billera, SP Holmes, and K Vogtmann. 2001. Geometry of the space of phylogenetic trees.
Advances in Applied Mathematics 27:733-767.
[DOI:10.1006/aama.2001.0759](https://doi.org/10.1006/aama.2001.0759)

M Kuhner and J Felsenstein. 1994. A simulation comparison of phylogeny
algorithms under equal and unequal evolutionary rates. Molecular Biology
and Evolution 11(3):459-468. 
[DOI:10.1093/oxfordjournals.molbev.a040126](https://doi.org/10.1093/oxfordjournals.molbev.a040126)

M Owens and JS Provan. 2011. A fast algorithm for computing geodesic distances
in tree space. IEEE/ACM Transactions on Computational Biology and Bioinformatics
8:2-12. [DOI:10.1109/TCBB.2010.3](https://doi.org/10.1109/TCBB.2010.3)

