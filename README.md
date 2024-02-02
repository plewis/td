# td

Calculates the Kuhner-Felsenstein distance between trees:

M Kuhner and J Felsenstein. 1994. A simulation comparison of phylogeny
algorithms under equal and unequal evolutionary rates. Molecular Biology
and Evolution 11(3):459-468. 
[DOI:10.1093/oxfordjournals.molbev.a040126](https://doi.org/10.1093/oxfordjournals.molbev.a040126)

It assumes you want tree distances between some reference tree and at 
least one other tree.

If you want distances between the first tree in trees.tre and all the others
Note that reftree assumes trees are indexed starting with 1 (not 0)

./td --treefile trees.tre --reftree 1

If you want distances between the first tree in true.tre and all trees in
trees.tre, use this command line. That is, if reffile is specified, then
reftree indexes a file in reffile; otherwise, reftree indexes a tree in
treefile.

./td --treefile trees.tre --reffile true.tre --reftree 1
