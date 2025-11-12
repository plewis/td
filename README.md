# op

Calculates the Billera-Holmes-Vogtmann geodesic distance 
(using the Owens-Provan algorithm) distance between trees.

If you want geodesic distances between the first tree in 
trees.tre and all the others,

op --treefile trees.tre --refdist --noisy

Including --noisy results in a lot of output, so you may 
want to only use this when trees.tre contains only two trees.

If you want geodesic distances between every pair of trees,

op --treefile trees.tre --pairwise

If you want to save the tree that is a fraction 0.7 of
the distance between the first and second tree in trees.tre,

op --treefile trees.tre --lambda 0.7

If you want to save trees in a format that can be read directly
by GTP (Java program by Owen and Provan),

op --treefile trees.tre --saveforgtp

If you want to calculate the Frechet mean tree,
specify a file name prefix (e.g. "meantree") and the parameters that
determine the stopping criterion:
* k is the maximum number of iterations
* n  is the number of recent iterations to compare
* epsilon is the maximum amount by which every pairwise comparison 
of recent iterations must differ in order to stop. 

The following example uses the default values for k (1000000), 
n (10), and epsilon (0.00001), which are similar to those 
used by Brown and Owen (2020):

op --treefile trees.tre --frechetmean --prefix meantree

You can also specify k, n, or epsilon explicitly:

op --treefile trees.tre --prefix meantree --frechetmean --frechet-e 0.001 --frechet-n 20 --frechet-k 10000

The resulting output file will, in this case, be named "meantree.R" 
and will contain comments (lines beginning with #) specifying the 
mean tree newick description. tree length, variance, and the 
95% HPD interval lower and upper bounds. The file can also 
be executed in R to plot the kernel density of tree distances 
with 95% HPD shaded.

Other options can be seen by asking for help:

./op --help

Literature Cited

LJ Billera, SP Holmes, and K Vogtmann. 2001. Geometry of the space of phylogenetic trees.
Advances in Applied Mathematics 27:733-767.
[DOI:10.1006/aama.2001.0759](https://doi.org/10.1006/aama.2001.0759)

DG Brown and M Owen. 2020. Mean and variance of phylogenetic trees. Systematic Biology 69:139-154. [DOI:10.1093/sysbio/syz041](https://doi.org/10.1093/sysbio/syz041)

M Owens and JS Provan. 2011. A fast algorithm for computing geodesic distances
in tree space. IEEE/ACM Transactions on Computational Biology and Bioinformatics
8:2-12. [DOI:10.1109/TCBB.2010.3](https://doi.org/10.1109/TCBB.2010.3)

