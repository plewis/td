#!/bin/bash

DIR="$( cd "$( dirname "$0" )" && pwd )"
cd $DIR

td --reffile true-species-tree.tre --treefile species.trees --reftree 1 --skip 2 --outfile separate.txt --quiet
td --treefile species.trees --reftree 1 --skip 1 --outfile combined.txt --quiet


