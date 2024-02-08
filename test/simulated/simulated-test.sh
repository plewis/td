#!/bin/bash

DIR="$( cd "$( dirname "$0" )" && pwd )"
cd $DIR

td --reffile ref-tree.tre --treefile other-tree.tre --reftree 1 --skip 0 --outfile dists.txt --debug


