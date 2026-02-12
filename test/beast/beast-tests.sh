#!/bin/bash

DIR="$( cd "$( dirname "$0" )" && pwd )"
cd $DIR

op --reftree true-species-tree.tre --refrooted yes --treefile species.trees  --rooted yes --prefix ref
op --treefile species.trees --rooted yes --frechetmean --prefix mean
