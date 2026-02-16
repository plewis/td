#!/bin/bash

DIR="$( cd "$( dirname "$0" )" && pwd )"
cd $DIR

#td --treefile unrooted.nex --reftree 1 --skip 0 --outfile dists.txt --debug
op --reftree unrooted.nex --refskip 0 --refrooted yes --treefile unrooted.nex --skip 1 --rooted yes --prefix dists --noisy


