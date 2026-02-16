#!/bin/bash

DIR="$( cd "$( dirname "$0" )" && pwd )"
cd $DIR

op --reftree unrooted.nex --refskip 0 --treefile unrooted.nex --skip 1 --prefix dists --noisy


