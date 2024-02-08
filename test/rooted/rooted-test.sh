#!/bin/bash

DIR="$( cd "$( dirname "$0" )" && pwd )"
cd $DIR

td --treefile rooted.nex --reftree 1 --skip 0 --outfile dists.txt --debug


