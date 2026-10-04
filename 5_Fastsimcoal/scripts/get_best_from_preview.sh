#!/usr/bin/env bash
# Parse an easySFS preview file (easySFS.py ... --preview -a > pop1_pop2.preview)
# and print, for each population, the projection that maximizes the number of
# segregating sites: population names, then "<projection> <n_sites>" per population.
#
# Usage: scripts/get_best_from_preview.sh pop1_pop2.preview
grep -B 1 "(2" $1 | grep -v '(2'| grep -v "\-\-"
grep "(2," $1| head -1| tr '\t' '\n'| sed 's/(//'| sed 's/)//'| sed 's/ //'| sed '/^$/d'| tr ',' '\t'| sort -n -k2,2|tail -1
grep "(2," $1| tail -1| tr '\t' '\n'| sed 's/(//'| sed 's/)//'| sed 's/ //'| sed '/^$/d'| tr ',' '\t'| sort -n -k2,2| tail -1
