#!/bin/bash

# generate regression config file

while read p;

do

dtype=$(echo $p| sed 's/.rds//g' | rev | cut -f2 -d_ | rev)
org=$(echo $p | rev | cut -f1 -d/ | rev | cut -f1 -d_)
taxlevel=$(echo $p | rev | cut -f1 -d/ | rev | cut -f4 -d_)
filepath=$p
algorithm=$(echo $p | rev | cut -f1 -d/ | rev | cut -f2 -d_)
cutoffs=$( echo $p| sed 's/.rds//g' | rev | cut -f1 -d_ | rev)
dataframedescr=$(echo $p | rev | cut -f1 -d/ | rev | cut -f3 -d_)

echo -e "$dtype\t$org\t$taxlevel\t$filepath\t$algorithm\t$cutoffs\t$dataframedescr"

done<all_file_locs > regression_config_tax


