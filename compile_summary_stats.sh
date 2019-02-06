#!/bin/bash

for f in */*oneline_summary_gen*; do cat $f; echo; done > summary_lines.csv
cp selcor*migr002/rep_1_*summary_header.csv temp_header.txt
echo >>temp_header.txt
cat temp_header.txt summary_lines.csv > summary_stats.csv
rm summary_lines.csv
rm temp_header.txt