#!/bin/bash

perl -e '$count=0; $len=0; while(<>) {s/\r?\n//; s/\t/ /g; if (s/^>//) { if ($. != 1) {print "\n"} s/|$/\t/; $count++; $_ .= "\t";} else {s/ //g; $len += length($_)} print  $_;} print "\n"; warn "\nConverted $count FASTA records in $. lines to tabular format\nTotal sequence length: $len\n\n";'  $1 > TAB.fasta

cut -f 2,3 TAB.fasta > ${1}.TAB
rm TAB.fasta
