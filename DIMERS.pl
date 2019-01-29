#!/usr/bin/perl

use strict;
use warnings;

use autodie;

my $HOME_DIR = "/scratch/asl47/PDB";

open (my $f_homomers, "<", "${HOME_DIR}/dimers.txt");

my @structures = <$f_homomers>;
shift @structures;
while(my $line = shift @structures) {
    chomp $line;
    last unless $line =~ /\S/;
    {
      local @ARGV = split(' ', $line, 3);
      splice(@ARGV,1,1);    
      do './dimer_analysis.pl';
    }

}


