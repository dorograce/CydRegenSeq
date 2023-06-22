#!/usr/bin/perl

# simple script to parse tabbed BLAST output and print query, subject, eval
# report bugs to Joseph Ryan <jfryan@ufl.edu>

use strict;
use warnings;

MAIN: {
    my $tab = $ARGV[0] or die "usage: $0 TABBED_BLAST\n";
    open IN, $tab or die "cannot open $tab:$!";
    my %hits = ();
    while (my $line = <IN>) {
        my @f = split /\t/, $line;
        next if ($hits{$f[0]});
        $hits{$f[0]} = $f[1];
        print "$f[0]\t$f[1]\t$f[10]\n";
    }
}

