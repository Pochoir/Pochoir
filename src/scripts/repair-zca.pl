#!/usr/bin/perl -w

# Fix the ZCA metadata in a Cilk assembly listing by removing the old comdat
# system (.gnu.linkonce).

# usage: ./repair-zca.pl foo.s > foo-fixed.s

while (<>) {
    $_ =~ s/\.section \.gnu\.linkonce\.d\.\.itt_notify_tab([^,]*),/.section .itt_notify_tab,/;
    print "$_";
}

