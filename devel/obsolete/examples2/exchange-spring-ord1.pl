#!/usr/bin/perl -w

# parse exchange-spring-ord1.output, make gnuplot

while(my $line=<>)
  {
    $line=~/M_Fe/ or next;
    my($pos,$dy,$fe)=(0,0,0);

    $line=~/^(\S+)/ and $pos=$1;
    $line=~/M_Dy(?:[^,]*,){2} (\S+),/ and $dy=$1;
    $line=~/M_Fe(?:[^,]*,){2} (\S+),/ and $fe=$1;

    print "$pos $fe $dy\n";
  }
