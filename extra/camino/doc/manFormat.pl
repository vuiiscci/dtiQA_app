#!/usr/bin/perl -w

#
# (PAC) - formats man page text to a fixed column width

use strict;
use Text::Wrap;

$Text::Wrap::columns = 90; 

if ($#ARGV < 0) {
   print "    manFormat.pl <man page> \n";
   exit 1;
}

my $file = $ARGV[0];

# Read lines into array
open TEXT , "<$file"; 


while (my $line = <TEXT>) {

  my $paragraph = "";

  # skip header lines
  if ($line =~ m/^\.(SH|TP|TH|B)/ || $line =~ m/^ /) {
    $paragraph = $line;
  } 
  else {

    $paragraph = $paragraph . $line;
   
    PARA: while ($line = <TEXT>) {
        last PARA if ($line =~ m/^\.(SH|TP|TH|B)/ || $line =~ m/^ /);

        $line =~ s/[ ]{2,}/ /g;
        $line =~ s/ +$//;

        $paragraph = $paragraph . $line;
    }

    # Get rid of newlines placed in the paragraph bodies
    $paragraph =~ s/([^\n])\n([^\n])/$1 $2/g; 

    $paragraph = wrap("", "", $paragraph);

    # Add next header line 
    if ($line) {
      $paragraph = $paragraph . "\n$line";
    }

  }
    
  print $paragraph;

}

close TEXT;
