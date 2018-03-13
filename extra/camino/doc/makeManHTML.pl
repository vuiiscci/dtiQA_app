#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);

$Bin =~ m|(.*)/doc|;

my $caminoDir = $1;

`mkdir -p ${caminoDir}/doc/manhtml/man1`;

my @manPages = `ls ${caminoDir}/man/man1`; 


my $header = qq{
<HTML><HEAD><TITLE>Camino man pages</TITLE>
</HEAD><BODY>
<H1>Camino man page index</H1>
};

open INDEX, ">${caminoDir}/doc/manhtml/index.html";  

print INDEX $header;

foreach my $manFile (@manPages) {

  chomp $manFile;

  `man2html -r -p ${caminoDir}/man/man1/$manFile > ${caminoDir}/doc/manhtml/man1/${manFile}.html`;

  $manFile =~ m/(.*)\.1$/;

  my $programName = $1;

  print INDEX "<p><a href=\"man1/${manFile}.html\">$programName</a></p>\n";

}

my $time = `date`;

my $footer = qq{

This document was created by
<A HREF="http://manpages.ubuntu.com/manpages/lucid/en/man1/man2html.1.html">man2html</A>,
using the manual pages.<BR>
Time: $time
</BODY>
</HTML>


};

print INDEX $footer;

close INDEX; 
