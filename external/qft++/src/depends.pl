#! /usr/bin/perl
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#
# A simple script that uses the compiler to automatically generate the 
# dependencies of qft++ source code.
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
open(outFile,">$ARGV[0]");
while(<STDIN>){
    chop; # remove end-of-line charcters
    s/\\/ /g; # remove \'s 
    if(/(\w+)\.o:(.*)/){print outFile "objects/$1.o depends/$1.d: $2"}
    else{print outFile;}
}
close(outFile);
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
