#!/usr/bin/perl

if(($#ARGV + 1) != 2)
{
	print "Usage: contigLength.pl <INPUT_FILE> <OUTPUT_FILE>\n";
}

$filenameIn = $ARGV[0];
$filenameOut = $ARGV[1];

open $FILEIN, "< $filenameIn";
open $OUTFILE, "> $filenameOut";

while($_ = <$FILEIN>) {
	
	$tmp = substr $_ ,0 ,1;	

	if($tmp eq ">")
	{
		if($len != 0)
		{
			print $OUTFILE join "", $len, "\n"; 
		}
		
		$len = 0;
	}
	
	else
	{
		$len = $len + length($_) - 1;
	}
}

print $OUTFILE join "", $len, "\n";
