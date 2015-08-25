#!/usr/bin/perl
use warnings;

#Command Line
#./ReformatForDeconvolution.pl outputpoolRIce 90

$num= $ARGV[1];
$pool = 'pool' . $num . ':';

open (IN, "<$ARGV[0]") || die "Error";
while (!eof(IN))
{
        chomp ( $_=<IN>);
	@field = split /$pool/, $_;
	if ($field[1] !~ /[0-9]/){
		print "$field[0]" . "$pool $num\n";
	}
	else{
		print "$_\n";
	}

}
close(IN);
