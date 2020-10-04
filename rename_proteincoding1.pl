#!/usr/bin/perl -w
use strict;

open FH,shift or die;

while(<FH>){
	chomp;
	my @field=split "\t";
	my @id=split(";",$field[8]);
	my $f8=$id[1]."; ".$id[0].";\n";
	print (join(qq{\t},@field[0..7],$f8));
	# $_=~/gene_name "(.+?)";/;
	# my $genename=$1;
	# $_=~s/gene_id ".+?"/gene_id "$genename"/;
	# print $_;
}
