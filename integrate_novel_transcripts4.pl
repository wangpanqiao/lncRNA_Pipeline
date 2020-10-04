#!/usr/bin/perl -w
use strict;

open FH,"novel.longRNA.CPAT.out" or die;   ####CPAT默认把trans_id改成全大写
<FH>;     #####过滤第一行,相当于readline FH; # skip the first line
####相当于直接定义标量变量运行一遍： my $skipfirst = <FH>; 

my %info;
while(<FH>){
	chomp;
	my @field=split "\t";
	if($field[5]>=0.364){ ###第6列，第一列没有header列名
		$info{$field[0]}{CPAT}="TUCP";  ###trans_ID大写转小写,如果本来就是大写,相当于没变化
	}else{
		$info{$field[0]}{CPAT}="lncRNA";
	}
}

open FH,"novel.longRNA.PLEK.out" or die;

while(<FH>){
	chomp;
	my @field=split "\t";
	# $field[2]=~/>(.+?) /;  ###提取id
	# my $id=$1;
	my @name=split(/\s/,$field[2]);
	my $id=$name[0];
	$id=~s/>//g;
	if($field[0] eq "Coding"){
		$info{$id}{PLEK}="TUCP";
	}else{
		$info{$id}{PLEK}="lncRNA";
	}
}

open FH,"novel.longRNA.CNCI.out" or die;
readline FH; # skip the first line
while(<FH>){
	chomp;
	my @field=split "\t";
	# $field[0]=~/^(.+?) /;   ###提取ID
	# my $id=$1;
	my @name=split(/\s/,$field[0]);
	my $id=$name[0];
	if($field[1] eq "coding"){
		$info{$id}{CNCI}="TUCP";
	}else{
		$info{$id}{CNCI}="lncRNA";
	}
}

open FH,"novel.longRNA.CPC2.out" or die;
readline FH; # skip the first line
while(<FH>){
	chomp;
	my @field=split "\t";
	my $id=$field[0];
	if($field[7] eq "coding"){
		$info{$id}{CPC2}="TUCP";
	}else{
		$info{$id}{CPC2}="lncRNA";
	}
}

open FH,"novel.longRNA.exoncount.txt" or die;

while(<FH>){
	chomp;
	my @field=split "\t";
	$info{$field[0]}{EXONCOUNT}=$field[1];
}

print "id\tCPAT\tPLEK\tCNCI\tCPC2\tEXONCOUNT\tintegrate\n";
foreach my $id (sort keys %info){
	print $id."\t".$info{$id}{CPAT}."\t".$info{$id}{PLEK}."\t".$info{$id}{CNCI}."\t".$info{$id}{CPC2}."\t".$info{$id}{EXONCOUNT}."\t";
	if($info{$id}{CPAT} eq "lncRNA" && $info{$id}{PLEK} eq "lncRNA" && $info{$id}{CNCI} eq "lncRNA" && $info{$id}{CPC2} eq "lncRNA"){
		print "lncRNA\n";
	}else{
		print "TUCP\n";
	}
	
}
