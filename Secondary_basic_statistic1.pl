#!/usr/bin/perl -w
#since CPAT arbitrarily transforms gene names into upper case, we apply 'uc' function to keep the genenames' consistency.  
use strict;
open OUT,">basic_charac.txt" or die;

open FH,"all_lncRNA_for_classifier.gtf" or die;

my %class;
my %g2t;
my %t2t;
my %trans_len;
my %exon_num;
while(<FH>){
	chomp;
	my @field=split "\t";
	$_=~/gene_id "(.+?)"/;
	my $gid=$1;
	$_=~/transcript_id "(.+?)"/;
	my $tid=uc($1);   ####uc是大写，表示upper convert，lc是小写，表示lower conver
	$t2t{$tid}=$1;
	$class{$tid}=$field[1];
	$g2t{$tid}=$gid;
	my $len=$field[4]-$field[3];
	$trans_len{$tid}=(exists $trans_len{$tid})?$trans_len{$tid}+$len:$len;
	$exon_num{$tid}=(exists $exon_num{$tid})?$exon_num{$tid}+1:1;
}
open FH,"protein_coding.final.gtf" or die;

while(<FH>){
	chomp;
	my @field=split "\t";
	# if ($field[2] ne "exon"){ 
		# next;
	# }
	next unless $field[2] eq 'exon';
	$_=~/gene_id "(.+?)"/;
	# my $gid=uc($1);   ###uc是大写，表示upper convert，lc是小写，表示lower conver
	my $gid=$1;   ###
	$_=~/transcript_id "(.+?)"/;
	my $tid=uc($1);
	$t2t{$tid}=$1;
	$class{$tid}="protein_coding";
	$g2t{$tid}=$gid;
	my $len=$field[4]-$field[3];
	$trans_len{$tid}=(exists $trans_len{$tid})?$trans_len{$tid}+$len:$len;
	$exon_num{$tid}=(exists $exon_num{$tid})?$exon_num{$tid}+1:1;
}

my %lin_class;
open IN,"lncRNA_classification.txt" or die;                 #change the file name
while(<IN>){
	chomp;
	my @data = split /\t/,$_;
	my $type = $data[1];
	$type=~ s/\s/_/g;
	$lin_class{$data[0]} = (exists $lin_class{$data[0]})?$lin_class{$data[0]}.",".$type:$type;
}


open FH,"lncRNA.final.CPAT.out" or die;

<FH>;

while(<FH>){
    chomp;
    my @field=split "\t";
    my $tid=uc($field[0]);   ####CPAT结果trans_id被转成大写
    my $old_tid=$t2t{$tid};
    my $class;
    if (defined($lin_class{$old_tid})){
        $class = $lin_class{$old_tid};
    }else{
        $class = 'NA';
    }
    print OUT $g2t{$tid}."\t".$old_tid."\t".$class{$tid}."\t".$field[5]."\t".$trans_len{$tid}."\t".$exon_num{$tid}."\t".$class."\n";
}
    
open FH,"protein_coding.final.CPAT.out" or die;

<FH>;
            
while(<FH>){
    chomp;
    my @field=split "\t";
    my $tid=uc($field[0]);   ####CPAT结果trans_id被转成大写
    my $old_tid=$t2t{$tid};
    my $class;
    if (defined($lin_class{$old_tid})){
        $class = $lin_class{$old_tid};
    }else{
        $class = 'protein_coding';
    }
    print OUT $g2t{$tid}."\t".$old_tid."\t".$class{$tid}."\t".$field[5]."\t".$trans_len{$tid}."\t".$exon_num{$tid}."\t".$class."\n";
 }