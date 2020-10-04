#!/usr/bin/perl -w
use strict;

my %know_lnc;
open FH,"known.lncRNA.bed" or die;
while(<FH>){
	chomp;
	my @field=split "\t";
	$know_lnc{$field[0].'\t'.$field[1].'\t'.$field[2].'\t'.$field[5].'\t'.$field[7]} = $field[3];
}


my %genecode;
open FH,"gencode.v25.annotation.chrX.gtf_mod.gtf" or die;
while(<FH>){
	chomp;
	my @field=split "\t";
	# $_=~/gene_name "(.+?)"/; ##修改
	$_=~/gene_id "(.+?)"/;
	my $gene_name=$1;
	my $loc = $field[0].'\t'.($field[3]-1).'\t'.$field[4].'\t'.$field[6].'\t'.$field[2];
	foreach my $location (keys %know_lnc){
		if($location eq $loc){
			$genecode{$know_lnc{$loc}} = $gene_name;
		}
	}
}
open FH,"lncipedia_4_0.chrX.gtf_mod.gtf" or die;
while(<FH>){
	chomp;
	my @field=split "\t";
	$_=~/gene_id "(.+?)"/;
	my $gene_name=$1;
	my $loc = $field[0].'\t'.($field[3]-1).'\t'.$field[4].'\t'.$field[6].'\t'.$field[2];
	foreach my $location (keys %know_lnc){
		if($location eq $loc){
			$genecode{$know_lnc{$loc}} = $gene_name;
		}
	}
}

my %exon;
my %gene;

open FH,"known.lncRNA.bed" or die;
while(<FH>){
	chomp;
	my @field=split "\t";
	my $chr=$field[0];
	my $genename="known:".$field[3];
	my $start=$field[1];
	my $end=$field[2];
	my $strand=$field[5];
	$_=~/transcript_id "(.+?)"/;
	my $transid=$1;
	my $loc=$chr."\t".($start+1)."\t".$end."\t".$strand;
	$exon{$genename}{$transid}{$loc}=1;
	$gene{$genename}{CHR}=$chr;
	$gene{$genename}{STRAND}=$strand;
	if(exists $gene{$genename}{START}){
		$gene{$genename}{START}=$gene{$genename}{START}<$start?$gene{$genename}{START}:$start;
	}else{
		$gene{$genename}{START}=$start;
	}
	if(exists $gene{$genename}{END}){
		$gene{$genename}{END}=$gene{$genename}{END}>$end?$gene{$genename}{END}:$end;
	}else{
		$gene{$genename}{END}=$end;
	}
	
}

open FH,"novel.lncRNA.stringent.filter.bed" or die;   ###利用BEDOPS的gtf2bed转换
while(<FH>){
	chomp;
	my @field=split "\t";
	my $chr=$field[0];
	my $genename="novel:".$field[3];    ###为gene PASA_cluster添加novel,添加内容，对应novel
	# my $genename="lncRNA:".$field[3];    ###为gene PASA_cluster添加lncRNA,添加内容，对应lncRNA
	my $start=$field[1];
	my $end=$field[2];
	my $strand=$field[5];
	$_=~/transcript_id "(.+?)"/;   ####bed文件需要过滤exon和transcript的
	my $transid=$1;
	my $loc=$chr."\t".($start+1)."\t".$end."\t".$strand;  ###
	$exon{$genename}{$transid}{$loc}=1;     ###exon  novel:PASA align_id:179956|asmbl_35  C01 194329  194719 + =1
	$gene{$genename}{CHR}=$chr;
	$gene{$genename}{STRAND}=$strand;
	if(exists $gene{$genename}{START}){
		$gene{$genename}{START}=$gene{$genename}{START}<$start?$gene{$genename}{START}:$start;
	}else{
		$gene{$genename}{START}=$start;
	}
	if(exists $gene{$genename}{END}){
		$gene{$genename}{END}=$gene{$genename}{END}>$end?$gene{$genename}{END}:$end;
	}else{
		$gene{$genename}{END}=$end;
	}
	
}
open OUT,">lncRNA.for_anno.bed" or die;
foreach my $k (keys %gene){
	print OUT $gene{$k}{CHR}."\t".$gene{$k}{START}."\t".$gene{$k}{END}."\t".$k."\t.\t".$gene{$k}{STRAND}."\n";
}

`sort-bed lncRNA.for_anno.bed > lncRNA.for_anno.srt.bed`;

###BEDOPS的closest输入必须为sorted。离得最近gene，重叠dist为0
`closest-features --dist lncRNA.for_anno.srt.bed  gencode.protein_coding.gene.bed > lncRNA.for_anno.srt.neighbour.txt`;

open FH,"lncRNA.for_anno.srt.neighbour.txt" or die;   ###bed文件没有exon信息
my %map;    ###
my %genename2gene;
my $naidx=0;
# my $num=0;
# my $sum1=0;
# my $sum2=0;
# my $sum3=0;
# my $sum4=0;
my %lnciscodegene;
my %genename2gene1;
while(<FH>){
	chomp;
	my @field=split /\|/;
	my ($chr,$start,$end,$geneid,$tmp2,$strand)=split "\t",$field[0];   ###novel:PASA_cluster_28,geneid为PASA的geneid
	if($strand ne "+" and $strand ne "-"){
		#print $field[0]."\n";
		next;
	}
	my $up_gene="";
	my $up_dist=999999;
	my $up_strand;
	my $down_gene="";
	my $down_dist=999999;
	my $down_strand;
	# my $closet_genename="";
	if($field[1] ne "NA"){
		# $field[1]=~/gene_name "(.+?)"/;
		$field[1]=~/ID=(.+?)$/;         ####因为gencode.protein_coding.gene.bed来自gff  [9]列是ID=,而不是gtf的gene_name ""
		$up_gene=$1;
		$up_dist=abs($field[2]);
		my @tmp=split "\t",$field[1];
		$up_strand=$tmp[5];
	}
	if($field[3] ne "NA"){
		# $field[3]=~/gene_name "(.+?)"/;
		$field[3]=~/ID=(.+?)$/;   ###同上
		$down_gene=$1;
		$down_dist=abs($field[4]);
		my @tmp=split "\t",$field[3];
		$down_strand=$tmp[5];
	}
	if($field[1] eq "NA" and $field[3] eq "NA"){  ####BEDOPS_closest-features两侧没有coding基因
		$naidx++;
		my $genename="NA-$naidx";
		$map{$genename}{$geneid}="NA";
		#print $genename."\n";
	}else{	
	if($up_dist < $down_dist){
		my $genename;
		if($up_dist==0){
			if($strand ne $up_strand){
				$genename=$up_gene."-AS";  ###up_strand  gene  chh11G006650   邻近基因name
				$genename2gene{$genename} = $up_gene;     ####邻近基因name-AS=邻近基因
				if(exists $map{$genename}{$geneid}){
					$map{$genename}{$geneid}=$map{$genename}{$geneid} > $up_dist?$up_dist:$map{$genename}{$geneid};
				}else{
					$map{$genename}{$geneid}=$up_dist;   ###map{邻近基因name}{PASA_gene_id}=距离
				}
				
			}else{
				$genename=$up_gene."-INside";  ###up_strand  gene  chh11G006650   自己添加
				$genename2gene1{$genename} = $up_gene;
				if(exists $map{$genename}{$geneid}){
					$lnciscodegene{$genename}{$geneid}=$lnciscodegene{$genename}{$geneid} > $up_dist?$up_dist:$lnciscodegene{$genename}{$geneid};
				}else{
					#print "LALALALALA"."\n";   ###有很多PASA属于chc01G00320基因内转录
					$lnciscodegene{$genename}{$geneid}=$up_dist;
				}
			}
		}else{
			$genename="LINC-".$up_gene;
			$genename2gene{$genename} = $up_gene;
			if(exists $map{$genename}{$geneid}){
					$map{$genename}{$geneid}=$map{$genename}{$geneid} > $up_dist?$up_dist:$map{$genename}{$geneid};
				}else{
					$map{$genename}{$geneid}=$up_dist;
				}
		}
	}else{
		my $genename;
		if($down_dist==0){
			if($strand ne $down_strand){
				$genename=$down_gene."-AS";
				$genename2gene{$genename} = $down_gene;
				if(exists $map{$genename}{$geneid}){
					$map{$genename}{$geneid}=$map{$genename}{$geneid} > $down_dist?$down_dist:$map{$genename}{$geneid};
				}else{
					$map{$genename}{$geneid}=$down_dist;
				}
			}else{   ###自己添加
				$genename=$down_gene."-INside";
				$genename2gene1{$genename} = $down_gene;
				if(exists $map{$genename}{$geneid}){
					$lnciscodegene{$genename}{$geneid}=$lnciscodegene{$genename}{$geneid} > $down_dist?$down_dist:$lnciscodegene{$genename}{$geneid};
				}else{
					#print "HAHAHAHAHAHAA"."\n";   ######有很多PASA属于chc01G00320基因内转录
					$lnciscodegene{$genename}{$geneid}=$down_dist;
				}
			}
		}else{
			$genename="LINC-".$down_gene;
			$genename2gene{$genename} = $down_gene;
			if(exists $map{$genename}{$geneid}){
					$map{$genename}{$geneid}=$map{$genename}{$geneid} > $down_dist?$down_dist:$map{$genename}{$geneid};
				}else{
					$map{$genename}{$geneid}=$down_dist;
				}
		}
	}
	}
}
# my length=0;
# while (($key, $value) = each %map) {
# length++;
# }

###自己添加，统计哈希
my $length = keys %map;
print STDERR "map Length : $length\n";  ###498

# my $size=scalar keys%map;
# print STDERR "top-level map $size\n";
# foreach  my $key (keys %map)
# {
    # $size=scalar keys%{$map{$key}};
    # print STDERR "second-level map $size\n";

# }
# foreach my $key (keys %map)
# {
    # foreach my $subkey (keys %{$map{$key}})
    # {
        # print STDERR "$key\t$subkey\n";
    # }
# }

my $size=scalar keys%genename2gene;
print STDERR "top-level genename2gene $size\n";   ###492

my %MSTRG2genename;
open OUT1,">lncRNA.final.v2.gtf" or die;
open OUT2,">lncRNA.final.v2.map" or die;
foreach my $genename (keys %map){    ##第一层key是临近基因name-AS chc可以重复
	my $tmp1=$map{$genename};
	my %tmp1=%$tmp1;    ##第二层哈希  PASA gene id = 距离
	my $gindex=0;
	my $out_gene = $genename2gene{$genename};   ###genename2gene哈希中< map哈希   ，等于邻近基因name，去掉AS
	foreach my $id (sort {$tmp1{$a}<=>$tmp1{$b}} keys %tmp1){
		my ($class,$cuffid)=split ":",$id;    ###分成 对应novel PASA_cluster_28   或对应lncRNA
		my $dist=$tmp1{$id};
		my $tmp2=$exon{$id};    ####需要检查%exon(哈希),把transcript加入了
		my %tmp2=%$tmp2;
		$gindex++;
		my $geneid=$genename."-".$gindex;
		$MSTRG2genename{$cuffid} = $geneid;
		my $tindex=0;
		print OUT2 $geneid."\t".$cuffid."\n";
		foreach my $tid (keys %tmp2){
			my $tmp3=$tmp2{$tid};
			my %tmp3=%$tmp3;
			$tindex++;
			# my $transid=$geneid.":".$tindex;
			my $transid=$geneid."\.".$tindex;  ###修改：为.
			foreach my $exon (keys %tmp3){
				my @loc=split "\t",$exon;
				print OUT1 $loc[0]."\t$class\texon\t".$loc[1]."\t".$loc[2]."\t.\t".$loc[3]."\t.\tgene_id \"$geneid\"; transcript_id \"$transid\"; dist \"$dist\"; "."closet_gene \"".$out_gene."\";"."\n";
			}

		}
	}
}

###自己添加，输出lncRNA in codinggene
my %MSTRG2genename1;
open OUT4,">lncRNAiscoding.v2.gtf" or die;
open OUT5,">lncRNAiscoding.v2.map" or die;
foreach my $genename (keys %lnciscodegene){
	my $tmp1=$lnciscodegene{$genename};
	my %tmp1=%$tmp1;
	my $gindex=0;  ###gene index
	my $out_gene = $genename2gene1{$genename};   ###genename2gene哈希中< map哈希
	foreach my $id (sort {$tmp1{$a}<=>$tmp1{$b}} keys %tmp1){   ###同一个基因附近 novel:PASA 的个数
		my ($class,$cuffid)=split ":",$id;
		my $dist=$tmp1{$id};
		my $tmp2=$exon{$id};  ## exon  novel:PASA  align_id:179956|asmbl_35  C01 194329  194719 + =1
		my %tmp2=%$tmp2;
		$gindex++;
		my $geneid=$genename."-".$gindex;
		$MSTRG2genename1{$cuffid} = $geneid;   ###PASA_gene=邻近基因-linc-1
		my $tindex=0;   ###trans  index
		print OUT5 $geneid."\t".$cuffid."\n";
		foreach my $tid (keys %tmp2){
			my $tmp3=$tmp2{$tid};
			my %tmp3=%$tmp3;
			$tindex++;
			my $transid=$geneid."\.".$tindex;
			foreach my $exon (keys %tmp3){
				my @loc=split "\t",$exon;
				print OUT4 $loc[0]."\t$class\texon\t".$loc[1]."\t".$loc[2]."\t.\t".$loc[3]."\t.\tgene_id \"$geneid\"; transcript_id \"$transid\"; dist \"$dist\"; "."closet_gene \"".$out_gene."\";"."\n";
			}

		}
	}
}

open OUT3,">lncRNA.mapping.file" or die;
foreach my $mstr(sort(keys %genecode)){
	if(defined($MSTRG2genename{$mstr})){
	print OUT3 $mstr."\t".$genecode{$mstr}."\t".$MSTRG2genename{$mstr}."\n";
	}else{
		#print $mstr."\n";
	}
}