#!/usr/bin/perl -w
use strict;

my $kmer=$ARGV[0];
my $line="";
my @col=();
my $tmp="";
my $N=0;
my $i="";
my %sample=();
my %sex=();

my $female=0;
my $male=0;
my $exist=0;
my $femaleT=0;
my $maleT=0;
my %score=();
my %recored=();


$female=0;
$male=0;
open IN,"sex.txt";
while($line=<IN>){
	chomp($line);
	@col=split(/\t/,$line);
	$sex{$col[0]}=$col[1];
	#print $col[0],"\t",$sex{$col[0]},"\n";
}
close IN;

$N=0;
open IN, "list";
while($line=<IN>){
        $N++;
        ($tmp)=($line=~/(.+) :/);
        #print $tmp,"\t",$N,"\n";
        $sample{$N}=$tmp;
	if($sex{$tmp}==1){$male++}
	elsif($sex{$tmp}==2){$female++}
}
close IN;

$femaleT=$female*0.1;
$maleT=$male*0.9;

$tmp="count.".$kmer.".matrix";
open RECORD,">$tmp";
print RECORD "kmer\ttotal\tmale($male)\tfemale($female)\n";
$tmp="./matrix_".$kmer;
open IN, "kmtricks aggregate --run-dir $tmp --matrix kmer -t 10 --format text --cpr-in --sorted |" or die "$!";
while($line=<IN>){
	chomp($line);
	#print STDERR $line,"\n";
	@col=split(/\s/,$line);
	#print STDERR  $col[0],"\t",$col[1],"\n";
	$female=0;$male=0;$exist=0;
	foreach $i (1 .. $#col){
		if($sex{$sample{$i}}==2){
			if($col[$i]>0){
				$female++;
				$exist++;
			}
		}
		elsif($sex{$sample{$i}}==1){
			if($col[$i]>0){
				$male++;
				$exist++;
			}
		}
	}
	if($female<$femaleT and $male>$maleT){
		print RECORD $col[0],"\t",$exist,"\t",$male,"\t",$female,"\t";
		for $i (1 .. $#col){
			print RECORD $col[$i]," ";
		}
		print RECORD "\n";
	}
}
close IN;
close RECORD;

