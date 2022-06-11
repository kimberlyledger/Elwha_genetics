#/usr/bin/perl -w
use strict;

my$keepLine=0;
my$sample;
my$sampleNum=0;
print"phaseID\tphaseSeq\n";
while(my$line=<>){
	chomp $line;
	if($line=~/^BEGIN GENOTYPES/){
		$keepLine=1;
		next;
	}
	
	if($line=~/^END GENOTYPES/){
		$keepLine=0;
		next;
	}
	
	if($keepLine==1){
		if($line=~/#/){
			$sample=$line;
			$sample=~s/^# //;
			$sampleNum=0;
		}else{
			$sampleNum++;
			my$phaseID=$sample."_".$sampleNum;
			my$phaseSeq=$line;
			$phaseSeq=~s/\s//g;
			print"$phaseID\t$phaseSeq\n";
		}
	}
}