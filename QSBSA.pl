#!/usr/bin/perl -s 
use Cwd 'abs_path';
use Cwd;

if (not defined($f)) { 
	print "\n###################################################  Qiao-shan's Bulk Segregation Analysis Program  ###################################################\n\n";
	print "Usage: perl /PATH/QSBSA.pl -f=<VCF_FILE> [-m=off] [-d=300] [-n=1] [-g=off] [-w=20000] [-o=SNPCount] [-t=off] [-v=/home/CAM/qlin/BSA/LF10_vcf.txt] \n\n";
	print "Options:\n\n";
	print "-m\tUse mapping to SL9 method to find candidates. Set -m=on to turn on this mode. Default:off\n";
	print "Filter by divergence\n";
	print "-d\tNumber of base pairs (distance between SNPs). Default:300\n";
	print "-n\tMaximum number of SNPs within the given distance. Default:1\n";
	print "-g\tGreedy mode (if the number of SNPs within any window of the given distance more than n, don't count any of them).\n";
	print "\tSet -g=on to turn on this mode. Default:off\n";
	print "SNP counting in windows\n";
	print "-w\tWidth of window in which valid SNPs are counted. Default:20000\n";
	print "Output\n";
	print "-o\tOutput file prefix. Default:SNPCount\n";
	print "-t\tText output only (no picture output).\n";
	print "\tSet -t=on to turn on. Default:off\n";
	print "Mutants compare filter\n";
	print "-v\tCompare vcf file with other mutants' vcf files to filter common SNPs.\n";
	print "\tSet (e.g. -v=VCF_files.txt) to provide a list of vcf files with their paths (each file in a line). Set -v=off to turn off this function. Default:/home/CAM/qlin/BSA/LF10_vcf.txt\n\n";
	exit;
}

my %length; #record contig length
my %positions; #record SNP positions on a contig
my %legalSNPs; #record legal SNP positions on a contig
my $max = 0; #record the maximum #SNP in a window
$d = 300 if not defined($d);
$n = 1 if not defined($n);
$w = 20000 if not defined($w);
$o = "SNPCount" if not defined($o);
$g = "off" if not defined($g);
$t = "off" if not defined($t);
$m = "off" if not defined($m);
$v = "/home/CAM/qlin/BSA/LF10_vcf.txt" if not defined($v);

open IN, "$f" or die "Can't open $f: $!";
while ($_ = <IN> and $_ =~ /^##/) { 
	$length{$1} = $2 if $_ =~ /^##contig=<ID=(.+),length=(.+)>$/;
}
while (<IN>) {
	@_ = split(/\t/,$_);
	$positions{$_[0]} .= "$_[1]\t";
}
close IN;

if ($v ne "off") {
	open IN, "$v" or die "Can't open $v: $!";
	my $dir = getcwd;
	foreach my $file (<IN>) {
		next if $file eq "" or $file =~ /^(\s)*$/ or $file =~ /^$dir\//;
		open VCF, "$file" or die "Can't open $file: $!";
		while (<VCF>) {
			next if $_ =~ /^#/;
			@_ = split(/\t/,$_);
			$positions{$_[0]} =~ s/\t$_[1]\t/\t/ if $positions{$_[0]} =~ /\t$_[1]\t/;
			$positions{$_[0]} =~ s/^$_[1]\t// if $positions{$_[0]} =~ /^$_[1]\t/;
		}
		close VCF;
	}

	open IN, "$f" or die "Can't open $f: $!";
	open OUT, ">$o.vcf";
	while (<IN>) {
		next if $_ =~ /^#/;
		@_ = split(/\t/,$_);
		print OUT "$_" if $positions{$_[0]} =~ /\t$_[1]\t/;
		print OUT "$_" if $positions{$_[0]} =~ /^$_[1]\t/;
	}
	close IN;
	close OUT;
}

if ($m eq "on") {
	foreach my $s (keys %positions){
		$positions{$s} =~ s/\t$//;
		my @pos = split(/\t/,$positions{$s});
		my @inner;
		for (my $i = 0; $i < scalar(@pos)-$n; $i++) {   $pos[$i+$n]-$pos[$i] > $d ? push(@inner, 1) : push(@inner, 0)   }
		my $pos;
		if ($g eq "on") {
	   		$pos = '1' x scalar(@pos);
			for (my $i = 0; $i < scalar(@inner); $i++) {    substr($pos, $i, $n+1) = '0' x ($n+1) if $inner[$i]==0    }
	   	}else{
			$pos = '0' x scalar(@pos);
			for (my $i = 0; $i < scalar(@inner); $i++) {	substr($pos, $i, $n+1) = '1' x ($n+1) if $inner[$i]==1    }
	   	}
		my @validSites = split(//,$pos);
	   	for (my $i = 0; $i < scalar(@pos); $i++) {    $legalSNPs{$s} .= "$pos[$i]\t" if $validSites[$i]==1    }
	}

	open OUT, ">$o.txt";
	print OUT "contig\tposition\tcount\n";
	foreach my $s (sort keys %legalSNPs) {
		$legalSNPs{$s} =~ s/\t$//;
		my @pos = split(/\t/,$legalSNPs{$s});
		for (my $i = 0; $i < $length{$s}; $i+=$w) {
			my $count = 0;
			foreach (@pos) {	$count++ if ($_ <= $i+$w && $_ >= $i+1)    }
			$max = $count if $count > $max;
			$tmps = $s;
			$tmps =~ s/^SL9_pseudoscaffold_// if $tmps =~ /^SL9_pseudoscaffold_/;
			print OUT $tmps,"\t",$i+1,"\t$count\n";
		}
	}
	close OUT;

	my $abs = abs_path($0);
	$abs =~ s/\/QSBSA.pl$//;
	$abs =~ s/ /\\ /g;
	$abs =~ s/\'/\\'/g;
	system("Rscript $abs/QSBSA.r $o $max") if $t eq "off";
}
