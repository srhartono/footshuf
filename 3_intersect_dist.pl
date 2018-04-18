#!/usr/bin/perl

use strict; use warnings; use Getopt::Std;
use vars qw($opt_v $opt_q $opt_f);
getopts("vqf");

my ($N,$GN,$YW,$CY,$RD,$PR,$LGN,$LCY,$LRD) = ("\e[0m","\e[0;32m","\e[1;33m","\e[0;36m","\e[0;31m","\e[0;35m","\e[1;32m","\e[1;36m","\e[1;31m","\e[1;35m");
my ($shuffledFile, $footPeakFile, $outDir, $gene, $type, $strand, $window, $threshold, $convType, $curroutdir, $label) = @ARGV;
check_sanity($shuffledFile, $footPeakFile);
$type = "NA" if not defined $type or $type eq "NOTYPE";
$gene = "UNK" if not defined $gene;
$strand = "UNK" if not defined $strand;
$window = "UNK" if not defined $window;
$threshold = "UNK" if not defined $threshold;
$convType = "UNK" if not defined $convType;
# Check sanity and get input file name and folder, defined intersect filename
my ($shuffledFolder, $shuffledFilename) = getFilename($shuffledFile, "folderfull");
my ($footPeakFolder, $footPeakFilename) = getFilename($footPeakFile, "folderfull");
my ($footPeakFilename2) = $footPeakFile =~ /($footPeakFilename.+)\.PEAK.PEAKS$/;
$footPeakFilename = defined $footPeakFilename2 ? $footPeakFilename2 : $footPeakFilename;
#my $statFile = "$footPeakFolder/$shuffledFilename\_$footPeakFilename.STAT";
#mkdir "$outDir" if not -d $outDir;
#mkdir "$outDir/TEMP" if not -d "$outDir/TEMP/";
my $statFile = "$curroutdir/$label\_$shuffledFilename\_$footPeakFilename.STAT";
print "statFile = $statFile\n";
open (my $outStat, ">", $statFile) or die "Failed to write to $statFile: $!\n";
print $outStat "#GENE\tINTTYPE\tPEAK\tODDSRATIO\tORIGPERCINTMEAN\tSHUFPERCINTMEAN\tSHUFPERCINTSD\tPVAL\tGENEFULLNAME\tSTRAND\tWINDOW\tTHRESHOLD\tCONVTYPE\tSHUFFLEDFILE\tFOOTLOOPPEAKFILE\tTOTALORIGPEAKS\tNUMBEROFSHUFFLE\tTOTALFOOTLOOPPEAKS\n";

# Parse info from input: total number of peaks from each file and total shuffled

# (Extra) Get the start site of footPeakFile

system("bedtools_bed_change.pl -q -a -x 0 -y 1 -i $footPeakFile -o $curroutdir/$label\_$footPeakFilename.BEGTEMP") == 0 or print "Failed to run bedtools_bed_change: $!\n" and die;
system("bedtools_bed_change.pl -q -c -x 0 -y 1 -i $footPeakFile -o $curroutdir/$label\_$footPeakFilename.MIDTEMP") == 0 or print "Failed to run bedtools_bed_change: $!\n" and die;
system("bedtools_bed_change.pl -q -b -x 0 -y 1 -i $footPeakFile -o $curroutdir/$label\_$footPeakFilename.ENDTEMP") == 0 or print "Failed to run bedtools_bed_change: $!\n" and die;
for (my $h = 0; $h < 4; $h++) 
{
my $footPeakFileTEMP = $h == 0 ? "$footPeakFile" : $h == 1 ? "$curroutdir/$label\_$footPeakFilename.BEGTEMP" : $h == 2 ? "$curroutdir/$label\_$footPeakFilename.MIDTEMP" : "$curroutdir/$label\_$footPeakFilename.ENDTEMP";
my $intersectOutFile = $h == 0 ? "$curroutdir/$label\_$shuffledFilename\_$footPeakFilename.INT" : $h == 1 ? "$curroutdir/$label\_$shuffledFilename\_$footPeakFilename.BEGINT" : $h == 2 ? "$curroutdir/$label\_$shuffledFilename\_$footPeakFilename.MIDINT" : "$curroutdir/$label\_$shuffledFilename\_$footPeakFilename.ENDINT";
my $intType = $h == 0 ? "WHOLE" : $h == 1 ? "BEG" : $h == 2 ? "MID" : "END";
# Run bedtools
system("cut -f1-4 $footPeakFileTEMP | bedtools intersect -wao -a $shuffledFile -b - > $intersectOutFile") == 0 or print "Failed to bedtools intersect: $!\n" and die;

# parse all peaks from shuffledFile.bed
my ($data, $origPeakNumber, $totalShuffles, $totalFootLoopPeaks, $maxInt) = parse_shuffledFile($shuffledFile, $footPeakFile);
my $peaks;
($data, $peaks) = parse_intersectFile($data, $intersectOutFile, $maxInt);

my $stat = class_perPeaks();
my %shufpeak;
$stat->{origmean} = asperc(mean($data->{orig}) / $totalFootLoopPeaks);
$stat->{origsd}   = asperc(sdev($data->{orig}) / $totalFootLoopPeaks);
#print "orig mean=$stat->{origmean}, sdev=$stat->{origsd}\n";
for (my $i = 0; $i < $totalShuffles; $i++) {
	foreach my $peak (sort keys %{$data->{shuf}}) {
		my $shufVal = $data->{shuf}{$peak}[$i];
#		$data->{shuf}{$peak1}[$count1-1] ++ if $count1 ne 0;
		push(@{$data->{shufVal}[$i]}, $shufVal);
		push(@{$data->{perPeak}{$peak}{shuf}}, $shufVal);
		$data->{perPeak}{$peak}{less} ++ if $shufVal <= $data->{perPeak}{$peak}{orig};
		$data->{perPeak}{$peak}{more} ++ if $shufVal >= $data->{perPeak}{$peak}{orig};
		my $add = $shufVal <= $data->{perPeak}{$peak}{orig} ? "" : "$LRD";
#		print "$shufVal\n" if $peak eq "POS5";
		#print "\tpeak=$peak, orig=$data->{perPeak}{$peak}{orig}, shufVal=$add$shufVal$N\n" if $peak eq "POS5"; #$i % 100 == 0;
	}
	$stat->{shufmean}[$i] = asperc(mean($data->{shufVal}[$i]) / $totalFootLoopPeaks);
	$stat->{shufsd}[$i]   = asperc(sdev($data->{shufVal}[$i]) / $totalFootLoopPeaks);
	$stat->{ntot} ++;
	$stat->{less} ++ if $stat->{shufmean}[$i] <= $stat->{origmean};
	$stat->{more} ++ if $stat->{shufmean}[$i] >= $stat->{origmean};
	#print "\t$i\tless=$stat->{less} more=$stat->{more} orig=$stat->{origmean} shuf=$stat->{shufmean}[$i], sd=$stat->{shufsd}[$i]\n" if $i % 100 == 0;
}
my $genefullname = $type eq "NA" ? $gene : "$gene\-$type";
foreach my $peak (sort keys %{$data->{perPeak}}) {
	$data->{perPeak}{$peak}{shufmean} = myformat(mean($data->{perPeak}{$peak}{shuf}));
	$data->{perPeak}{$peak}{shufsd}   = myformat(sdev($data->{perPeak}{$peak}{shuf}));
	$data->{perPeak}{$peak}{origmeanperc} = myformat($data->{perPeak}{$peak}{orig} / $totalFootLoopPeaks);
	$data->{perPeak}{$peak}{shufmeanperc} = myformat(mean($data->{perPeak}{$peak}{shuf})/$totalFootLoopPeaks);
	$data->{perPeak}{$peak}{shufsdperc}   = myformat(sdev($data->{perPeak}{$peak}{shuf})/$totalFootLoopPeaks);
	$data->{perPeak}{$peak}{p}        = $data->{perPeak}{$peak}{more} < $data->{perPeak}{$peak}{less} ? (1+$data->{perPeak}{$peak}{more}) / (1+$stat->{ntot}) : (1+$data->{perPeak}{$peak}{less}) / (1+$stat->{ntot});
	$data->{perPeak}{$peak}{p}        = myformat($data->{perPeak}{$peak}{p},1e6) if $data->{perPeak}{$peak}{p} > (1/1e6);
	$data->{perPeak}{$peak}{fold}     = myformat(($data->{perPeak}{$peak}{origmeanperc} + 10)/($data->{perPeak}{$peak}{shufmeanperc} + 10),1000);
#	print	"$peak\tfold=$data->{perPeak}{$peak}{fold}\torig=$data->{perPeak}{$peak}{orig} ($data->{perPeak}{$peak}{origmeanperc} \%)\tshuf=$data->{perPeak}{$peak}{shufmean} ($data->{perPeak}{$peak}{shufmeanperc} \%)\tsd=$data->{perPeak}{$peak}{shufsd} ($data->{perPeak}{$peak}{shufsdperc} \%)\tp=$data->{perPeak}{$peak}{p}\n";
#	print	"$peak\tfold=$data->{perPeak}{$peak}{fold}\torig=$data->{perPeak}{$peak}{origmeanperc} \%\tshuf=$data->{perPeak}{$peak}{shufmeanperc} \%\tsd=$data->{perPeak}{$peak}{shufsdperc} \%\tp=$data->{perPeak}{$peak}{p}\n";
#print $outStat "#GENE\tPEAK\tFOLD\tORIG\tSHUF\tSD\tP\tSTRAND\tWINDOW\tTHRESHOLD\tCONVTYPE\tFILE1\tFILE2\n";
	print	$outStat "$gene\t$intType\t$peak\t$data->{perPeak}{$peak}{fold}\t$data->{perPeak}{$peak}{origmeanperc}\t$data->{perPeak}{$peak}{shufmeanperc}\t$data->{perPeak}{$peak}{shufsdperc}\t$data->{perPeak}{$peak}{p}\t$genefullname\t$strand\t$window\t$threshold\t$convType\t$shuffledFile\t$footPeakFile\t$origPeakNumber\t$totalShuffles\t$totalFootLoopPeaks\n";
}
@{$stat->{shufmean}} = sort {$a <=> $b} @{$stat->{shufmean}};
my $p = $stat->{more} < $stat->{less} ? (1+$stat->{more}) / (1+$stat->{ntot}): (1+$stat->{less}) / (1+$stat->{ntot});
$p = int($p*1e6+0.5)/1e6;
my $shufmean = myformat(mean($stat->{shufmean})*100+0.5)/100;
my $shufmed = myformat($stat->{shufmean}[int(@{$stat->{shufmean}}/2)]*100+0.5)/100;
my $shufsd = myformat(mean($stat->{shufsd})*100+0.5)/100;
my $fold = myformat($stat->{origmean} / $shufmean);

print	$outStat "$gene\t$intType\tAGGREGATE\t$fold\t$stat->{origmean}\t$shufmean\t$shufsd\t$p\t$genefullname\t$strand\t$window\t$threshold\t$convType\t$shuffledFile\t$footPeakFile\t$origPeakNumber\t$totalShuffles\t$totalFootLoopPeaks\n";
}
print "Output = $statFile\n";
sub check_sanity {
	die "\nUsage: $YW$0 <${LCY}peak1.bed$N> <${LGN}peak2.bed$N>
- ${LCY}peak1.bed$N: intersectFile from 1_random_peak.pl
- ${LGN}peak2.bed$N: intersectFile PEAK.PEAKS from footPeak.pl
\n\n" unless defined $shuffledFile and -e $shuffledFile and -s $footPeakFile > 0 and defined $footPeakFile and -e $footPeakFile and -s $footPeakFile > 0;

}

sub parse_shuffledFile {
	my ($shuffledFile, $footPeakFile) = @_;
	my ($origPeakNumber, $totalShuffles) = (0,0);
	my $data;
	open (my $in1, "cat $shuffledFile|") or print "Failed to cat $shuffledFile: $!\n" and die;
	while (my $line = <$in1>) {
		chomp($line);
		my ($chr, $beg, $end, $name) = split("\t", $line);
		my ($gene, $peak1, $count1) = split(",", $name);
		if ($count1 eq 0) {
			$data->{orig}{$peak1} = 0;
			$origPeakNumber ++;
		}
		else {
			$data->{shuf}{$peak1}[$count1-1] = 0;
			$totalShuffles ++;
		}
	}
	close $in1;

	$totalShuffles /= $origPeakNumber;
	my ($totalFootLoopPeaks) = `wc -l $footPeakFile` =~ /^(\d+)/;
	my $maxIntOrig = $origPeakNumber * $totalShuffles * $totalFootLoopPeaks;
	my $maxIntShuf = $origPeakNumber * $totalFootLoopPeaks;
	#print "Total Expected Peaks: $origPeakNumber\nTotal FootLoop Peaks: $totalFootLoopPeaks\nNumber of Shuffles: $totalShuffles\n";
	#print "MaxIntOrig: $maxIntOrig\nMaxIntShuf: $maxIntShuf\n" if not defined $opt_q;
	
	return($data, $origPeakNumber, $totalShuffles, $totalFootLoopPeaks, $maxIntOrig+$maxIntShuf);
}

sub parse_intersectFile {
	my ($data, $intersectOutFile, $maxInt) = @_;
	my $linecount = 0;
	open (my $in2, "<", $intersectOutFile) or die;
	while (my $line = <$in2>) {
		chomp($line);
		$linecount ++;
		print "Done $linecount\n" if $linecount % 1000000 == 0;	
		my ($chr, $beg, $end, $name, $chr2, $beg2, $end2, $name2, @others) = split("\t", $line);
		my $int = @others == 0 ? die "Corrupted file at line=$line\n" : $others[@others-1];
		next if ($beg2 eq "-1" and $int == 0);
		my ($gene, $peak1, $count1) = $name =~ /^(.+),(.+),(\d+)$/;
		my $peak2 = $name2;
		$data->{orig}{$peak1} += $int if $count1 eq 0;
		$data->{shuf}{$peak1}[$count1-1] += $int if $count1 ne 0;
		#print "$linecount, chr=$chr, beg=$beg, end=$end, name=ame, gene=$gene, peak1=$peak1, count1=$count1, line=$line\n";
		#die if $linecount == 10;
	}
	close $in2;

	my $perc = asperc($linecount / $maxInt);
#	print "intersectFile total line: $linecount/$maxInt ($perc \%)\n" if not defined $opt_q;

	foreach my $peak (sort keys %{$data->{orig}}) {
		$data->{perPeak}{$peak} = class_perPeaks();
		#print "$peak\t$data->{orig}{$peak}\n";
		$data->{perPeak}{$peak}{orig} = myformat($data->{orig}{$peak});
	}
	return($data);
}

sub asperc {
	my ($val, $digit) = @_;
	$digit = 100 if not defined $digit;
	return int($val * $digit * 100+0.5)/$digit;
}
sub myformat {
	my ($val, $digit) = @_;
	$digit = 100 if not defined $digit;
	return int($val * $digit+0.5)/$digit;
}
sub class_perPeaks {
	my %data = (
		orig => 0,
		sdev => 0,
		p    => 0,
		more => 0,
		less => 0
	);
	return \%data;
}

sub mean {
	my ($values) = @_;
	my $mean = 0;
	return $mean if not defined $values;
	if ($values =~ /ARRAY/) {
		return(0) if @{$values} == 0;
		for (my $i =0 ; $i < @{$values}; $i++) {
			$mean += $values->[$i] / @{$values};
		}
	}
	elsif ($values =~ /HASH/) {
		my $total = (keys %{$values});
		return 0 if $total == 0;
		foreach my $key (keys %{$values}) {
#			print "key=$key, value=$values->{$key}, total=$total, meancurr=$mean\n";
			$mean += ($values->{$key} / $total);
		}
	}
	else { die "Cannot determine type of values=$values=\n";}
	return($mean);
}

sub sdev {
	my ($values) = @_;
	my $sdev = 0; my @sd;
	return $sdev if not defined $values;
	my $mean = mean($values);
	if ($values =~ /ARRAY/) {
		return(0) if @{$values} == 0;
		for (my $i =0 ; $i < @{$values}; $i++) {
	      $sdev += @{$values} == 1 ? (($values->[$i] - $mean)**2) : (($values->[$i] - $mean)**2) / (@{$values}-1);
		}
	}
	elsif ($values =~ /HASH/) {
		my $total = (keys %{$values});
		return 0 if $total == 0;
		foreach my $key (keys %{$values}) {
			push(@sd, $values->{$key});
	      $sdev += $total == 1 ? (($values->{$key} - $mean)**2) : (($values->{$key} - $mean)**2) / ($total-1);
#			print "key=$key, value=$values->{$key}, total=$total, meancurr=$mean, sdevcurr=$sdev\n";
		}
	}
	else { die "Cannot determine type of values=$values=\n";}
	$sdev = sqrt($sdev);
	return($sdev);
}

sub getFilename {
    my ($fh, $type) = @_;
   die "INPUT <@_> /usr/local/bin/Perl/mitochy.pm: getFilename <fh> <type (folder, full, folderfull, all)\n" if not defined $fh or not -e $fh;
   my $fh0 = $fh;
   $fh =~ s/\/+/\//g;
   my ($first, $last) = (0,0);
   if ($fh =~ /^\//) {
      $first = 1;
      $fh =~ s/^\///;
   }
   if ($fh =~ /\/$/) {
      $last = 1;
      $fh =~ s/\/$//;
   }
   # Split folder and fullname
   my (@splitname) = split("/", $fh);
   my $fullname = pop(@splitname);
   my @tempfolder = @splitname;
   my $folder;
   if (@tempfolder == 0) {
      $folder = "./";
   }
   else {
      $folder = join("\/", @tempfolder) . "/";
   }
   #$folder = "./" if $folder eq "/";
#  print "\nFolder = $folder\nFile = $fh, fullname=$fullname\n\n";
   # Split fullname and shortname (dot separated)
        @splitname = split(/\./, $fullname);
        my $shortname = $splitname[0];
   if ($first == 1) {$folder = "/$folder";}
   if ($last == 1) {$fullname = "$fullname/";}
   #print "\nFh=$fh0\nFolder=$folder\nShortname=$shortname\nFullname=$fullname\n\n";
        return($shortname)          if not defined($type);
        return($folder, $fullname)     if defined($type) and $type =~ /folderfull/;
        return($folder, $shortname)    if defined($type) and $type =~ /folder/;
        return($fullname)        if defined($type) and $type =~ /full/;
        return($folder, $fullname, $shortname)  if defined($type) and $type =~ /all/;
}


__END__
for (my $i = 0; $i < $totcount; $i++) {
	$data->[$i] = 0 if not defined $data->[$i];
	tot ++;
	$data->{more} ++ if $data->[$i] >= $orig;
	$less ++ if $data->[$i] <= $orig;
	push(@shuf, $data->[$i] / $total);
}
$data->{more} /= tot;
$less /= tot;
my $shufmean = int(mean(@shuf)*10000+0.5)/100;
my $data->{origmean} = int($orig / $total * 10000+0.5)/100;

__END__
my ($less, $data->{more}, tot) = (0,0,0);
for (my $i = 0; $i < $totcount; $i++) {
	$data->[$i] = 0 if not defined $data->[$i];
	tot ++;
	$more ++ if $data->[$i] >= $orig;
	$less ++ if $data->[$i] <= $orig;
	push(@shuf, $data->[$i] / $total);
}
$more /= tot;
$less /= tot;
my $p = $more < $less ? $more : $less;
my $shufmean = int(mean(@shuf)*10000+0.5)/100;
my $data->{origmean} = int($orig / $total * 10000+0.5)/100;
print "orig=$orig, total=$total, mean=$data->{origmean}, shuf=$shufmean, p=$p\n";



__END__
my ($shuffledFolder, $shuffledFilename) = mitochy::getFilename($shuffledFile, "folder");


my %data;
open (my $in1, "<", $shuffledFile) or die "Cannot read from $shuffledFile: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /^#/;
	my @arr = split("\t", $line);
}
close $in1;


open (my $out1, ">", "$shuffledFilename.out") or die "Cannot write to $shuffledFilename.out: $!\n";
close $out1;

