#!/usr/bin/perl

use strict; use warnings; use Getopt::Std;
use vars qw($opt_v $opt_f $opt_w $opt_t $opt_g $opt_s $opt_c $opt_o $opt_i $opt_a $opt_b $opt_x $opt_h $opt_F);
getopts("vfw:t:g:s:cho:i:a:b:xfF");

my ($N,$GN,$YW,$CY,$RD,$PR,$LGN,$LCY,$LRD) = ("\e[0m","\e[0;32m","\e[1;33m","\e[0;36m","\e[0;31m","\e[0;35m","\e[1;32m","\e[1;36m","\e[1;31m","\e[1;35m");
my ($script, $dir1, $dir2) = ($opt_i, $opt_a, $opt_b);
die "\nusage: $YW$0$N [options] $LCY-o$N <Output Dir> $LCY-i$N <3_intersect.pl> $LCY-a$N <dir from 1_random_peak.pl> $LCY-b$N <footPeak dir>

This script will intersect predicted peaks with footloop peaks 
if predicted peaks are Pos strand (e.g. CALM3) then footloop peaks intersect will only be Pos strand, and vice versa

-x: DRY RUN (don't run the sbatch)

Options (without these it'll just intersect all files in -b with -a)
-w (e.g. 10 for 10bp) intersect files with this window
-t (e.g. 0.65) intersect files with this threshold
-g (e.g. CALM3) intersect files with this gene name

C or G context options (default: use all)
-c to use CG or GC files
-h to use CH or GH files

\n" unless defined $opt_i and -e $opt_i and defined $opt_a and -e $opt_a and defined $opt_b and -e $opt_b;

die "-t must be number ($opt_t)\n" if defined $opt_t and $opt_t !~ /^\d+\.?\d*e?\-?\d*$/;
die "-w must be number ($opt_w)\n" if defined $opt_w and $opt_w !~ /^\d+$/;
die "-o (Output Dir) isn't defined!\n" if not defined $opt_o;
die "-o output dir $LCY$opt_o$N already exist! (use -f to overwrite at your own risk)\n\n" if -d $opt_o and not defined $opt_f;



my %data; my %strand;
$script = "perl $script";
my @dir1 = <$dir1/*.bed>;
my ($labelz) = `cat $dir1/.LABEL`; chomp($labelz);
die "There is not .bed file in -a $dir1\n" if @dir1 == 0;
foreach my $file1 (@dir1) {
	my $file1TEMP = $file1;
	($file1TEMP) = getFilename($file1TEMP, 'full');
	my ($gene, $type, $number) = $file1TEMP =~ /^$labelz\_(\w+)_(\w)_(\d+).bed$/;
		($gene, $number) = $file1TEMP =~ /^$labelz\_([A-Za-z0-9_]+)_(\d+).bed$/ if $file1TEMP =~ /^$labelz\_([A-Za-z0-9_]+)_(\d+).bed$/;
		($gene, $number) = $file1TEMP =~ /^$labelz\_([A-Za-z0-9]+)_(\d+).bed$/ if $file1TEMP =~ /^$labelz\_([A-Za-z0-9]+)_(\d+).bed$/;
		$type = "NOTYPE" if $file1TEMP =~ /^$labelz\_([A-Za-z0-9_]+)_(\d+).bed$/;
		$gene = uc($gene);

	#print "file1=$file1TEMP gene=$gene, type=$type, number=$number\n";
	print "Label=$labelz, File = $file1TEMP, gene=$gene, number=$number, type=$type\n";
	print "Warning: gene $gene has multiple names\n" if defined $data{$gene}{$file1};
	$data{$gene}{$file1} = $type;
	my @pos = `grep ",POS0," $file1`;
	my @neg = `grep ",NEG0," $file1`;
	$strand{$gene} = @pos > 0 ? "Pos" : @neg > 0 ? "Neg" : "Unk";
}
mkdir $opt_o if not -d $opt_o;
if (not -d $opt_o) {die "Cannot create directory $LCY$opt_o$N: $!\n"}
mkdir "$opt_o/TEMP" if not -d "$opt_o/TEMP";
if (not -d "$opt_o/TEMP") {die "Cannot create directory $LCY$opt_o/TEMP$N: $!\n"}


my @outputs;
my $exist = 0;
my @dir2 = `ls -Ra $dir2`;
my $dir2curr;
my %file2;
my $count = 0;
open (my $outshell, ">", "$opt_o/$labelz\_run.sh") or die;
print $outshell "#!/bin/bash -l\n#SBATCH -p high --mem 16000\n";
my $print = ""; my %skipped;
for (my $i = 0; $i < @dir2; $i++) {
	chomp($dir2[$i]);
	if ($dir2[$i] =~ /:$/) {
		$dir2curr = $dir2[$i]; $dir2curr =~ s/:$//;
	}
	else {
		my $file2 = $dir2curr . "/" . $dir2[$i];
		next if $file2 !~ /.PEAKS$/;
		# create directory similar to footpeak one
		my @dirs2 = split("/", $dir2curr);
		my $curroutdir = "$opt_o/TEMP/";
		for (my $k = 0; $k < @dirs2; $k++) {
			$curroutdir .= "$dirs2[$k]/";
		#	print "Making $curroutdir\n" if not -d $curroutdir;
			mkdir $curroutdir if not -d $curroutdir;
		}
		my ($label, $gene, $strand, $window, $threshold, $convType) = $dir2[$i] =~ /^(.+)_gene(.+)_(Pos|Neg|Unk)_(\d+)_(\d+\.?\d*)_(CH|GH|CG|GC)/;
		die "Cannot parse from file in footpeak $LCY$dir2[$i]$N\n" if not defined $threshold or not defined $gene;
		$gene = "NA" if not defined $gene;
		$gene = uc($gene);
		if (defined $opt_g) {
			next if $gene eq "NA";
			next if $gene ne $opt_g;
		}
		my $geneStrand = $strand{$gene};
		if (not defined $data{$gene} and not defined $skipped{$gene}) {
			$skipped{$gene} = 1;
			$print .= "gene $gene from FOOTPEAK is skipped as this gene does not have any Robert's predicted peaks (no bed file exists in $dir1)\n";
			next;
		}
		next if not defined $geneStrand;
		$strand = "NA" if not defined $strand;
		$window = "NA" if not defined $window;
		$threshold = "NA" if not defined $threshold;
		$convType = "NA" if not defined $convType;
		($convType) = $convType =~ /^(\w+)\./ if $convType =~ /\./;
		next if $geneStrand ne $strand; 
		next if $strand eq "Pos" and $convType !~ /^(CH|CG)$/;
		next if $strand eq "Neg" and $convType !~ /^(GH|GC)$/;
#		die "gene=$gene, strand=$strand, window=$window, threshold=$threshold, $convType = $convType\n";
		if (defined $opt_w) {
			next if $window eq "NA";
			next if $window != $opt_w;
		}
		if (defined $opt_t) {
			next if $threshold eq "NA";
			next if $threshold != $opt_t;
		}
		if (defined $opt_h) {
			next if $convType eq "NA";
			next if $convType !~ /^(GH|CH)$/;
		}
		if (defined $opt_c) {
			next if $convType eq "NA";
			next if $convType !~ /^(GC|CG)$/;
		}
		#print "$file2\n\tgene=$gene, strand=$strand, window=$window, thresold=$threshold, convType=$convType\n" if $file2 =~ /.PEAKS$/;
		foreach my $file1 (sort keys %{$data{$gene}}) {
			my $type = $data{$gene}{$file1};
			my ($folder1, $fileName1) = getFilename($file1, "folderfull");
			my ($folder2, $fileName2) = getFilename($file2, "folderfull");
			$fileName2 =~ s/.PEAK.PEAKS$//;
         my ($fileName22) = $file2 =~ /($fileName2.+)\.PEAK.PEAKS$/;
         $fileName2 = defined $fileName22 ? $fileName22 : $fileName2;
#         my $output = "$folder2/$fileName1\_$fileName2.STAT";
			my ($shuffledFolder, $shuffledFilename) = getFilename($file1, "folderfull");
			my ($footPeakFolder, $footPeakFilename) = getFilename($file2, "folderfull");
         my $output = "$curroutdir/$labelz\_$shuffledFilename\_$footPeakFilename.STAT";
			#print "2: $output\n";
			push(@outputs, $output);
			print $outshell "$script $file1 $file2 $opt_o $gene $type $strand $window $threshold $convType $curroutdir $labelz\n" if defined $opt_F or (not defined $opt_F and not -e $output) or (not defined $opt_F and -e $output and -s $output == 0);
#			system("$script $file1 $file2 $gene $type $strand $window $threshold $convType");
			$exist ++ if -e $output;
			$count ++;
		}
	}
}
print "\n" . $print . "\n";
open (my $out3, ">", "$opt_o/$labelz\_STAT.sh") or die;
print $out3 "#!/bin/bash\n";
my $good = 0;
for (my $i = 0; $i < @outputs; $i++) {
#	print "2b: $outputs[$i] does nto exist\n" if not -e $outputs[$i];
#	print "2c: $outputs[$i] Exist!\n" if -e $outputs[$i];
#	next if not -e $outputs[$i];
	print $outshell "cat $outputs[$i] > $opt_o/$labelz\_STAT.tsv\n" if $good == 0;
	print $outshell "grep -vP '^#' $outputs[$i] >> $opt_o/$labelz\_STAT.tsv\n" if $good != 0;
	print $out3 "cat $outputs[$i] > $opt_o/$labelz\_STAT.tsv\n" if $good == 0;
	print $out3 "grep -vP '^#' $outputs[$i] >> $opt_o/$labelz\_STAT.tsv\n" if $good != 0;
	$good ++;
}
close $out3;

if (not defined $opt_x) {
	my ($job) = `sbatch $opt_o/$labelz\_run.sh` =~ /job (\d+)$/;
	print "\nRunning: sbatch $opt_o/$labelz\_run.sh\nJob ID: $job\n\n";
}
else {
	print "\n-R is on: NOT RUNNING: sbatch $opt_o/$labelz\_run.sh\n\n";
}
if ($exist < $count or defined $opt_F) {
	my $remaining = $count - $exist;
	if (not defined $opt_F) {
		print "This will run all necessary intersects that haven't been run yet ($LGN$remaining$N/$LCY$count$N runs left)\n";
		print "Use -f when using this script to force run all intersects (doesn't matter if the output file exists or not)\n\n";
	}
}
else {
	print "All intersects ($LGN$exist$N/$LCY$count$N) has been run!\n";
}
print "

Statistics from $exist intersect files will be concatenated into$LCY $opt_o/$labelz\_STAT.tsv$N
If once $opt_o/$labelz\_run.sh is done but there's no $opt_o/$labelz\_STAT.tsv, then manually run $LCY$opt_o/$labelz\_STAT.sh$N to concatenate all intersect stats.\n\n";

#If total intersected is $exist < total=$count, then run the $opt_o.sh and it'll intersect whatever not exist yet
#Add -f when using this script if you want to redo everything (and not just intersect whatever not exist).
#Concatenated Output from $exist files: $opt_o

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
        return($fullname)        if defined($type) and $type =~ /fullname/;
        return($folder, $fullname, $shortname)  if defined($type) and $type =~ /all/;
}
__END__
my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");


my %data;
open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /^#/;
	my @arr = split("\t", $line);
}
close $in1;


open (my $out1, ">", "$fileName1.out") or die "Cannot write to $fileName1.out: $!\n";
close $out1;



__END__
		if (defined $opt_s) {
			next if $strand eq "NA";
			if ($opt_s eq "Rel") {
				next unless ($strand eq "Pos" and $convType = /^(CH|CG)$/) or ($strand eq "Neg" and $convType = /^(GH|GC)$/);
			}
			else {
				next if $strand ne $opt_s;
			}
		}
		if (defined $opt_c) {
			next if defined $opt_s and $opt_s eq "Rel";
			next if $convType eq "NA";
			next if $convType ne $opt_c;
		}

__END__
die "-s must be either Pos, Neg, Unk, or Rel (case sensitive)\n" if defined $opt_s and $opt_s !~ /^(Pos|Neg|Unk|Rel)$/;
die "-c must be either CH/CG/GH/GC (case sensitive)\n" if defined $opt_c and $opt_c !~ /^(CH|CG|GH|GC)$/;

-s (Pos/Neg/Unk/Rel) get files with this strand 
	-> Rel = get files with relevant strand and convType only
		which are Pos_CH, Pos_CG, Neg_GH, Neg_GC files.
		The reason being, C will be expected conversion in a Pos strand read and not G,
		so we only use Pos_CH and Pos_CG and not use Pos_GH and Pos_GC.
		Similarly, G will be expected conversion in a Neg strand read and not C
		so we only use Neg_GH and Neg_GC and not Neg_CH and Neg_CG.
-c (CH/CG/GH/GC) get files with this conversion type only (disabled if -s is Rel)
