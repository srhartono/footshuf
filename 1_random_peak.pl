#!/usr/bin/perl

use strict; use warnings; use Getopt::Std;
use vars qw($opt_v $opt_i $opt_r $opt_o $opt_b $opt_s $opt_f $opt_l);
getopts("vi:r:o:bs:fl:");

my ($N,$GN,$YW,$CY,$RD,$PR,$LGN,$LCY,$LRD) = ("\e[0m","\e[0;32m","\e[1;33m","\e[0;36m","\e[0;31m","\e[0;35m","\e[1;32m","\e[1;36m","\e[1;31m","\e[1;35m");

my ($input1, $label, $randNum, $outDir, $backupBoolean) = check_sanity();
my $seed_number = defined($opt_s) ? $opt_s : "NO SEED DEFINED";
if (defined $opt_s) {
	srand($seed_number);
}
print "
-s seed number   : $seed_number
-i Input         : $input1
-l Label         : $label
-r Random Number : $randNum
-o Output Dir    : $outDir
-b Backup files  : $backupBoolean;
";

my @files = <$outDir/*>;
=comment
if ($backupBoolean eq "TRUE" and -d $outDir) {
	my $backupDir = "$outDir/TEMPBackup" . int(rand(1000));
	makedir($backupDir);
	system("mv $outDir/*.* $backupDir") if @files > 0 or (@files == 1 and -d "$outDir/BackupDir");
	system("mv $outDir/BackupDir $backupDir") if -d "$outDir/BackupDir";
	system("mv $backupDir $outDir/BackupDir/");
}
=cut
# Gene count is gene in order of appearance in the input file. Will be put in the output name in case of same gene name in same file
my ($GENE, $GENE_COUNT) = (0,0); 
my ($peaks, $browserCoor);
my $logFile = "$input1\_logFile.txt";
open (my $outLog, ">", $logFile) or die "Cannot write to $logFile: $!\n";

LOG($outLog, "\n\n\e[0;33m ======================================== \e[0m\n\n");
LOG($outLog, "Parsing \e[1;36m$input1\e[0m\n\n");

my @input1 = split("/", $input1);
open (my $in1, "<", $input1) or LOG($outLog, "Failed to read from $input1: $!\n") and die;
open (my $outLabel, ">", "$outDir/.LABEL") or die "Failed to write to $outDir/.LABEL: $!\n";
print $outLabel $label;
close $outLabel;
my $linecount = 0;
while (my $line = <$in1>) {
	chomp($line);
	$linecount ++;
	my @arr = split("\t", $line);
	if ($line =~ /^browser posi/i) {
		$GENE_COUNT ++;
		my ($CHR0, $BEG0, $END0) = $line =~ /position (.+):([,\d]+)\-([,\d]+)$/i;
		my $line2 = $line;
		$line = <$in1>; chomp($line);
		$line = uc($line);
		($GENE) = $line =~ /^#(.+)$/;
		die "Cannot read gene name under line=$line2:\nline=$line\n" if not defined $GENE;
		@{$peaks->{$GENE_COUNT}{bcoor}} = ($CHR0, $BEG0, $END0);
		$peaks->{$GENE_COUNT}{strandc}{Pos} = 0;
		$peaks->{$GENE_COUNT}{strandc}{Neg} = 0;
#		$line = <$in1>; chomp($line);
#		LOG($outLog, "Line does not contain track_name! ($line)\n") if $line !~ /track name/;
#		($GENE) = $line =~ /description="HG19_(\w+(\-\w+)*)_\w+|DNA"/;
#		LOG($outLog, "Undefined gene at line=$line\n") and die if not defined $GENE;
	}
	elsif (@arr >= 3) {
		my ($chr, $beg, $end, $name, $zero, $strand) = split("\t", $line);
		die "Undefined col6 strand at line\n\n$line\n\n" if not defined $strand;
		my $strandz = $strand eq "+" ? "Pos" : $strand eq "-" ? "Neg" : die "Strand has to be either + or - (current: $strand)\n\n$line\n\n";
		$peaks->{$GENE_COUNT}{strandc}{$strandz} ++;
		my ($CHR0, $BEG0, $END0) = @{$peaks->{$GENE_COUNT}{bcoor}};
		LOG($outLog, "\n\n\e[0;31mERROR!!\e[0m: Undefined chr, beg, end, or name at line=$line\n") and die if not defined ($chr . $beg . $end . $name);
		LOG($outLog, "\n\n\e[0;31mERROR!!\e[0m: Current track's CHR is $LGN$CHR0$N  but current line's chr is $LCY$chr$N at line=$line\n") and die if $chr ne $CHR0;
		($beg, $end) = ($beg - $BEG0, $end - $BEG0);
		push(@{$peaks->{$GENE_COUNT}{coor}}, "$chr,$beg,$end,$name,$strand");
		$peaks->{$GENE_COUNT}{gene} = $GENE;
	}
}
close $in1;

LOG($outLog, "Parsed $LGN" . scalar(keys %{$peaks}) . "$N genes from -i $LCY$input1$N:\n");

foreach my $GENE_COUNT (sort {$a <=> $b} keys %{$peaks}) {
	my $GENE = $peaks->{$GENE_COUNT}{gene};
	my $total = @{$peaks->{$GENE_COUNT}{coor}};
	my ($CHR0, $BEG0, $END0) = @{$peaks->{$GENE_COUNT}{bcoor}};
	my $strandPos = $peaks->{$GENE_COUNT}{strandc}{Pos};
	my $strandNeg = $peaks->{$GENE_COUNT}{strandc}{Neg};
	$peaks->{$GENE_COUNT}{strand} = $strandPos >= $strandNeg ? "Pos" : "Neg";
	LOG($outLog, "\e[0;33m$GENE_COUNT\e[0m. \e[0;36m$GENE\e[0m (\e[0;32m$total\e[0m total peaks)\n");
	for (my $i = 0; $i < @{$peaks->{$GENE_COUNT}{coor}}; $i++) {
		my ($chr, $beg, $end, $name, $strand) = split(",", $peaks->{$GENE_COUNT}{coor}[$i]);
		$beg += $BEG0;
		$end += $BEG0;
		LOG($outLog, "  \e[0;32m$GENE_COUNT.$i\e[0m: chr=$chr,beg=$beg,end=$end,name=$name,strand=$strand\n");
	}
}

LOG($outLog, "\n\n\e[0;33m ======================================== \e[0m\n\n");

my $geneCount = 0;
foreach my $GENE_COUNT (sort {$a <=> $b} keys %{$peaks}) {
	my $gene = $peaks->{$GENE_COUNT}{gene};
	my $coor = $peaks->{$GENE_COUNT}{coor};
	my $bcoor = $peaks->{$GENE_COUNT}{bcoor};
	my $strand = $peaks->{$GENE_COUNT}{strand};
	next if not defined $coor;
	next if not defined $bcoor;
	next if not defined $gene;
	#next if $gene !~ /CALM3/;
	LOG($outLog, "$gene:\n\t- " . join("\n\t- ", @{$coor}) . "\n",1);
	LOG($outLog, "$geneCount. Shuffled $LGN" . scalar(@{$coor}) . "$N peaks $PR$randNum$N times from gene $LCY$gene$N\n");
	shuf($gene, $GENE_COUNT, $coor, $bcoor, $outDir, $outLog, $label, $strand);
	$geneCount ++;
}
LOG($outLog, "Log file: $logFile\n");

sub check_sanity {
	die "\nUsage: $YW$0$N $LCY-i$N <BED4 (4+ columns)> $PR-l$N <label [alphanumeric]> $LGN-s$N <seed number [0]> $LGN-r$N <number of randoms [def: 1000]> $PR-o$N <output dir>\n\n" unless defined $opt_i;
	die "\nError: -i input1 does not exist!\n" if not -e $opt_i;
	die "\nError: -i input1 is emptyt!\n" if -s $opt_i == 0;
	$opt_r = 1000 if not defined $opt_r;
	$opt_o = defined $opt_o ? $opt_o : "$opt_i\_RANDOM";
	die "\nError: -r number of random must be positive integer!\n" if $opt_r !~ /^\d+(e\d+)?$/;
	die "\nError: -o output dir $LCY$opt_o$N cannot be the same as input ($opt_i)\n" if $opt_o eq $opt_i;
	die "\nError: -o output dir $LCY$opt_o$N already exist! (use -f to override at your own risk)\n\n" if -d $opt_o and not defined $opt_f;
	makedir($opt_o) if not -d $opt_o;
	my $backupBoolean = defined $opt_b ? "TRUE" : "FALSE";
	die "Seed number must be positive integer! (currently: $LGN$opt_s$N)\n" if defined $opt_s and $opt_s !~ /^\d+$/;
	die "Label must be defined!\n" if not defined $opt_l;
	die "Label must be alphanumeric\n" if defined $opt_l and $opt_l !~ /^[A-Za-z0-9]+$/;
	return($opt_i, $opt_l, $opt_r, $opt_o, $backupBoolean);
}

sub LOG {
	my ($outLog, $text, $quiet) = @_;
	print STDERR $text if not defined $quiet;
	print $outLog $text;
}
sub shuf {
	my ($id, $id_count, $coor, $bcoor, $outDir, $outLog, $label, $strand) = @_;
	LOG($outLog, "Undefined \$coor or \$bcoor->{'GENE'} id=$id\n") and die if not defined $coor or not defined $bcoor;
	my ($CHR, $BEG0, $END0) = @{$bcoor};
	my ($gene, $type) = $id =~ /^(\w+)(\-\w)?$/;
	$type = "" if not defined $type; $type =~ s/\-//; $type = "_$type" if $type ne "";
	my $BEG = 1; my $END = $END0 - $BEG0; my $LEN0 = $END0 - $BEG0;
	LOG($outLog, "CHR=$CHR, $BEG-$END, BEG=$BEG0, END=$END0, LEN=$LEN0, GENE=$gene, TYPE=$type, GENE_COUNT=$id_count, BEG=$BEG, END=$END\n",1);
	open (my $out1, ">", "$outDir/$label\_$gene$type\_$id_count\_$strand.bed") or die "Cannot write to $gene.bed: $!\n";
	my $meanIter; my $print;
	for (my $i = 0; $i < $randNum; $i++) {
		($print, $meanIter) = shuffle_coor("$gene$type\_$id_count\_$strand", $coor, $CHR, $BEG0, $END0, $i, $meanIter, $outLog);
		print $out1 $print;
	}
	close $out1;
	$meanIter = mean($meanIter);
	LOG($outLog, "\tITER=$meanIter\n",1)
}
sub shuffle_coor {
	my ($gene, $coorArr, $CHR, $BEG0, $END0, $i, $meanIter, $outLog) = @_;
	my $BEG = 0; my $END = $END0 - $BEG0; my $LEN0 = $END0 - $BEG0;
	my %used; my $print0 = ""; my $print = "";
	my $coor;
	for (my $j = 0; $j < @{$coorArr}; $j++) {
		my ($chr, $begOrig, $endOrig, $name, $strand) = split(",", $coorArr->[$j]);
		die if defined $coor->{$name};
		$coor->{$name}{coor} = $coorArr->[$j];
		$coor->{$name}{len} = $endOrig - $begOrig;
	}
	while (1) {
		my $checkIter = 0;
		foreach my $name (sort {$coor->{$b}{len} <=> $coor->{$a}{len}} keys %{$coor}) {
			my ($chr, $begOrig, $endOrig, $name, $strand) = split(",", $coor->{$name}{coor});
			
			my $lenOrig = $endOrig - $begOrig;
			my $begShuf = int(rand($END-1-$lenOrig));
			my $endShuf = $begShuf + $lenOrig;
			LOG($outLog, "\t$chr,ORIG=$begOrig-$endOrig,$name,$strand,SHUF=$begShuf-$endShuf\n",1);
			my $check = 1; my $iter = 0;
			while ($check == 1) {
				$iter ++;
				# check if begShuf and endShuf intersect with any currently shuffled peak, then redo
				$check = 0;
				foreach my $coorUsed (sort {$used{$b}[5] <=> $used{$a}[5]} keys %used) {
					my ($begUsed, $endUsed) = ($used{$coorUsed}[3], $used{$coorUsed}[4]);
					$check = 1 if intersect($begUsed, $endUsed, $begShuf, $endShuf) == 1;
					last if $check == 1;
				}
				# If yes, then redo randoming begShuf and endShuf
				if ($check == 1) {
					$begShuf = int(rand($END-$lenOrig));
					$endShuf = $begShuf + $lenOrig;
				}
	
				# otherwise, store it
				else {
					@{$used{$name}} = ($name, $begOrig, $endOrig, $begShuf, $endShuf, $lenOrig);
					last;
				}
				last if $iter == 200;
			}
			push(@{$meanIter}, $iter);
			# if this goes on too much without result then just stop at iter 200 and restart
			if ($iter == 200) {
				%used = ();
				$checkIter = 1;
				last;
			}
		}
		last if $checkIter == 0;
	}
	foreach my $coorUsed (sort {$used{$a}[1] <=> $used{$b}[2]} keys %used) {
		my $name      = $used{$coorUsed}[0];
		my $begOrig   = $used{$coorUsed}[1] + $BEG0;
		my $endOrig   = $used{$coorUsed}[2] + $BEG0;
		my $begShuf   = $used{$coorUsed}[3] + $BEG0;
		my $endShuf   = $used{$coorUsed}[4] + $BEG0;
		my $countShuf = $i + 1;
		$print0 .= "$CHR\t$begOrig\t$endOrig\t$gene,$name,0\n" if $i == 0;
		$print .= "$CHR\t$begShuf\t$endShuf\t$gene,$name,$countShuf\n";
	}
#	print $print0 . $print; die;
	return ($print0 . $print, $meanIter);
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
      foreach my $key (keys %{$values}) {
#        print "key=$key, value=$values->{$key}, total=$total, meancurr=$mean\n";
         $mean += ($values->{$key} / $total);
      }
   }
   else { die "Cannot determine type of values=$values=\n";}
   return($mean);
}

sub shuffle {
	my ($value, $times) = @_;
	my @value = @{$value};
	$times = @value if not defined($times) or $times !~ /^\d+$/;
	for (my $i = 0; $i < $times; $i++) {
		my ($ran1, $ran2) = ( int(rand(@value)), int(rand(@value)) );
		my ($val1, $val2) = ( $value[$ran1]    , $value[$ran2]     );
		($value[$ran1], $value[$ran2]) = ($val2, $val1);
	}
	return(@value);
}

sub intersect {
my ($start1, $end1, $start2, $end2) = @_;
die "Died at mitochy::intersect: start1 ($start1) can't be bigger than end1 ($end1)\n" if $start1 > $end1;
die "Died at mitochy::intersect: start2 ($start2) can't be bigger than end2 ($end2)\n" if $start2 > $end2;
return(1) if ($start1 >= $start2 and $start1 <= $end2) or ($start2 >= $start1 and $start2 <= $end1);
return(0);
}

sub makedir {
   my ($folder, $isFile) = @_;
   my $log = "";
   #$folder = getFullpath($folder);
   my @folders = split("/", $folder);
   die if @folders == 0;
   my $curr_folder = "";
   $log .= "FOLDER $LCY$folder$N is already exist!\n" and return $log if -e $folder;
   for (my $i = 0; $i < @folders; $i++) {
      last if $i == @folders-1 and defined $isFile;
      $curr_folder .= "$folders[$i]/";
      next if $curr_folder =~ /^\/$/;
      next if $curr_folder =~ /^(\/|\/home\/)$/;
      $log .= "$i: Undefined folder to add. Current folder=$LGN$curr_folder$N\n" and return $log if not defined $folders[$i];
      next if -d $curr_folder or -e $curr_folder;
      system("mkdir $curr_folder") == 0 or die "Failed to mkdir $curr_folder: $!\n";
      $log .= "$i. $curr_folder\n";
   }
   $log .= "FOLDER $LGN$curr_folder$N is made!\n" if -e $curr_folder;
   return ($log);
}

__END__


my %data;
open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /^#/;
	my @arr = split("\t", $line);
}
close $in1;


close $out1;

