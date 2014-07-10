#!/usr/bin/perl -w

use strict "vars";

# Get arguments
(@ARGV == 0) || usage(1, "No argument accepted.\n");

my $baseDir = "/afs/cern.ch/user/r/robutti/Workspace/private/CMSSW_6_2_0/src/ERobutti/HoughTransChecks";

# Scan parameters
my $maxNBins = 200000000;
my @nBinsDefault = (8, 40, 64, 12, 25);
my $nBinSteps = 7;
my $maxBinScale = 2.;

# Create directory if it does not exist
my $scanDir = "$baseDir/ScanHTBinPars_2";
-d $scanDir || mkdir $scanDir;
system "rm $scanDir/*";

for (my $iPar = 0; $iPar < 5; $iPar++) {
  my @nBins = @nBinsDefault;
  for (my $iStep = 0; $iStep < $nBinSteps; $iStep++) {
#    my $part = $maxBinScale**((2.*$iStep + 1 - $nBinSteps)/($nBinSteps - 1));
#    print "$part\n";
    $nBins[$iPar] = int ($nBinsDefault[$iPar]*($maxBinScale**((2.*$iStep + 1 - $nBinSteps)/($nBinSteps - 1))));
    my $iJob = join '-', @nBins;
    my $logFile = "$scanDir/houghCheck2Steps_$iJob.log";
    open LOGFILE, ">$logFile";
    print LOGFILE "Running HoughCheck2Steps with parameters:\n".
                   "nBinsDoca = $nBins[0]\n".
                   "nBinsSqrtK = $nBins[1]\n".
                   "nBinsPhi = $nBins[2]\n".
                   "nBinsZ0 = $nBins[3]\n".
                   "nBinsEta = $nBins[4]\n".
                   "phiBinOverlap = 0.25\n".
                   "etaBinOverlap = 0.25\n";
    # Customize configuration file
    my $pyConfigFile = "$scanDir/HoughCheck2Steps_cff_$iJob.py";
    print LOGFILE "Configuration file: $pyConfigFile\n";
    open PYTEMPLATE, "<$baseDir/test/HoughCheck2Steps_scanTemplate_cff.py";
    open PYJOB, ">$pyConfigFile";
    while (my $line = <PYTEMPLATE>) {
      if ($line =~ /nBins.*\((.*)\)/) {
        $line =~ s/$1/$nBins[0], $nBins[1], $nBins[2], $nBins[3], $nBins[4]/;
      } elsif ($line =~ /phiBinOverlap.*\((.*)\)/) {
        $line =~ s/$1/0.25/;
      } elsif ($line =~ /etaBinOverlap.*\((.*)\)/) {
        $line =~ s/$1/0.25/;
      } elsif ($line =~ /outRootFile.*\((.*)\)/) {
        $line =~ s/$1/\'$scanDir\/houghCheck_2Steps_$iJob\.root\'/;
      }
      print PYJOB $line;
    }
    system "export PYCONFIGFILE=$pyConfigFile; bsub -q 1nd -o $logFile source $baseDir/test/runHoughCheck2Steps.sh";
  }
}


# Help
sub usage {
  my($exit, $message) = @_;

  print STDERR $message if defined $message;
  print STDERR <<INLINE_LITERAL_TEXT; #'

Usage: $0

INLINE_LITERAL_TEXT
#'
  exit($exit) if defined $exit;
}
