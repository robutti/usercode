#!/usr/bin/perl -w

use strict "vars";

# Get arguments
(@ARGV == 0) || usage(1, "No argument accepted.\n");

my $baseDir = "/afs/cern.ch/user/r/robutti/Workspace/private/CMSSW_6_1_0/src/HoughTest/HoughTransChecks";

# Scan parameters
my $maxNBins = 200000000;
my @nBinsDoca = (4, 8, 16);
my @nBinsSqrtK = (20, 40, 80);
my @nBinsPhi = (64, 128, 256);
my @nBinsZ0 = (3, 6, 12);
my @nBinsEta = (25, 50, 100);
my @phiBinOverlap = (0.25);
my @etaBinOverlap = (0.25);

# Create directory if it does not exist
my $scanDir = "$baseDir/ScanHTBinPars";
-d $scanDir || mkdir $scanDir;
system "rm $scanDir/*";

my @pBin = (0, 0, 0, 0, 0, 0, 0);
foreach my $nDoca (@nBinsDoca) {
  $pBin[0]++;
  @pBin[1..6] = (0, 0, 0, 0, 0, 0);
  foreach my $nSqrtK (@nBinsSqrtK) {
    $pBin[1]++;
    @pBin[2..6] = (0, 0, 0, 0, 0);
    foreach my $nPhi (@nBinsPhi) {
      $pBin[2]++;
      @pBin[3..6] = (0, 0, 0, 0);
      foreach my $nZ0 (@nBinsZ0) {
        $pBin[3]++;
        @pBin[4..6] = (0, 0, 0);
        foreach my $nEta (@nBinsEta) {
          $pBin[4]++;
          @pBin[5..6] = (0, 0);
          next if ($nDoca*$nSqrtK*$nPhi*$nZ0*$nEta > $maxNBins);
          foreach my $phiOverlap (@phiBinOverlap) {
            $pBin[5]++;
            $pBin[6] = 0;
            foreach my $etaOverlap (@etaBinOverlap) {
              $pBin[6]++;
              my $iJob = join '', @pBin;
              my $logFile = "$scanDir/houghCheck2Steps_$iJob.log";
              print "Running HoughCheck2Steps with parameters:\n".
                    "nBinsDoca = $nDoca\n".
                    "nBinsSqrtK = $nSqrtK\n".
                    "nBinsPhi = $nPhi\n".
                    "nBinsZ0 = $nZ0\n".
                    "nBinsEta = $nEta\n".
                    "phiBinOverlap = $phiOverlap\n".
                    "etaBinOverlap = $etaOverlap\n";
              # Customize configuration file
              my $pyConfigFile = "$scanDir/HoughCheck2Steps_cff_$iJob.py";
              print "Configuration file: $pyConfigFile\n";
              open PYTEMPLATE, "<$baseDir/test/HoughCheck2Steps_scanTemplate_cff.py";
              open PYJOB, ">$pyConfigFile";
              while (my $line = <PYTEMPLATE>) {
                if ($line =~ /nBins.*\((.*)\)/) {
                  $line =~ s/$1/$nDoca, $nSqrtK, $nPhi, $nZ0, $nEta/;
                } elsif ($line =~ /phiBinOverlap.*\((.*)\)/) {
                  $line =~ s/$1/$phiOverlap/;
                } elsif ($line =~ /etaBinOverlap.*\((.*)\)/) {
                  $line =~ s/$1/$etaOverlap/;
                } elsif ($line =~ /outRootFile.*\((.*)\)/) {
                  $line =~ s/$1/\'$scanDir\/houghCheck_2Steps_$iJob\.root\'/;
                }
                print PYJOB $line;
              }
              system "export PYCONFIGFILE=$pyConfigFile; bsub -q 1nd -o $logFile source /afs/cern.ch/user/r/robutti/Workspace/private/CMSSW_6_1_0/src/HoughTest/HoughTransChecks/test/runHoughCheck2Steps.sh";
            }
          }
        }
      }
    }
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
