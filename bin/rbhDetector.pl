#!/usr/bin/perl

#==============================================================================#
#
# PROGRAM: Reciprocal Best Hit Detector v.1.0
#
# DATE LAST MODIFIED:
#
# PROGRAMMER: Blair J. Rossetti
#
# SUMMARY: The purpose of this Perl program is to determine reciprocal best
#          BLAST hits from forward and reverse BLAST outputs. Any ID may have
#          multiple reciprocal best BLAST hits. A reciprocal best hit is defined
#          as a set of two sequences that are each other's best hit in forward
#          and reverse BLASTs. Reciprocal BLAST hits are printed to a user-
#          defined output file. This program is initiated from the command
#          prompt.
#
# INPUT: There are two input files required for this Perl program: a tabular
#        BLAST output of best hits for a forward BLAST and a tabular BLAST
#        output of best hits for a reverse BLAST. All input files must be placed
#        in the \input subdirectory of rBLAST1.0. Here are example input files:
#
#        Example Forward Tabular BLAST File - forwardBlast.txt
#        _______________________________________________________________
#       |query_001 id03 37.07 259 148 5 138 382 1574 1831 2e-037 154    |
#       |query_003 id57 45.75 365 174 9 1   343 1    363  3e-074 275    |
#       |query_003 id22 45.75 365 174 9 1   343 1    363  3e-074 275    |
#       |query_005 id38 50.23 240 138 4 16  397 1403 1213 3e-021 232    |
#       |_______________________________________________________________|
#
#        Example Reverse Tabular BLAST File - reverseBlast.txt
#        _______________________________________________________________
#       |id03 query_001 37.07 259 148 5 138 382 1554 1831 2e-034 153    |
#       |id57 query_003 34.56 352 173 8 3   321 1    343  3e-072 272    |
#       |id22 query_003 43.43 312 173 7 2   354 1    341  3e-070 269    |
#       |id38 query_001 50.23 259 138 4 16  382 176  1213 3e-017 143    |
#       |_______________________________________________________________|
#
#
# COMMAND PROMPT: This Perl program is initiated and controlled from the
#                 command prompt or terminal using command line arguments. To
#                 run this Perl program, open the command prompt or terminal and
#                 change prompt to the directory containing this Perl program
#                 (C:\...\rBLAST1.0\bin). The following are the command line
#                 arguments and an example command prompt entry:
#
#        rbhDetector 1.0 arguments:
#
#           -f   Forward tabular BLAST-output file name [File In]
#           -r   Reverse tabular BLAST-output file name [File In]
#           -o   Output file name [File Out]
#
#        Example Command Prompt -
#        _______________________________________________________________
#       |                                                               |
#       |C:\rBlast1.0\bin> perl rbhDetector.pl -f forwardBlast.txt      |
#       |-r reverseBlast.txt -o rbhOutput.txt                           |
#       |_______________________________________________________________|
#
#       The example command prompt entry shown above initiates the Reciprocal
#       Best Hit Filter v.1.0 program and specifies three arguments. The
#       "-f forwardBlast.txt" argument tells the program which file contains
#       the forward BLAST output of best hits. The "-r reverseBlast.txt"
#       argument tells the program which file contains the reverse BLAST output
#       of best hits. The "-o rbhOutput.txt" argument is a user-defined file
#       name for the output reciprocal BLAST file.
#
# OUTPUT: The output for this Perl program is a tab-delimited formatted text
#         file containing the best reciprocal BLAST hits. All output files are
#         saved to the \output subdirectory of rBLAST1.0. Here is an example
#         output file:
#
#        Example Output File - rbhOutput.txt
#        _______________________________________________________________
#       |query_001 id03 2e-037 154 2e-034 153                           |
#       |query_003 id57 3e-074 275 3e-072 272                           |
#       |query_003 id22 3e-074 275 3e-070 269                           |
#       |_______________________________________________________________|
#
#        Tabular Output Columns:
#           (0)  Forward query name
#           (1)  Reverse query name
#           (2)  Forward E-value
#           (3)  Forward bit score
#           (4)  Reverse E-value
#           (5)  Reverse bit score
#
#==============================================================================#

#******************************************************************************#
#---------------------------PROGRAM CODE BEGINS HERE---------------------------#
#******************************************************************************#

#_______________________________________________________________________________
# Pragma and Module List:

use strict;           #Pragma to restrict unsafe constructs
use warnings;         #Pragma to control optional warnings
use Getopt::Long;     #Module to extend processing of command line options
use File::Spec;       #Module to portably perform operations on file names

#_______________________________________________________________________________
# Data Structure Declaration List:

#------Array------
my @multiIDs;         #Will hold the query IDs of all queries with multiple RBHs

#----Variables----
my $commandLine;      #Will hold user input from command line
my $inputOne;         #Will hold file name for forward input file containing
                      # tabular BLAST best hit data for forward BLAST
my $inputTwo;         #Will hold file name for reverse input file containing
                      # tabular BLAST best hit data for reverse BLAST
my $outputFile;       #Will hold user-defined file name for output file
my $forwardPath;      #Will hold relative path for forward BLAST file
my $reversePath;      #Will hold relative path for reverse BLAST file
my $outputPath;       #Will hold relative path for output file
my $hitCount;         #Will hold the total number of reciprocal best hits
my $multiCount;       #Will hold the number of queries with multiple RBHs


#_______________________________________________________________________________
#Section 1: Retrieve Command Line Options

$commandLine = "@ARGV";                #Place command line arguments in variable

#-----------------------
#Determine if command line contains all necessary arguments
if (!($commandLine =~ m/-{1,2}f(orward)?/i)
    or !($commandLine =~ m/-{1,2}r(everse)?/i)
    or !($commandLine =~ m/-{1,2}o(utput)?/i))
{
    print "\n\nrbhDetector 1.0 arguments:\n\n";

    print "   -f   Forward tabular BLAST output file name [File In]\n";
    print "   -r   Reverse tabular BLAST output file name [File In]\n";
    print "   -o   Output file name [File Out]\n\n\n";
    print "   Example Command Prompt Arguments:\n\n";
    print "   C:\\rBlast1.0\\bin>perl rbhDetector.pl -f forwardBlast.txt\n";
    print "   -r reverseBlast.txt -o rbhOutput.txt\n\n";
    exit;
}

GetOptions('forward|f=s' => \$inputOne,   #Retrieve command line arguments
           'reverse|r=s' => \$inputTwo,
           'output|o=s' => \$outputFile);

#_______________________________________________________________________________
#Section 2: Set Relative Paths

$forwardPath = File::Spec->catfile("../input", $inputOne);    #Create relative
$reversePath = File::Spec->catfile("../input", $inputTwo);    # paths for input
$outputPath = File::Spec->catfile("../output", $outputFile);  # and output files

#_______________________________________________________________________________
#Section 3: File Checking

checkFile($forwardPath, $reversePath);       #Call checkFile subroutine to check
                                             # existance and readablity of files

#_______________________________________________________________________________
#Section 4: RBH Detection

                                             #Call rBlastFilter subroutine to
                                             # determine reciprocal best hits
($hitCount, $multiCount, @multiIDs) = rbhDetect($forwardPath,
                                                $reversePath,
                                                $outputPath);


print "\n\nProcess Complete...\n\n\n";
print "Reciprocal Best Hit Detector 1.0 Console Report:\n\n";
print "* Number of RBHs = $hitCount\n";
print "* Number of Multiple RBHs = $multiCount\n";

#-----------------------
#Determine if multiple best hits exist
if (!($multiCount == 0))
{
    print "* IDs with Multiple RBH:\n";
    print @multiIDs;
}

print "\n";

#_______________________________________________________________________________
#===============================================================================
#Section 5: Subroutines

#----------\
# checkFile \
#-------------------------------------------------------------------------------
# SUMMARY: Checks files for existance and readability.
#
# ARGUMENTS: One or more file name or path.
#
# RETURNS: Nothing.
#
sub checkFile
{
  #_________________________________
  # Data Structure Declaration List:

  #------Array------
  my @files;          #Will hold all files to be checked

  #----Variables----
  my $file;           #Will hold each file in @files

  #_________________________________
  #Section 1: Retrieve File

  @files = @_;                            #Place arguments into an array

  #_________________________________
  #Section 2: Check File

  #-----------------------
  #Examine one file at a time
  foreach $file (@files)
  {
      #-----------------------
      #Determine if variable is undefined
      if (!defined($file))
      {
          print "\nError: One or more input files were not assigned\n";
          exit;
      }
      #-----------------------
      #Determine if file exists
      if (!(-e $file))
      {
          print "\nError: $file does not exist\n";
          exit;
      }
      #-----------------------
      #Determine if file is readable
      if (!(-r $file))
      {
          print "\nError: $file is not readable\n";
          exit;
      }
  }#End of foreach loop
}
#_______________________________End_of_checkFile_______________________________#



#----------\
# rbhDetect \
#-------------------------------------------------------------------------------
# SUMMARY: Determines reciprocal best BLAST hits from forward and reverse BLAST
#          outputs. Any ID may have multiple reciprocal best BLAST hits. A
#          reciprocal best hit is defined as a set of two sequences that are
#          each others best hit in forward and reverse BLASTs. Reciprocal BLAST
#          hits are printed to a user-defined output file.
#
# ARGUMENTS: (0) Forward tabular BLAST file name or path, (1) Reverse tabular
#            BLAST file name or path, (2) Output file name or path
#
# RETURNS: The total number of best reciprocal BLAST hits
#
sub rbhDetect
{
  #_________________________________
  # Data Structure Declaration List:

  #------Array------
  my @reverse;        #Will hold every line from tabular reverse BLAST output
  my @blastLineF;     #Will hold each column from a forward BLAST line
  my @blastLineR;     #Will hold each column from a reverse BLAST line
  my @multiIDs;

  #----Variables----
  my $forwardPath;    #Will hold relative path for forward BLAST file
  my $reversePath;    #Will hold relative path for reverse BLAST file
  my $outputPath;     #Will hold relative path for output file
  my $queryColumn;    #Will hold column position of query ID in tabular BLAST
                      # output
  my $subjectColumn;  #Will hold column position of subject ID in tabular BLAST
                      # output
  my $evalColumn;     #Will hold column position of E-values in tabular BLAST
                      # output
  my $scoreColumn;    #Will hold column position of bit score in tabular BLAST
                      # output
  my $queryF;         #Will hold forward BLAST query ID
  my $queryR;         #Will hold reverse BLAST query ID
  my $subjectF;       #Will hold forward BLAST subject ID
  my $subjectR;       #Will hold reverse BLAST subject ID
  my $evalF;          #Will hold forward BLAST E-value
  my $evalR;          #Will hold reverse BLAST E-value
  my $scoreF;         #Will hold forward BLAST bit score
  my $scoreR;         #Will hold reverse BLAST bit score
  my $forwardLine;    #Will hold each line from tabular forward BLAST output
  my $reverseLine;    #Will hold each line from tabular reverse BLAST output
  my $hitCount;       #Will hold the total number of reciprocal best hits
  my $multiCount;     #Will hold the number of queries with multiple RBHs
  my $lastQueryF;     #Will hold the previous forward BLAST query ID
  my $lastQueryR;     #Will hold the previous reverse BLAST query ID
  my $i;              #Will hold index number

  #_________________________________
  #Section 1: Set Constants

  $scoreColumn = 11;          #Set bit score column position
  $evalColumn = 10;           #Set E-value column position
  $subjectColumn = 1;         #Set subject column position
  $queryColumn = 0;           #Set query column position

  #_________________________________
  #Section 2: Retrieve File

  $forwardPath = shift;            #Shift off first argument
  $reversePath = shift;            #Shift off second argument
  $outputPath = shift;             #Shift off third argument

  #_________________________________
  #Section 3: Open/Create Files

  open(FORWARD, "<", "$forwardPath")                    #Open file for reading
      or die "Error: Unable to open $forwardPath [$!]\n";

  open(REVERSE, "<", "$reversePath")                    #Open file for reading
      or die "Error: Unable to open $reversePath [$!]\n";

  open(OUTPUT, ">", "$outputPath")                      #Open file for writing
      or die "Error: Unable to open $outputPath: [$!]\n";

  #_________________________________
  #Section 4: Determine RBH

  @reverse = <REVERSE>;                             #Place reverse BLAST line
                                                    # into array

  $hitCount = 0;                                    #Initialize variable
  $multiCount = 0;                                  #Initialize variable
  $lastQueryF = "";                                 #Initialize variable
  $lastQueryR = "";                                 #Initialize variable

  #-----------------------
  #Examine each line of forward BLAST file
  while($forwardLine = <FORWARD>)
  {
      chomp($forwardLine);                          #Remove return character

      @blastLineF = split("\t", $forwardLine);      #Separate BLAST line into
                                                    # array

      $queryF = $blastLineF[$queryColumn];          #Extract forward query ID
      $subjectF = $blastLineF[$subjectColumn];      #Extract forward subject ID
      $evalF = $blastLineF[$evalColumn];            #Extract forward E-value
      $scoreF = $blastLineF[$scoreColumn];          #Extract forward bit score

      #-----------------------
      #Examine each line of reverse BLAST file
      foreach $reverseLine (@reverse)
      {

          chomp($reverseLine);                      #Remove return character

          @blastLineR = split("\t", $reverseLine);  #Separate BLAST line into
                                                    # array

          $queryR = $blastLineR[$queryColumn];      #Extract reverse query ID
          $subjectR = $blastLineR[$subjectColumn];  #Extract reverse subject ID
          $evalR = $blastLineR[$evalColumn];        #Extract reverse E-value
          $scoreR = $blastLineR[$scoreColumn];      #Extract reverse bit score

          #-----------------------
          #Compare forward and reverse BLAST hits
          if(($queryF eq $subjectR) and ($subjectF eq $queryR))
          {
              chomp($scoreF);                       #Remove return character
              chomp($scoreR);                       #Remove return character

              print OUTPUT "$queryF\t$queryR\t$evalF\t$scoreF\t$evalR\t$scoreR\n";

              #-----------------------
              #Determine if multiple RBHs exist
              if (($queryF eq $lastQueryF) and ($i == 0))
              {
                  $multiIDs[$multiCount] = "$queryF\n";
                  $multiCount++;
                  $i = 1;
              }
              elsif (($queryR eq $lastQueryR) and ($i == 0))
              {
                  $multiIDs[$multiCount] = "$queryR\n";
                  $multiCount++;
                  $i = 1;
              }
              else
              {
                  $i = 0;
              }

              $lastQueryR = $queryR;


              $hitCount++;                          #Increase hit count by one
                                                    # unit

              last;                                 #End search after finding
                                                    # match
          }

          @blastLineR = ();                         #Clear array

      }#End of foreach loop

      $lastQueryF = $queryF;

      @blastLineF =();                              #Clear array

  }#End of while loop

  #_________________________________
  #Section 5: Return

  return($hitCount, $multiCount, @multiIDs);

}
#_______________________________End_of_rbhDetect_______________________________#