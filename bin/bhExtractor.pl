#!/usr/bin/perl

#==============================================================================#
#
# PROGRAM: Best Hit Extractor v.1.0
#
# PROGRAMMER: Blair J. Rossetti
#
# DATE LAST MODIFIED:
#
# SUMMARY: The purpose of this Perl program is to determine the best BLAST
#          hit(s) for each BLAST query sequence when given a list of multiple
#          BLAST hits. Best BLAST hits are determined first by highest bit score
#          and subsequently by lowest E-value. Multiple best hits are provided
#          as output for each query sequence if those hits share both identical
#          bit score and identical E-value. A best BLAST hit must have an
#          E-value lower than a user-defined E-value cutoff (default = 1E-10).
#          Best BLAST hits are printed to a user-defined output file. The user
#          can also invoke an optional function to print a file of hit IDs. The
#          hit ID file is formatted for use in Sequence Extractor. This program
#          is initiated from the command prompt.
#
# INPUT: There is one input file required for this Perl program: a tabular
#        (-m 8) BLAST output file. The tabular BLAST output file is created by
#        BLAST when the -m 8 option is specified during a BLAST run. The
#        tabular BLAST output file contains twelve tab-separated columns
#        containing the following data:
#
#        Tabular BLAST Output Columns:
#           (0)  Query name
#           (1)  Subject name
#           (2)  Percent identities
#           (3)  Aligned length
#           (4)  Number of mismatched positions
#           (5)  Number of gap positions
#           (6)  Query sequence start
#           (7)  Query sequence end
#           (8)  Subject sequence start
#           (9)  Subject sequence end
#           (10) E-value
#           (11) Bit score
#
#        All input files must be placed in the \input subdirectory of rBLAST1.0.
#        Here is an example input file:
#
#        Example Tabular BLAST File - blastOutput.txt
#        _______________________________________________________________
#       |query_001 id03 37.07 259 148 5 138 382 1574 1831 2e-037 154    |
#       |query_001 id43 33.33 270 151 7 143 385 336  603  2e-033 140    |
#       |query_001 id32 34.75 282 146 7 139 382 876  1157 2e-032 137    |
#       |query_001 id96 36.60 265 149 8 137 386 195  455  7e-029 125    |                                                |
#       |query_002 id05 22.97 209 121 6 160 361 245  420  6e-005 45.8   |
#       |query_003 id57 45.75 365 174 9 1   343 1    363  3e-074 275    |
#       |query_003 id22 45.75 365 174 9 1   343 1    363  3e-074 275    |
#       |query_003 id79 33.65 312 147 4 5   306 121  382  5e-038 154    |
#       |_______________________________________________________________|
#
# COMMAND PROMPT: This Perl program is initiated and controlled from the
#                 command prompt or terminal using command line arguments. To
#                 run this Perl program, open the command prompt or terminal and
#                 change prompt to the directory containing this Perl program
#                 (C:\...\rBLAST1.0\bin). The following are the command line
#                 arguments and an example command prompt entry:
#
#        bhExtractor 1.0 arguments:
#
#           -i   Tabular BLAST-output file name [File In]
#           -o   Output file name [File Out]
#           -e   Expectation value (E) [Real]
#             default = 1E-10
#           -b   Best BLAST hit ID File [File Out] Optional
#
#        Example Command Prompt -
#        _______________________________________________________________
#       |                                                               |
#       |C:\rBlast1.0\bin> perl bhExtractor.pl -i blastOutput.txt       |
#       |-o bestHits.txt -e 1E-6 -b hitIDs.txt                          |
#       |_______________________________________________________________|
#
#       The example command prompt entry shown above initiates the Best Hit
#       Finder v.1.0 program and specifies four arguments. The
#       "-i blastOutput.txt" argument tells the program which file contains
#       the tabular BLAST data. The "-o bestHits.txt" argument is a user-
#       defined file name for the output best hit file. The "-e 1E-6" argument
#       tells the program that the user would like the maximum E-value for a
#       best hit is 1E-6. The "-b hitIDs.txt" is an optional argument that tells
#       the program to create a file, with the name "hitIDs.txt", containing the
#       hit IDs for all best hits.
#
# OUTPUT: There are two output files for this Perl program: a file containing
#         the best BLAST hits in tabular BLAST output format and an optional
#         text file containing the hit IDs for all best hits. All  output files
#         are saved to the \output subdirectory of rBLAST1.0. Here are example
#         output files:
#
#        Example Output File - bestHits.txt
#        _______________________________________________________________
#       |query_001 id03 37.07 259 148 5 138 382 1574 1831 2e-037 154    |
#       |query_003 id57 45.75 365 174 9 1   343 1    363  3e-074 275    |
#       |query_003 id22 45.75 365 174 9 1   343 1    363  3e-074 275    |
#       |_______________________________________________________________|
#
#        Example Hit ID File - hitIDs.txt
#        _______________________________________________________________
#       |id03                                                           |
#       |id57                                                           |
#       |id22                                                           |
#       |_______________________________________________________________|
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

#-----Arrays------
my @multiIDs;         #Will hold the query IDs for those with multiple best hits

#----Variables----
my $blastFile;        #Will hold file name for tabular BLAST output file
my $outputFile;       #Will hold user-defined file name for output file
my $eval;             #Will hold cutoff E-value
my $idFile;           #Will hold fule name for best BLAST hit ID list file
my $idFileOption;     #Will hold true or false option if an ID file name has
                      # been defined by the user
my $commandLine;      #Will hold user input from command line
my $blastPath;        #Will hold relative path for tabular BLAST file
my $outputPath;       #Will hold relative path for output file
my $idPath;           #Will hold relative path for best BLAST hit ID list file
my $hitCount;         #Will hold the total number of best hits
my $numBlastHits;     #Will hold the total number of hits in the BLAST file
my $multiCount;       #Will hold the total number of queries with multiple best
                      # hits

#_______________________________________________________________________________
#Section 1: Retrieve Command Line Options

$commandLine = "@ARGV";                #Place command line arguments in variable

#-----------------------
#Determine if command line contains all necessary arguments
if (!($commandLine =~ m/-{1,2}i(nput)?/i)
    or !($commandLine =~ m/-{1,2}o(utput)?/i))
{
    print "\n\nbhExtractor 1.0 arguments:\n\n";

    print "   -i   Tabular BLAST output file name [File In]\n";
    print "   -o   Output file name [File Out]\n";
    print "   -e   Expectation value (E) [Real]\n";
    print "     default = 1E-10\n";
    print "   -b   Best BLAST hit ID File [File Out] Optional\n\n\n";
    print "   Example Command Prompt Arguments:\n\n";
    print "   C:\\rBlast1.0\\bin>perl bhExtractor.pl -i blastOutput.txt\n";
    print "   -o bestHits.txt -e 1E-6 -b hitIDs.txt\n\n";
    exit;
}

$eval = 1E-10;                           #Set default cutoff E-value

GetOptions('input|i=s' => \$blastFile,   #Retrieve command line arguments
           'output|o=s' => \$outputFile,
           'evalue|e=f' => \$eval,
           'bestHit|b=s' => \$idFile);

#_______________________________________________________________________________
#Section 2: Set Relative Paths

$blastPath = File::Spec->catfile("../input", $blastFile);     #Create relative
$outputPath = File::Spec->catfile("../output", $outputFile);  # paths for input
                                                              # and output files
#-----------------------
#Set ID file relative path
if (defined($idFile))
{
    $idPath = File::Spec->catfile("../output", $idFile);
    $idFileOption = "true";                                   #Set ID file
}                                                             # option to true
else
{
    $idPath = "false";
    $idFileOption = "false";
}

#_______________________________________________________________________________
#Section 3: Check Files

checkFile($blastPath);                       #Call checkFile subroutine to check
                                             # existance and readablity of files

#_______________________________________________________________________________
#Section 4: Find Best BLAST Hits

                                             #Call bhExtract subroutine to
                                             # determine best BLAST hits
($numBlastHits, $hitCount, $multiCount, @multiIDs) = bhExtract($blastPath,
                                                               $eval,
                                                               $idPath,
                                                               $idFileOption,
                                                               $outputPath);


print "\n\nProcess Complete...\n\n\n";
print "Best Hit Extractor 1.0 Console Report:\n\n";
print "* Number of Hits in BLAST Output: $numBlastHits\n";
print "* Number of Best Hits: $hitCount\n";
print "* Number of Queries with Multiple Best Hits = $multiCount\n";

#-----------------------
#Determine if multiple best hits exist
if (!($multiCount == 0))
{
    print "* Queries with Multiple Best Hits:\n";
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
# bhExtract \
#-------------------------------------------------------------------------------
# SUMMARY: Determines the best BLAST hit(s) for each query sequence listed in a
#          tabular BLAST output file. Multiple best hits are allowed for each
#          query sequence. Best hits are first determined by highest bit score
#          and subsequently by lowest E-value. Only those best hits with an E-
#          value lower than a specified cutoff value are printed to the output
#          text file. If the ID file option is specified as true, an ID file
#          contain the IDs of the best hits will be created. This subroutine
#          require three other subroutines: blastStatSort, compare, and
#          printBestHits.
#
# ARGUMENTS: (0) BLAST input file name or path, (1) Expectation value cutoff,
#            (2) ID file name or path, (3) ID file option, (4) Output file name
#            or path
#
# RETURNS: (0) The number of hits in BLAST file, (1) The total number of best
#          hits
#
sub bhExtract
{
  #_________________________________
  # Data Structure Declaration List:

  #-----Arrays------
  my @line;           #Will hold each column of tabular BLAST output
  my @blastLines;     #Will hold tabular BLAST lines
  my @sortedLines;    #Will hold sorted BLAST hits
  my @multiIDs;       #Will hold the query IDs for those with multiple best hits

  #----Variables----
  my $queryColumn;    #Will hold column position of query ID in tabular BLAST
                      # output
  my $evalColumn;     #Will hold column position of E-values in tabular BLAST
                      # output
  my $blastPath;      #Will hold relative path for tabular BLAST file
  my $outputPath;     #Will hold relative path for output file
  my $eval;           #Will hold cutoff E-value
  my $idPath;         #Will hold relative path for best BLAST hit ID list
                      # file
  my $idFileOption;   #Will hold true or false option if an ID file name
                      # has been defined by the user
  my $line;           #Will hold each line of tabular BLAST file
  my $i;              #Will hold index number
  my $lastQuery;      #Will hold BLAST ID of previous query
  my $query;          #Will hold BLAST ID of query
  my $numBlastHits;   #Will hold the total number of hits in the BLAST file
  my $hitCount;       #Will hold the total number of best hits
  my $hitCountSub;    #Will hold a subset of the total number of best hits
  my $multiCount;     #Will hold the total number of queries with multiple best
                      # hits
  my $multiCountSub;  #Will hold a subset of the total number of queries with
                      # multiple best hits
  my $multiID;        #Will hold the query ID for those with multiple best hits

  #_________________________________
  #Section 1: Set Constants

  $queryColumn = 0;           #Set query column position
  $evalColumn = 10;           #Set E-value column position

  #_________________________________
  #Section 2: Retrieve Arguments

  $blastPath = shift;                    #Shift off first argument
  $eval = shift;                         #Shift off second argument
  $idPath = shift;                       #Shift off third argument
  $idFileOption = shift;                 #Shit off fourth argument
  $outputPath = shift;                   #Shift off fifth argument

  #_________________________________
  #Section 3: Open/Create Files

  open(BLAST, "<", "$blastPath")                    #Open file for reading
      or die "Error: Unable to open $blastPath [$!]\n";

  open(OUTPUT, ">", "$outputPath")                  #Open file for writing
      or die "Error: Unable to open $outputPath: [$!]\n";

  #-----------------------
  #Create ID file
  if ($idFileOption eq "true")
  {
      open(IDFILE, ">", "$idPath")                  #Open file for writing
          or die "Error: Unable to open $idPath [$!]\n";
  }

  #_________________________________
  #Section 4: Build and Sort Best Hit Arrays

  $line = "";                                  #Initialize variable
  $lastQuery = "";                             #Initialize variable
  $numBlastHits = 0;                           #Initialize variable
  $hitCount = 0;                               #initialize variable
  $multiCount = 0;                             #initialize variable
  $i = 0;                                      #Set index to 0

  #-----------------------
  #Compare each line of BLAST file
  while ($line = <BLAST>)
  {

      $numBlastHits++;                         #Increase number of blast hits by
                                               # one unit

      #-----------------------
      #Skip blank line in input
      if (length($line) == 1)
      {
          next;
      }

      chomp($line);                            #Remove return character

      @line = split("\t", $line);              #Separate BLAST line into array

      #-----------------------
      #Skip line if E-value does not meet cutoff
      if ($line[$evalColumn] > $eval)
      {
          next;
      }

      $query = $line[$queryColumn];            #Extract query ID from BLAST
                                               # line array                                            ,

      #-----------------------
      #Build array of BLAST lines for each query
      if (($lastQuery eq $query) or ($lastQuery eq ""))
      {
          $blastLines[$i] = $line;             #Place BLAST line in array
          $i++;                                #Increase index by one
          $lastQuery = $query;
      }
      else
      {
	  @sortedLines = blastStatSort(@blastLines); #Call blastStatSort
                                                     # subroutine to sort BLAST
                                                     # lines by bit score and
                                                     # E-values

                                                     #Call bhPrint subroutine to
                                                     # print the best BLAST hits
          ($hitCountSub, $multiCountSub, $multiID) = bhPrint($outputPath,
                                                             $idPath,
                                                             $idFileOption,
                                                             @sortedLines);


          $hitCount = $hitCount + $hitCountSub;      #Add subset to total

          @blastLines = ();                          #Clear array
          @sortedLines = ();                         #Clear array
          $blastLines[0] = "$line";                  #Place line in new array
          $i = 1;                                    #Set index to 1
          $lastQuery = $query;
      }

  }#End of while loop

  @sortedLines = blastStatSort(@blastLines);         #Call blastStatSort
                                                     # subroutine to sort BLAST
                                                     # lines by bit score and
                                                     # E-values

                                                     #Call bhPrint subroutine to
                                                     # print the best BLAST hits
  ($hitCountSub, $multiCountSub, $multiID) = bhPrint($outputPath,
                                                     $idPath,
                                                     $idFileOption,
                                                     @sortedLines);

  $multiIDs[$multiCount] = "$multiID\n";             #Place query ID in array

  $hitCount = $hitCount + $hitCountSub;              #Add subset to total
  $multiCount = $multiCount + $multiCountSub;        #Add subset to total

  #-----------------------
  #Determine if IDFILE is open
  if ($idFileOption eq "true")
  {
      close(IDFILE);
  }

  close(BLAST);
  close(OUTPUT);

  #_________________________________
  #Section 5: Return

  return($numBlastHits, $hitCount, $multiCount, @multiIDs);

}
#_______________________________End_of_bhExtract_______________________________#



#--------------\
# blastStatSort \
#-------------------------------------------------------------------------------
# SUMMARY: Sort tabular BLAST output lines based on highest bit score and
#          subsequently on lowest E-value.
#
# ARGUMENTS: (0) Array of tabular BLAST output lines
#
# RETURNS: Sorted array of tabular BLAST output lines
#
sub blastStatSort
{
  #_________________________________
  # Data Structure Declaration List:

  #-----Arrays------
  my @unsortedHits;   #Will hold unsorted BLAST hits
  my @sortedHits;     #Will hold sorted BLAST hits

  #_________________________________
  #Section 1: Retrieve Arguments

  @unsortedHits = @_;                         #Place arguments into an array

  #_________________________________
  #Section 2: Sort BLAST hits

  @sortedHits = sort compare(@unsortedHits);  #Sort BLAST hits using the compare
                                              # subroutine

  #_________________________________
  #Section 3: Return

  return(@sortedHits)

}
#_____________________________End_of_blastStatSort_____________________________#



#--------\
# compare \
#-------------------------------------------------------------------------------
# SUMMARY: Compares the bit scores and E-values in an array of tabular BLAST
#          output lines and sorts the array.
#
# ARGUMENTS: (0) Array of unsorted tabular BLAST output lines
#
# RETURNS: A sorted array of tabular BLAST lines based on highest bit score
#          and subsequently E-values.
#
sub compare
{
  #_________________________________
  # Data Structure Declaration List:

  #-----Arrays------
  my @A;              #Will hold the first tabular BLAST line
  my @B;              #Will hold the second tabular BLAST line

  #----Variables----
  my $evalColumn;     #Will hold column position of E-values in tabular BLAST
                      # output
  my $scoreColumn;    #Will hold column position of bit score in tabular BLAST
                      # output
  my $aScore;         #Will hold the bit score of first hit
  my $bScore;         #Will hold the bit score of second hit
  my $aEval;          #Will hold the E-value of first hit
  my $bEval;          #Will hold the E-value of second hit

  #_________________________________
  #Section 1: Set Constants

  $scoreColumn = 11;          #Set bit score column position
  $evalColumn = 10;           #Set E-value column position

  #_________________________________
  #Section 2: Sort BLAST hits

  @A = split(" ", $a);        #Separate first BLAST line into columns
  @B = split(" ", $b);        #Separate second BLAST line into columns

  $aScore = $A[$scoreColumn]; #Extract bit score for first hit
  $bScore = $B[$scoreColumn]; #Extract bit score for second hit

  $aEval = $A[$evalColumn];   #Extract E-value for first hit
  $bEval = $B[$evalColumn];   #Extract E-value for second hit

  #_________________________________
  #Section 3: Return

  return ($bScore <=> $aScore || $aEval <=> $bEval);   #Sort by highest bit
                                                       # score then by lowest
                                                       # E-value
}
#________________________________End_of_compare________________________________#



#--------\
# bhPrint \
#-------------------------------------------------------------------------------
# SUMMARY: Determines which BLAST hits are the best from a sorted list. Results
#          are printed to a user-defined output file. The user can also invoke
#          an optional function to print a file of hit IDs. The hit ID file is
#          prepared for use in Sequence Extractor.
#
# ARGUMENTS: (0) Output file name or path, (1) ID file name or path,
#            (2) ID file option, (3) Sorted tabular BLAST line array
#
#
# RETURNS: (0) The total number of best hits printed, (1) The total number of
#          query with multiple best hits (2) An array of the query IDs with
#          multiple best hits
#

sub bhPrint
{
  #_________________________________
  # Data Structure Declaration List:

  #-----Arrays------
  my @sortedLines;    #Will hold sorted BLAST hits
  my @columns;        #Will hold the each column of a tabular BLAST line

  #----Variables----
  my $evalColumn;     #Will hold column position of E-values in tabular BLAST
                      # output
  my $scoreColumn;    #Will hold column position of bit score in tabular BLAST
                      # output
  my $hitColumn;      #Will hold column position of hit ID in tabular BLAST
                      # output
  my $queryColumn;     #Will hold column position of query ID in tabular BLAST
                      # output
  my $outputPath;     #Will hold relative path for output file
  my $idPath;         #Will hold relative path for best BLAST hit ID list
                      # file
  my $idFileOption;   #Will hold true or false option if an ID file name
                      # has been defined by the user
  my $line;           #Will hold each line of tabular BLAST file
  my $i;              #Will hold index number
  my $score;          #Will hold bit score of each hit
  my $eval;           #Will hold E-value of each hit
  my $hitID;          #Will hold hit ID of each hit
  my $queryID;        #Will hold query ID of each hit
  my $lastScore;      #Will hold bit score of previous hit
  my $lastEval;       #Will hold E-value of previous hit
  my $lastQuery;      #Will hold query ID of prevous hit
  my $hitCount;       #Will hold the total number of best hits
  my $multiCount;     #Will hold the total number of queries with multiple best
                      # hits
  my $multiID;        #Will hold query ID with multiple best hits

  #_________________________________
  #Section 1: Set Constants

  $scoreColumn = 11;          #Set bit score column position
  $evalColumn = 10;           #Set E-value column position
  $hitColumn = 1;             #Set hit ID column position
  $queryColumn = 0;           #Set query column position

  #_________________________________
  #Section 2: Retrieve Arguments

  ($outputPath, $idPath, $idFileOption, @sortedLines) = @_;  #Place arguments
                                                             # into variables
                                                             # and an array

  #_________________________________
  #Section 3: Print Best Hits

  $i = 0;                                       #Set index to 0
  $lastScore = 0;                               #Initialize variable
  $lastEval = 1;                                #Initialize variable
  $lastQuery = "";                              #Initialize variable
  $hitCount = 0;                                #Initialize variable
  $multiCount = 0;                              #Initialize variable

  #-----------------------
  #Determine best hits
  foreach $line (@sortedLines)
  {
      @columns = split("\t", $line);            #Separate BLAST line into array
      $score = $columns[$scoreColumn];          #Extract score from BLAST line
      $eval = $columns[$evalColumn];            #Extract E-value from BLAST line
      $hitID = $columns[$hitColumn];            #Extract hit ID from BLAST line
      $queryID = $columns[$queryColumn];        #Extract query ID from BLAST
                                                # line

      #-----------------------
      #Print only best hits
      if ($score >= $lastScore and $eval <= $lastEval)
      {
          $hitCount++;                          #Increase hit count by one unit

          print OUTPUT "$sortedLines[$i]\n";

          #-----------------------
          #Determine if IDFILE is open
          if ($idFileOption eq "true")
          {
              print IDFILE "$hitID\n";
          }

          #-----------------------
          #Determine if query has multiple best hits
          if ($queryID eq $lastQuery)
          {
              $multiCount = 1;                  #Increase multiple best hit
                                                # count to 1
              $multiID = $queryID;
          }

          $i++;                                 #Increase index by one
          $lastScore = $score;
          $lastEval = $eval;
          $lastQuery = $queryID;
      }
      else
      {
          last;
      }

  }#End of foreach loop

  #_________________________________
  #Section 4: Return

  return($hitCount, $multiCount, $multiID);

}
#________________________________End_of_bhPrint________________________________#