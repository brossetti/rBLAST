#!/usr/bin/perl

#==============================================================================#
#
# PROGRAM: FASTA Filter v.1.0
#
# DATE LAST MODIFIED:
#
# PROGRAMMER: Blair J. Rossetti
#
# SUMMARY: The purpose of this Perl program is to extract a number of FASTA
#          formatted sequences based on sequence identifiers (IDs) listed in a
#          separate text file from a FASTA file containing multiple FASTA
#          sequences. The output file is a FASTA formatted text file containing
#          the available sequences for all specified identifiers. This program
#          is initiated from the command prompt.
#
# INPUT: There are two input files required for this Perl program: a FASTA
#        formatted file containing the sequences to be extracted and a text file
#        containing the list of IDs for which sequence extraction is desired.
#        The ID file must contain identification strings that are found in the
#        FASTA sequence header lines as per the example below. All input files
#        must be placed in the \input subdirectory of rBLAST1.0. Here are
#        example input files:
#
#        Example ID File - humanEST_IDs.txt
#        _______________________________________________________________
#       |GH270327.1                                                     |
#       |BM251187.1                                                     |
#       |BM251188.1                                                     |
#       |BM251189.1                                                     |
#       |_______________________________________________________________|
#
#        Example FASTA File - humanEST_FASTA.txt
#        _______________________________________________________________
#       |>gi|217678545|gb|GH270327.1|CHA-V-97 Hs cDNA example sequence  |
#       |ACTTTTCCCCCCCCCCCCCCCCCCCCCCTTTCGGCAGTTCAGGTTTTAAGGGAGGGCAACAGA|
#       |ACAATACTAGTAGAAACAGAATATGAAATGAAAAGTGGGGAAAGGCTTCATGAGTATGAATTT|
#       |ATTTCAATGTTACTTAGGCTGAATAAAGACAATGATCATTAAGTAATTCTAGAAAATAGAGAT|
#       |CCTTATAAGATTCCTTTGATTCCTGTTAATTGCAATCCAATACATAAATATTGAAGGCTTATT|
#       |TTTCGTATGACAAAGATGAAAATACTTCTGCCGGTTAACATTTTATAATTCTATGACCTTATG|
#       |AAGATATTTTTAAACCTTGACAGTTACATTGGGGAGATTTACATCAAGAACCAGCACTATGGA|
#       |GTTCTATCTTCAAAGTTTTACTACAAACCGGTGCCACTGGGTCATTTTTAGATTCAGACCATG|
#       |CCACAGTGAACCGGACAGTCATGTGCCACTTGGCAGTTTTCAAGGTATGAGTTTATTCTGTTT|
#       |TCCAGGGGACAAGCAGCCAGATGGTTAATTGCCCTTCCTTAAGCACAGACGCACTCACAGAAA|
#       |CCCTCCGAAAGAGGAGCCTCCACAGAGAGGTGACGTTTTTTTTTTTTG               |
#       |>gi|217678546|gb|GH270328.1|CHA-V-165 Hs cDNA example sequence |
#       |ACGAACACACANGTTTCCCAAAACCACAACTACTACCTTTAGAGTTCAGGCTACTTTCACATA|
#       |AGAGGACCCAAAACAACAAGACATCTGTGGAGGGAGGAGAACAAAACAGAGAATAACTTAGAC|
#       |TCTAGAAGCTATCTTCTACAGAAAAACATCTCCTTTCCAAATTTTCCCTTCTTTCAGCATAAT|
#       |AACTACTTTATTGAGGTATACTTGATACACAGGAAACTGCACATATTTAAAGCATGCAATTTG|
#       |TATCCAACCATGAAACAGTCCCCACCAGCAAGGTGATGTGCACACCCCTCATCCCTCAAAGTC|
#       |CATCCCTCCCTTCCACCACACTCTATGCCTCATCTCCAGCAAACACTGAGCTGCTTCCAATTT|
#       |TTTTAGAAAATGTTTCAAATCAGAAACTAAAGTCTATACTGAACATTATTGTTTTACATGTGC|
#       |GCTCCTCTAAACCAATCGTCATCCTGTCTCTAAAGTGATTAGAAATACTGGAAGTCGAATACT|
#       |ATGCTCCTGCAGAGTTGGCTATGTCAAGGCTGCATGCACAGCCTGTTCAGGAAAAGTGGATGG|
#       |GCACATCTCCAGCTTCCCAGGCAATCCGCTGGAGGCAAGCTACACACATTCTGAAAATAAATG|
#       |>gi|218927127|gb|BM251187.1| EST20021 Hs cDNA example sequence |
#       |ACTCTGAGGCTTGTAGGAGGGTAAAATAGAGACCCAGTAAAATTGTAATAAGCAGTGCTTGAA|
#       |TTTTCTATTAGACTATGGTGAGCTCAGGTGATTGATACTCCTGATGCGAGTAATACGGATGTG|
#       |TAGGGGATTTAGCGGGGTGATGCCTGTTGGGGGCCAGTGCCCTCCTAATTGGGGGGTAGGGGC|
#       |AGGCTCAGAAAAATCCTGCGAAGAAAAAAACTTCTGAGGTAATAAATAGGATTATCCCGTATC|
#       |GGTGGTGTGTGGTGGCCTTGGTATGTGCTTTCTCGTGTTACATCGCGCCATCATTGGTATATG|
#       |TAGGCCTAGTATGAGGAGCGTTATGGAGTGGAAGTGAAATCACATGGCTAGGCCGGAGGTCAT|
#       |CCCCTGTTAGGGGTCATGGGCTGGGTTTTACTATATGATAG                      |
#       |>gi|218927128|gb|BM251188.1| EST10048 Hs cDNA example sequence |
#       |ACTGCAATAACAAAATACAGCAATAAAACAACTGGACACTCCTAGGGGACACCAAAGATAAAG|
#       |TAGGCCAGAGAAACCCAACCTGTTGGCAATATGACGCTCTTTCCCAACTGGGTCTTGGTGAGA|
#       |CTGTCAGTGCATGTGCATAAATTGTAGACCAGGTCCCACTATGCTACTTCAGGATTCAGCCAG|
#       |GAGGTCCCTTGGTCCTTATTCATCTTGATATACTCATGGGATGTTTGGAATTAAGGAGCCCAA|
#       |TGAAGTCTTTCACTTCATATCCTACCCCATCAGTCTAAGAGCCCACCCAACAAGGGTAGCTAC|
#       |CTATAGGCTGCCTGATCCTGGAACCAGCCTGGGGCCCTGATTATGATCTCCACGGGGCTGTCA|
#       |TTTAAGACCCAGCG                                                 |
#       |>gi|218927129|gb|BM251189.1| EST20019 Hs cDNA example sequence |
#       |ACTGATCATTCTATTTCCCCCTCTATTGATCCCCACCTCCAAATATCTCATCAACAACCGACT|
#       |GACTAATCAAACTAACCTCAAAACAAATGATAACCATACACAACACTAAAGGACGAACCTGAT|
#       |TTAATCATTTTTATTGCCACAACTAACCTCCTCGGACTCCTGCCTCACTCATTTACACCAACC|
#       |CCTAGCCATGGCCATCCCCTTATGAGCGGGCGCAGTGATTATAGGCTTTCGCTCTAAGATTAA|
#       |TCTTACCACAAGGCACACCTACACCCCTTATCCCCATACTAGTTATTATCGAAACCATCAGCC|
#       |GCCCTGGCCGT                                                    |
#       |_______________________________________________________________|
#
# COMMAND PROMPT: This Perl program is initiated and controlled from the
#                 command prompt or terminal using command line arguments. To
#                 run this Perl program, open the command prompt or terminal and
#                 change prompt to the directory containing this Perl program
#                 (C:\...\rBLAST1.0\bin). The following are the command line
#                 arguments and an example command prompt entry:
#
#        fastaFilter 1.0 arguments:
#
#           -i   ID file name [File In]
#           -f   FASTA file name [File In]
#           -o   Output file name [File Out]
#           -m   Allow multiple FASTA sequences for each specified ID [T/F]
#             default = F
#
#        Example Command Prompt -
#        _______________________________________________________________
#       |                                                               |
#       |C:\rBlast1.0\bin>perl fastaFilter.pl -i humanEST_IDs.txt       |
#       |-f humanEST_FASTA.txt -o humanEST_output.txt -m F              |
#       |_______________________________________________________________|
#
#       The example command prompt entry shown above initiates the Sequence
#       Extractor v.1.0 program and specifies four arguments. The
#       "-i humanEST_IDs.txt" argument tells the program that the file
#       "humanEST_IDs.txt" contains the FASTA IDs for which sequences will be
#       extracted. The "-f humanEST_FASTA.txt" argument tells the program that
#       the file humanEST_FASTA.txt contains the FASTA sequences. The
#       "-o humanEST_output.txt" argument is a user-defined file name for the
#       output FASTA file. The "-m F" argument tells the program that there
#       exists only one corresponding sequence for every ID listed in
#       "humanEST_IDs.txt".
#
# OUTPUT: The output for this Perl program is a FASTA formatted text file
#         containing the FASTA formatted sequences specified by the ID file. All
#         output files are saved to the \output subdirectory of rBLAST1.0.
#         Here is an example output file:
#
#        Example Output File - humanEST_output.txt
#        _______________________________________________________________
#       |>gi|217678545|gb|GH270327.1| CHA-V-97 Hs cDNA example sequence |
#       |ACTTTTCCCCCCCCCCCCCCCCCCCCCCTTTCGGCAGTTCAGGTTTTAAGGGAGGGCAACAGA|
#       |ACAATACTAGTAGAAACAGAATATGAAATGAAAAGTGGGGAAAGGCTTCATGAGTATGAATTT|
#       |ATTTCAATGTTACTTAGGCTGAATAAAGACAATGATCATTAAGTAATTCTAGAAAATAGAGAT|
#       |CCTTATAAGATTCCTTTGATTCCTGTTAATTGCAATCCAATACATAAATATTGAAGGCTTATT|
#       |TTTCGTATGACAAAGATGAAAATACTTCTGCCGGTTAACATTTTATAATTCTATGACCTTATG|
#       |AAGATATTTTTAAACCTTGACAGTTACATTGGGGAGATTTACATCAAGAACCAGCACTATGGA|
#       |GTTCTATCTTCAAAGTTTTACTACAAACCGGTGCCACTGGGTCATTTTTAGATTCAGACCATG|
#       |CCACAGTGAACCGGACAGTCATGTGCCACTTGGCAGTTTTCAAGGTATGAGTTTATTCTGTTT|
#       |TCCAGGGGACAAGCAGCCAGATGGTTAATTGCCCTTCCTTAAGCACAGACGCACTCACAGAAA|
#       |CCCTCCGAAAGAGGAGCCTCCACAGAGAGGTGACGTTTTTTTTTTTTG               |
#       |>gi|218927127|gb|BM251187.1| EST20021 Hs cDNA example sequence |
#       |ACTCTGAGGCTTGTAGGAGGGTAAAATAGAGACCCAGTAAAATTGTAATAAGCAGTGCTTGAA|
#       |TTTTCTATTAGACTATGGTGAGCTCAGGTGATTGATACTCCTGATGCGAGTAATACGGATGTG|
#       |TAGGGGATTTAGCGGGGTGATGCCTGTTGGGGGCCAGTGCCCTCCTAATTGGGGGGTAGGGGC|
#       |AGGCTCAGAAAAATCCTGCGAAGAAAAAAACTTCTGAGGTAATAAATAGGATTATCCCGTATC|
#       |GGTGGTGTGTGGTGGCCTTGGTATGTGCTTTCTCGTGTTACATCGCGCCATCATTGGTATATG|
#       |TAGGCCTAGTATGAGGAGCGTTATGGAGTGGAAGTGAAATCACATGGCTAGGCCGGAGGTCAT|
#       |CCCCTGTTAGGGGTCATGGGCTGGGTTTTACTATATGATAG                      |
#       |>gi|218927128|gb|BM251188.1| EST10048 Hs cDNA example sequence |
#       |ACTGCAATAACAAAATACAGCAATAAAACAACTGGACACTCCTAGGGGACACCAAAGATAAAG|
#       |TAGGCCAGAGAAACCCAACCTGTTGGCAATATGACGCTCTTTCCCAACTGGGTCTTGGTGAGA|
#       |CTGTCAGTGCATGTGCATAAATTGTAGACCAGGTCCCACTATGCTACTTCAGGATTCAGCCAG|
#       |GAGGTCCCTTGGTCCTTATTCATCTTGATATACTCATGGGATGTTTGGAATTAAGGAGCCCAA|
#       |TGAAGTCTTTCACTTCATATCCTACCCCATCAGTCTAAGAGCCCACCCAACAAGGGTAGCTAC|
#       |CTATAGGCTGCCTGATCCTGGAACCAGCCTGGGGCCCTGATTATGATCTCCACGGGGCTGTCA|
#       |TTTAAGACCCAGCG                                                 |
#       |>gi|218927129|gb|BM251189.1| EST20019 Hs cDNA example sequence |
#       |ACTGATCATTCTATTTCCCCCTCTATTGATCCCCACCTCCAAATATCTCATCAACAACCGACT|
#       |GACTAATCAAACTAACCTCAAAACAAATGATAACCATACACAACACTAAAGGACGAACCTGAT|
#       |TTAATCATTTTTATTGCCACAACTAACCTCCTCGGACTCCTGCCTCACTCATTTACACCAACC|
#       |CCTAGCCATGGCCATCCCCTTATGAGCGGGCGCAGTGATTATAGGCTTTCGCTCTAAGATTAA|
#       |TCTTACCACAAGGCACACCTACACCCCTTATCCCCATACTAGTTATTATCGAAACCATCAGCC|
#       |GCCCTGGCCGT                                                    |
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
my @notFound;         #Will hold all sequence identifiers not found in FASTA
                      # file

#----Variables----
my $startTime;        #Will hold time for determining elapsed time of program
my $elapsedTime;      #Will hold elapsed time of program
my $idFile;           #Will hold file name for input file containing
                      # list of IDs to be extracted
my $fastaFile;        #Will hold file name for input FASTA file containing
                      # the sequences requiring extraction
my $outputFile;       #Will hold user-defined file name for output file
my $multiOption;      #Will hold true or false option for allowing multiple
                      # FASTA sequences for each specified ID
my $commandLine;      #Will hold user input from command line
my $fastaPath;        #Will hold relative path for input FASTA file
my $idPath;           #Will hold relative path for input ID file
my $outputPath;       #Will hold relative path for output file
my $numIDs;           #Will hold the total number of IDs in ID file
my $numSeqs;          #Will hold the total number of sequence printed to output
my $numMissing;       #Will hold the total number of sequences in ID file that
                      # that were not found in FASTA file

#_______________________________________________________________________________
#Section 1: Retrieve Command Line Options

$startTime = time;                     #Start elapsed timer

$commandLine = "@ARGV";                #Place command line arguments in variable

#-----------------------
#Determine if command line contains all necessary arguments
if (!($commandLine =~ m/-{1,2}i(dFile)?/i)
    or !($commandLine =~ m/-{1,2}f(astaFile)?/i)
    or !($commandLine =~ m/-{1,2}o(utput)?/i))
{
    print "\n\nfastaFilter 1.0 arguments:\n\n";

    print "   -i   ID file name [File In]\n";
    print "   -f   FASTA file name [File In]\n";
    print "   -o   Output file name [File Out]\n";
    print "   -m   Allow multiple FASTA sequences for each specified ID [T/F]";
    print "\n     default = F\n\n\n";
    print "   Example Command Prompt Arguments:\n\n";
    print "   C:\\rBlast1.0\\bin>perl fastaFilter.pl -i humanEST_IDs.txt\n";
    print "   -f humanEST_FASTA.txt -o humanEST_output.txt -m F\n\n";
    exit;
}

$multiOption = "false";                      #Set default to false

GetOptions('fastaFile|f=s' => \$fastaFile,   #Retrieve command line arguments
           'idFile|i=s' => \$idFile,
           'output|o=s' => \$outputFile,
           'multiple|m=s' => \$multiOption);

#-----------------------
#Set Multiple Sequence option
if ($multiOption =~ m/f(alse)?/i)
{
    $multiOption = "false";
}
elsif ($multiOption =~ m/t(rue)?/i)
{
    $multiOption = "true";
}

#_______________________________________________________________________________
#Section 2: Set Relative Paths

$fastaPath = File::Spec->catfile("../input", $fastaFile);     #Create relative
$idPath = File::Spec->catfile("../input", $idFile);           # paths for input
$outputPath = File::Spec->catfile("../output", $outputFile);  # and output files

#_______________________________________________________________________________
#Section 3: Check Files

checkFile($fastaPath, $idPath);              #Call checkFile subroutine to check
                                             # existance and readablity of files

#_______________________________________________________________________________
#Section 4: Sequence Extraction

($numIDs, $numSeqs, $numMissing, @notFound) = filter($fastaPath,    #Call filter
                                                     $idPath,       # subroutine
                                                     $outputPath,   # to extract
                                                     $multiOption); # target
                                                                    # sequences
                                                                    # based on
                                                                    # IDs

$elapsedTime = time - $startTime;                       #Determine elapsed time

print "\n\nProcess Complete...\n\n\n";
print "FASTA Filter 1.0 Console Report:\n\n";
print "* Number of Requested FASTA Sequences: $numIDs\n";
print "* Number of Sequences Found: $numSeqs\n";
print "* Number of Sequences Not Found: $numMissing\n";

#-----------------------
#Determine if missing sequences exist
if (!($numMissing == 0))
{
    print "* IDs With Missing Sequences:\n";
    print @notFound;
}

print "\n\nTime Elapsed: $elapsedTime sec\n";

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


#-------\
# filter \
#-------------------------------------------------------------------------------
# SUMMARY: Extracts FASTA formatted sequences from a text file containing
#          multiple FASTA sequences. FASTA sequences are extracted based on
#          sequence identifiers (IDs) listed in a text file. The output file
#          is a FASTA formatted text file containing the available sequences
#          for all specified identifiers.
#
# ARGUMENTS: (0) FASTA file name or path, (1) ID file name or path, (2) Output
#            file name or path, (3) Multiple Sequences option
#
# RETURNS: (0) Total number of request sequences, (1) Total number of sequences
#          found, (2) Total number of sequences not found, (3) Array of IDs
#          for which no sequence was found
#
sub filter
{
  #_________________________________
  # Data Structure Declaration List:

  #-----Arrays------
  my @ids;            #Will hold all sequence identifiers
  my @notFound;       #Will hold all sequence identifiers not found in FASTA
                      # file

  #----Variables----
  my $fastaPath;      #Will hold relative path for input FASTA file
  my $idPath;         #Will hold relative path for input ID file
  my $outputPath;     #Will hold relative path for output file
  my $multiOption;    #Will hold true or false option for allowing multiple
                      # FASTA sequences for each specified ID
  my $numIDs;         #Will hold the total number of IDs in ID file
  my $id;             #Will hold each sequence idenfitier in @ids
  my $line;           #Will hold each line of FASTA file
  my $numSeqs;        #Will hold the total number of sequence printed to output
  my $numMissing;     #Will hold the total number of sequences in ID file that
                      # that were not found in FASTA file
  my $i;              #Will hold index number
  my $j;              #Will hold index number

  #_________________________________
  #Section 1: Retrieve Arguments

  $fastaPath = shift;                          #Shift off first argument
  $idPath = shift;                             #Shift off second argument
  $outputPath = shift;                         #Shift off third argument
  $multiOption = shift;                        #Shift off fourth argument

  #_________________________________
  #Section 2: Open/Create Files

  open(FASTA, "<", "$fastaPath")                    #Open file for reading
      or die "\nError: Unable to open $fastaFile [$!]\n";

  open(IDFILE, "<", "$idPath")                      #Open file for reading
      or die "\nError: Unable to open $idFile [$!]\n";

  open(OUTPUT, ">", "$outputPath")                  #Open file for writing
      or die "\nError: Unable to open $outputFile [$!]\n";

  #_________________________________
  #Section 3: Retrieve ID Data

  @ids = <IDFILE>;      #Place all IDs into an array
  close (IDFILE);       #Close ID file

  #_________________________________
  #Section 4: Extract FASTA Sequences

  $numIDs = @ids;                                 #Determine the number of IDs
  $numSeqs = 0;                                   #Initialize variable
  $numMissing = 0;                                #Initialize variable

  #-----------------------
  #Prepare each sequence ID and search
  foreach $id (@ids)
  {
      $i = 0;                                     #Initialize index
      $j = 1;                                     #Initialize index
      $line = "";                                 #Initialize sequence

      chomp($id);                                 #Remove return character

      seek(FASTA, 0, 0);                          #Set/Reset read pointer

      #-----------------------
      #Search each line for ID
      while ($line = <FASTA>)
      {
          #-----------------------
          #Determine if line matches ID
          if (($line =~ m/^>.*$id/) and ($i == 0))
          {
              print OUTPUT "$line";
              $i = 1;                             #Set index to 1
              $j = 0;                             #Set index to 0
              $numSeqs++;                         #Increase number of sequences
                                                  # by one unit
          }
          elsif (($line =~ m/^>.*$id/) and ($i == 1))
          {
              #-----------------------
              #Stop search after first sequence match
              if ($multiOption eq "false")
              {
                  last;                           #Move to next ID
              }

              print OUTPUT "$line";

              $numSeqs++;                         #Increase number of sequences
                                                  # by one unit
          }
          elsif (($line =~ m/^>/) and ($i == 1))  #Match start of next sequence
          {
              $i = 0;                             #Set index to 0

             #-----------------------
             #Stop search after first sequence match
             if ($multiOption eq "false")
             {
                 last;                            #Move to next ID
             }
          }
          elsif ($i == 1)
          {
              print OUTPUT "$line";
          }
      }#End of while loop

      #-----------------------
      #Determine if ID was not found
      if ($j == 1)
      {
          $notFound[$numMissing] = "$id\n";       #Add missing ID to array
          $numMissing++;                          #Increase number of missing
                                                  # missing sequences
      }

  }#End of foreach loop

  close (FASTA);                                  #Close FASTA file
  close (OUTPUT);                                 #Close output file

  #_________________________________
  #Section 5: Return

  return($numIDs, $numSeqs, $numMissing, @notFound);

}
#_________________________________End_of_filter________________________________#