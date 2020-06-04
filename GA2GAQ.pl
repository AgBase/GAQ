#!/usr/bin/perl -w
##################################################################################
#  Teresia Buza, Cathy Gresham, Nan Wang
#  10/2009
#
#  Script Name: GA2GAQ.pl
#  Purpose: Translate a gene association file, gene_ontology obo file 
#          and an evidence rank file to GAQ score.
#
# 11/09/2010 CRG Change to 17 column and 18 column data file (annot_cross and gaf)
# 09/30/2011 CRG Change to use a flat file to pull GO Depth data since it is not easy to
#    install Gavin Sherlock's GO::TermFinder / GO::Node in a windows enviroment
#
#  Input:
#     files: 	gene_association file 
#  		gene ontology obo file
#  		evidence rank file
#    flags: ignore obsolete go ids
#           shortest or longest path to go_term
#  Output:
#     files:  GAQ summary file      
#             GAQ detail file      
#             GAQ error file (possibly)
#
#  Root Terms
# GO:0008150 : biological_process 
# GO:0005575 : cellular_component 
# GO:0003674 : molecular_function 
#
# 08/15/2012 CRG getting invalid substitution in uninit values $goid and $evd line 281. wrap defined stmt
#         Initial cause of this error is 19 column GAF file (input id and GO Term Name) Handle these files.
#
##################################################################################

use strict;
use warnings;
use File::stat;
use Time::localtime;
use File::Temp qw/ tempfile /;
use Getopt::Std;
use IO::File;

#######################
# Global variables
#######################
my %GO = ();
my ($pathLength, $oboFile, $rootNode, $n );

my $outfile = "goDepth.txt";
my $currDate = time();

my %GOName = ();
my %GOParents = ();
my %Ranks = ();
my %GODepths = ();
my %GOObs= ();
my $ex_obs = 1;

my $obofile = "ontology/go-basic.obo";
my $depthfile = "ontology/goDepth.txt";
my $evdfile = "ontology/evidence_rank.txt";
#
# write to out_depth is commented out in script but can
# be used to show the depth that was calculated for
# each GO ID
#
my $out_summary = "GAQ_Summary.txt";
my $out_detail = "GAQ_Detail.txt";
my $out_error = "GAQ_Error.txt";
my $input_file = "";
my $results={};

our ( $opt_h, $opt_o, $opt_i, $opt_g, $opt_b, $opt_s, $opt_e, $opt_r);
##########################
# subroutine definitions
##########################
sub sub_buildGO;
sub sub_getDepth;
sub sub_checkOpts;
sub sub_getRank;
sub sub_parseGA($);
sub sub_writeGAQ($$);

####################################
# Begin program
####################################

getopts('hbo:i:g:s:e:r:');

my $printhelp = defined($opt_h);

# check the passed in options
&sub_checkOpts;

#  build the GO, the Depth and Rank hashes
&sub_buildGO;
&sub_getDepth;
&sub_getRank;

# Parse the input gene association file
$results->{gaq_info} = &sub_parseGA($input_file);
if (!$results) {
  die "Parser error. Could not process input file\n. Exiting now.\n\n";
}

# Write the output detail and summary files
&sub_writeGAQ($results, $ex_obs);

#################################################################
#
# sub_checkOpts
#
# Desc: this subroutine will read in user supplied options
#
#################################################################
sub sub_checkOpts {
  if ($printhelp) {
     print STDERR <<END;
Usage:  perl $0 [-h] [-b] [-g input obo file] [-o output detail file] [-s output summary file] [-r output error file] [-e evd rank file]  -i input gene association file

Required parameters:
        -i input gene association file.

Optional parameters:
    -h displays this message 
    -b do not exclude obsolete go ids.
                Default is to exclude obsolete go ids.
    -g alternative name for the gene ontology obo file.
                Default name is ontology/go-basic.obo
    -e alternative name for the evidence rank file.
                Default name is ontology/evidence_rank.txt
    -o alternative name for the output detail file
                Default name is GA2GAQ_detail.txt
    -s alternative name for the output summary file
                Default name is GA2GAQ_summary.txt
    -r alternative name for the output error file
                Default name is GA2GAQ_Error.txt
                This file will only be created if errors are found.
Examples:
    Generate file of GAQ scores 
        % perl $0 -i gene_association.txt 

    Generate file specifying output names
        % perl $0 -i gene_association.txt -o gaq_detail.txt -s gaq_summary.txt

    Print Help message
        % perl $0 -h
END
    exit;
  }

  if (!$opt_i) {
    die "Input filename must be provided.\nExiting now.\n\n";
  }
  else {
    if (!-f $opt_i) {
      die "Input filename $opt_i must exist.\nExiting now.\n\n";
    }
    if (-z $opt_i) {
      die "Input filename $opt_i contains no records.\nExiting now.\n\n";
    }
    if (!-r $opt_i) {
     die "Input filename $opt_i cannot be read.\nExiting now.\n\n";
    }
    $input_file = $opt_i;
   } 

   if ($opt_g)     {
        $obofile = $opt_g;
  }
  if (!-f $obofile) {
      die "Ontology filename must exist.\nExiting now.\n\n";
   }
  elsif (-z $obofile) {
      die "Ontology filename contains no records.\nExiting now.\n\n";
  }
  elsif (!-r $obofile) {
     die "Ontology filename cannot be read.\nExiting now.\n\n";
  }
  my $st = stat($obofile) or die "Cannot stat $obofile: $!";
  my $mtime = $st->mtime;
  my $age = ($currDate - $mtime) / 86400; # convert seconds to days
  if ($age > 30)
  {
    print STDERR "Ontology file is more than 30 days old.  Might want to consider ftping new file\n";
  }

  if (!-f $depthfile) {
      die "GO Term depth file must exist.\nExiting now.\n\n";
   }
  elsif (-z $depthfile) {
      die "GO Term depth file  contains no records.\nExiting now.\n\n";
  }
  elsif (!-r $depthfile) {
     die "GO Term depth filen cannot be read.\nExiting now.\n\n";
  }
  $st = stat($depthfile) or die "Cannot stat $depthfile: $!";
  $mtime = $st->mtime;
  $age = ($currDate - $mtime) / 86400; # convert seconds to days
  if ($age > 30)
  {
    print STDERR "GO Depth file is more than 30 days old.  Might want to consider downloading new file from AgBase. \n";
  }

   if ($opt_e)     {
        $evdfile = $opt_e;
  }
  if (!-f $evdfile) {
      die "Evidence rank filename must exist.\nExiting now.\n\n";
   }
  elsif (-z $evdfile) {
      die "Evidence rank filename contains no records.\nExiting now.\n\n";
  }
  elsif (!-r $evdfile) {
     die "Evidence rank filename cannot be read.\nExiting now.\n\n";
  }

  if ($opt_o)     {
     $out_detail = $opt_o;
    }
 
  if ($opt_s)     {
     $out_summary = $opt_s;
    }

  if ($opt_r)     {
     $out_error = $opt_r;
    }
 
  if ($opt_b) {
      $ex_obs = 0;
   }
  else {
     $ex_obs = 1;
   } 
} #end of subroutine

################################################################
#
# sub_buildGO
#
# Desc: this subroutine will build a GO ID -Name/Obsolete  hash for
#   each Term id in the gene_ontology file
################################################################
sub sub_buildGO {
   my ($id, $alt_id, $term, $defCnt, $obsolete, $i);
   $id = $alt_id = $term = '';
   $obsolete = $defCnt = $i = 0;
   open(INFILE, "$obofile") || die "$!\n";
   while(my $line = <INFILE>){
        chomp($line);
        if($line =~ m/^\[Term\]/){
           $defCnt++;
           if ($i != 0) {
               if ($id ne '' && $term ne '') {
                 $GOName{$id} = $term; 
                 if ((defined $alt_id)  && ($alt_id ne '')) {
                     my @alts = split(/\|/,$alt_id);
                     foreach my $alt (@alts) {
                        $GOName{$alt} = $term;
                     }
                 } # defined alternate ides
                 if ((defined $obsolete)  && ($obsolete == 1)) {
                    $GOObs{$id} = 1;
                    if ((defined $alt_id)  && ($alt_id ne '')) {
                       my @alts = split(/\|/,$alt_id);
                        foreach my $alt (@alts) {
                          $GOObs{$alt} =1;
                        }
                   }
                 } # defined obsolete
               } # have an id and term
           } # id not first record
           $i++;
           $id = $alt_id =  $term =  '';
           $obsolete= 0;
        } # Term
        if (($line =~ m/^\[Typedef\]/) ||
           ($line =~ m/^\[Instance\]/)){
           if ($i != 0) {
               if ($id ne '' && $term ne '') {
                 $GOName{$id} = $term; 
                 if ((defined $alt_id)  && ($alt_id ne '')) {
                     my @alts = split(/\|/,$alt_id);
                     foreach my $alt (@alts) {
                        $GOName{$alt} = $term;
                     }
                 } # defined alternate ides
                 if ((defined $obsolete)  && ($obsolete == 1)) {
                    $GOObs{$id} = 1;
                    if ((defined $alt_id)  && ($alt_id ne '')) {
                       my @alts = split(/\|/,$alt_id);
                        foreach my $alt (@alts) {
                          $GOObs{$alt} =1;
                        }
                   }
                 } # defined obsolete
               } # have an id and term
           } # not the first record
           $id = $alt_id =  $term =  '';
           $obsolete= 0;
        }
        if ($line =~ m/^id:\s+(\S+)/)   { $id = $1; }
        if($line =~ m/^def:\s+(.*)/){
                $term= $1;
                $term=~ s/"//g;
        }
        if($line =~ m/^alt_id:\s+(.*)/){
          if ($alt_id) { $alt_id .= '|'.$1; }
          else { $alt_id = $1; }
        }
        if($line =~ m/^is_obsolete:\s+(.*)/){ $obsolete = 1; }
   } # while input data
close(INFILE);

# set the last GO ID if not already processed
if (($id ne '') && ($term ne '')) {
   $GOName{$id} = $term; 
   if ((defined $alt_id)  && ($alt_id ne '')) {
       my @alts = split(/\|/,$alt_id);
       foreach my $alt (@alts) {
          $GOName{$alt} = $term;
       }
   } # defined alternate ides
   if ((defined $obsolete)  && ($obsolete == 1)) {
       $GOObs{$id} = 1;
       if ((defined $alt_id)  && ($alt_id ne '')) {
           my @alts = split(/\|/,$alt_id);
           foreach my $alt (@alts) {
               $GOObs{$alt} =1;
           }
        }
    } # defined obsolete
} # have an id and term

} #end of subroutine


################################################################
#
# sub_getDepth
#
# Desc: this subroutine will build the depth of each GO term
#
################################################################
sub sub_getDepth {
  my ($line, @flds, $go_id, $go_depth);
  my $depth_cnt=0;
  open (IN,"<$depthfile");
  while (<IN>) {
      $line = $_;
      chomp $line;
      @flds = split(/\t/,$line);
      $go_id = $flds[0];
      $go_depth= $flds[1];
      if ((defined $go_id) && ($go_id ne ' ' ) &&
          ($go_id ne 'GO_ID') &&
          (defined $go_depth) && ($go_depth ne ' ') &&
          ($go_depth !~ m/ath/)
          ) {
              $GODepths{$go_id} = $go_depth;
          }
  }
  close IN;
  $depth_cnt += keys %GODepths;
  if ($depth_cnt == 0 ) {
     print STDERR "No GO Depth records found.\n";
     exit;
  }
} #end of subroutine

################################################################
#
# sub_getRank
#
# Desc: this subroutine will build the evidence hash
#
################################################################
sub sub_getRank{
  my ($line, @flds, $ecode, $rank, $rank_cnt);
  $rank_cnt=0;
  open (IN,"<$evdfile");
  while (<IN>) {
      $line = $_;
       chomp $line;
      @flds = split(/\t/,$line);
      $ecode = $flds[0];
      $rank = $flds[1];
      if ((defined $ecode) && ($ecode ne ' ' ) &&
          (defined $rank) && ($rank ne ' ') ) {
          if ($rank =~ /^[+-]?\d{1,3}$/)   {
             $Ranks{$ecode} = $rank;
          }
          else {
            print STDERR "Rank either non-numeric or too large for $ecode. Must be numeric and between 1 and 3 digits.\n";
         }
      }
  }
  close IN;
  $rank_cnt += keys %Ranks;  
  if ($rank_cnt == 0 ) {
     print STDERR "No evidence rank records found.\n";
     exit;
  }      
} #end of subroutine

################################################################
sub sub_parseGA($){
###############################################################
    my $ga_file = shift;
    my $ga = [];

    local $/ = "\n";

   my %seen;
   open(IN, "<$ga_file");
   while (my $line = <IN>){
      chomp($line);

      my $result = ();
     $result->{error} = [];
      my ($pid, $goid, $evd, $goname, $dbRef, $withFrom, $aspect, $date);
      my $line2 = '';
      my $continue=1;

        # ignore comment lines
     if ($line =~ /^!.*/) { $continue = 0; }
     if ($line =~ /^Database/) { $continue = 0; } 
     if ($line =~ /^Input_Accession/) { $continue = 0; } 
     if ($line =~ /^Input Accession/) { $continue = 0; } 
  
     my $size=1;
     my @tmp;
    if ($continue) {
      @tmp = split("\t", $line,-1);
      $size = $#tmp+1;
      if ($size < 15) {
         push(@{$result->{error}}, { type=> 'Not enough columns', line=>$line });
         $continue = 0;
      }
    } # continue
    #
    # try to accommodate both a 15 column ga file and the output from GoRetreiver
    #   which adds an extra column GOTermName in column 6
    #
    if ($continue) {
       if ($size== 19) {
           $pid = $tmp[2];
           $goid = $tmp[5];
           if (defined $goid && $goid ne '') {$goid =~ s/^\s+//; }
           if (defined $goid && $goid ne '') { $goid=~ s/\s+$//; }
           $evd = $tmp[8];
           if (defined $evd && $evd ne '') {$evd =~ s/^\s+//; }
           if (defined $evd && $evd ne '') { $evd=~ s/\s+$//; }
           $goname = '';
           $dbRef =  $tmp[7];
           $withFrom = $tmp[9];
           $aspect = $tmp[10];
           $date = $tmp[15];
        }
        elsif ($size== 17) {
            $pid = $tmp[1];
            $goid = $tmp[4];
            $evd = $tmp[6];
            $goname = '';
            $dbRef =  $tmp[5];
            $withFrom = $tmp[7];
            $aspect = $tmp[8];
            $date = $tmp[13];
      }
      else {
           $pid = $tmp[1];
           $goid = $tmp[4];
           $evd = $tmp[7];
           $goname = $tmp[5];
           $dbRef =  $tmp[6];
           $withFrom = $tmp[8];
           $aspect = $tmp[9];
           $date = $tmp[14];
      }
        $pid =~ s/\r|\n|\t|\f//g;
        $goid =~ s/^\s+|\s+$//;
        $goid =~ s/\r|\n|\t|\f//g;
        $evd =~ s/^\s+|\s+$//;
        $evd =~ s/\r|\n|\t|\f//g;
        $goname =~ s/\r|\n|\t|\f//g;
        $dbRef =~ s/\r|\n|\t|\f//g;
        $withFrom =~ s/\r|\n|\t|\f//g;
        $aspect =~ s/\r|\n|\t|\f//g;
        $date =~ s/\r|\n|\t|\f//g;

        if ($goid !~ m/^GO:\d{7}/) {
          push(@{$result->{error}}, { type=> 'Invalid GO ID', line=>$line });
          $continue = 0;
        }
    }  # continue -- line had at least 15 columns

    $result->{gaq} = [];
    if ($continue) {
       $line2 = "$pid$goid$goname$dbRef$withFrom$aspect$evd$date";
       if (grep { ! $seen{$_} ++ } $line2) {
         push(@{$result->{gaq}}, { pid => $pid,
                  goid => $goid,
                  goname => $goname,
                  dbRef => $dbRef,
                  withFrom => $withFrom,
                  aspect => $aspect,
                  evd => $evd,
                 date => $date });
        } #unique lines only
    } # continue
    push(@$ga,$result);
}  # while input lines
return $ga;
close IN;   
}           
################################################################
sub sub_writeGAQ($$){
###############################################################
            
my $results = shift;
my $ex_obs= shift;
       
my $errorsCnt=0; 
my %total = ();

open(GA2GAQ_SUM,"> $out_summary") or die "Can not open output summary text file\n";
open(GA2GAQ_DET,"> $out_detail") or die "Can not open output detail text file\n";
open(GA2GAQ_ERRORS,"> $out_error") or die "Can not open output error file\n";



print GA2GAQ_DET join ("\t", ('Gene_product_ID','GO_ID','GO_name',
    'DB:Reference','With(or)From','Aspect','Evidence', 'GAQ_score', 'Annotation_date')), "\n";


foreach my $result (@{$results->{gaq_info}}){
  foreach my $errors (@{$result->{error}}){
    chomp($errors->{line});
    $errors->{line} =~ s/\r\n//g;
    $errors->{line} =~ s/\n//g;
    print GA2GAQ_ERRORS "$errors->{line}  Error:$errors->{type}\n";
    $errorsCnt++;
  }
  foreach my $gaq(@{$result->{gaq}}){
        my $line2= join ("\t", (
           $gaq->{pid},
           $gaq->{goid},
           $gaq->{goname},
           $gaq->{dbRef},
           $gaq->{withFrom},
           $gaq->{aspect},
           $gaq->{evd},
           $gaq->{date}));
       #print "$line2<p>";

        my $go_id = $gaq->{goid};
        my $evd_code = $gaq->{evd};
        my $gname = $gaq->{goname};
        if ((!defined $gname) || ($gname eq ' ' ) || ($gname eq '' )) {
            if (exists $GOName{$go_id}) {
               $gname = $GOName{$go_id};
            }
        }

      my $continue = 1;

      my $rank = $Ranks{$evd_code};

      my $depth=-1;
      if (exists $GODepths{$go_id}) { 
         $depth = $GODepths{$go_id};
      }

        if ($GOObs{$go_id} ) {
           if ($ex_obs) {
                print GA2GAQ_ERRORS "$line2 GO ID is obsolete record will not be in output.\n";
                $continue = 0;
            }
          else {
                print GA2GAQ_ERRORS "$line2 GO ID is obsolete. GO Depth will be equal to 1.\n";
               $depth=1;
          }  
        }
  
       if ($continue) {
         if ((!defined $depth) || ($depth < 0)) {
           print GA2GAQ_ERRORS "$line2 GO Depth not found.\n";
           $continue = 0;
           $depth=0;
        }
       }

       if ($continue) {
         if ((!defined $rank) || ($rank < 0)) {
             print GA2GAQ_ERRORS "$line2 Evidence Code Rank not found\n";
             $continue = 0;
             $rank= 0;
         }
       }

       my $gaq_score = $depth*$rank;
       if ($continue) {
          print GA2GAQ_DET join ("\t", (
            $gaq->{pid},
            $gaq->{goid},
            $gname,
            $gaq->{dbRef},
            $gaq->{withFrom},
            $gaq->{aspect},
            $gaq->{evd},
            $gaq_score,
            $gaq->{date})), "\n";
          my $keys = $gaq->{pid};
          my $values = $gaq_score;
          $total{$keys} += $values;
        }

   } #foreach gaq record
  } #foreach result output record
 close GA2GAQ_DET;
 close GA2GAQ_ERRORS;

  print GA2GAQ_SUM "Gene_product_ID\tGAQ_scores\n";

# Summarize the GAQ scores
my $totalGAQ = 0;
my $count = 0;
foreach my $keys (keys %total) {
        print GA2GAQ_SUM "$keys\t $total{$keys}\n";
        $count++;
        $totalGAQ = $totalGAQ + $total{$keys};
}
my $mean = 0;
if ($count != 0) {
$mean = $totalGAQ/$count;
}

print GA2GAQ_SUM "\nSUMMARY\n";
print GA2GAQ_SUM "Total GAQ score:\t$totalGAQ\n";
print GA2GAQ_SUM "Number of non-redundant gene products:\t$count\n";
print GA2GAQ_SUM "Mean GAQ score:\t$mean\n";

close GA2GAQ_SUM;


} #end of writeGAQ subroutine


