
-------------------------GA2GAQ--------------------------------
GA2GAQ.pl is a standalone script 
that allows users to calculate a GAQ score for a gene 
association file. 

GAQ is a tool that allows users to generate scores to 
be used to quantitatively measure the quality of GO annotation of
a set of gene products.

GAQ scores include the breadth of GO annotation, the level of detail
of annotation and the type of evidence used to make the annotation. 

The scores generated can also be used by annotators to track changes 
in GO annotations over time.

The depth of each GO Term is defined in a text file that is updated weekly(at the same time the obo file is updated) in the downloadable zip file.  The GO term depth is generated using Gavin Sherlock's GO::TermFinder / GO::Node.  

The tool will determine the depth of each GO term and the rank of 
each evidence code associated with the annotation and returns a GAQ score 
as a product of depth and evidence code rank. 

The total GAQ score of each annotated gene product will also be calculated
and a summary will be generated showing the overall total GAQ scores, the 
number of gene products annotated and the average (mean) GAQ score of the whole set. 

More information about the tool can be found at 
http://nar.oxfordjournals.org/cgi/content/full/36/2/e12


The input file need to be a tab-delimited
text file from the AgBase GORetriever tool or a gene association 
formatted file.  

This tool was written in Perl therefore it should run
on systems that support Perl.

There is also a web based version of this tool available 
at http://agbase.msstate.edu/cgi-bin/GA2GAQ.cgi


I   -  Requirements and setup for GA2GAQ

II  -  How the tool works/What the tool does

III -  Computational requirements for GA2GAQ

IV  -  Sample Input Records

V   -  Sample Output Records

VI  -  Software use disclaimer

I - Requirements and setup for GA2GAQ
-------------------------------------

I - 1) Requirements 

     The GA2GAQ tool is a Perl script.

     Therefore the system must have Perl installed.
 
        a) Perl
           Specific Perl packages include (these should be standard packages):
 		File::stat;
                Time::localtime;
                File::Temp;
 		Getopt::Std;

        b) go-basic.obo 

           This text file will be required to run the script.

           The default location should be in a subdirectory called ontology/
           that exists in the directory where the perl script is placed.

           You will be able to override this default location if you have this file
           in another location by specifying the -g option at execution time.  OR
           you can always edit the script and change the file location.
        
           HOWEVER since the goDepth.txt file is generated in the zip file at the 
           same time as the obo file is downloaded I would not point to a new location for
           an obo file.

           To update the obo file or the GOdepth.txt file download the zip file again from AgBase.
          
           The files should be updated weekly.

           This text file is used to build a list of all the valid GO ids,
             whether the id is obsolete, the name of the GO Term (if the name
              is not provided in the input file(that is the format is not the 
             AgBase GORetriever output) and will calculate the depth of each term.

           If this file is older than 30 days then the script will remind you that
           you might want to download a new version of this file.

           You can download the GA2GAQ standalone zip file from    
 		http://agbase.msstate.edu/cgi-bin/GA2GAQ.cgi
    
       c) goDepth.txt
    
           This text file is required to run the script.                

           The default location should be in a subdirectory called ontology/
           that exists in the directory where the perl script is placed.

           This text file is build using the obo file from step b) and Ganin Sherlock's
           GO::TermFinder / GO::Node perl modules.  
   
           The GO Depths from this text file are the longest paths for a go term.

           If this file is older than 30 days then the script will remind you that
           you might want to download a new version of this file. 
    
           This file is updated weekly in the AgBase GA2GAQ standalone zip file.

           You can download the GA2GAQ standalone zip  file from    
 		http://agbase.msstate.edu/cgi-bin/GA2GAQ.cgi

        c) evidence_rank.txt

           This text file will be required to run the script.

           The default location should be in a subdirectory called ontology/
           that exists in the directory where the perl script is placed.

           You will be able to override this default location if you have this file
           in another location by specifying the -e option at execution time.  OR
           you can always edit the script and change the file location.

           This text file is used to calculate the GAQ score. 
           
           The file must contain the Evidence Code and a numeric rank (1-3 digits) 
           assigned to each code.  The 2 columns must be separated by a tab character.

           You can download this file* from    
 		http://agbase.msstate.edu/cgi-bin/GA2GAQ.cgi

          *Note: These ranks were determined to be the best at determining GAQ scores as 
           explained in the GAQ manuscript (supplementary data 1) which can be found at 
           http://nar.oxfordjournals.org/cgi/content/full/36/2/e12

        d) ontology/ subdirectory must be defined in the directory where
           you are running GA2GAQ OR you must define an alternate location 
           for the go-basic.obo file by specifying the -g option
           and you will also be required to list an alternative location for the
           evidence_rank.txt file by using the -e option.

       
I - 2) Setup

        Download the following files:
        ----------------------------
        a) GA2GAQ.zip

           This file can be downloaded from http://agbase.msstate.edu/cgi-bin/GA2GAQ.cgi

           The zip file will contain:
              GA2GAQ.pl
              GA2GAQ_README.TXT
              ontology/go-basic.txt
              ontology/goDepth.txt
              ontology/evidence_rank.txt
       
II  -  How the tool works/What the tool does

---------------------------------------------
Usage:  perl GA2GAQ.pl [-h] [-b] [-g input obo file] [-o output detail file] [-s output summary file] [-r output error file] [-e evd rank file] -i input gene association file

Required parameters:
     -i name to a gene association or GORetriever tab-delimited text file.

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
    -r alternative name for the output error file Default name is GA2GAQ_Error.txt
                This file will only be created if errors are found.
Examples:
    Generate file of GAQ scores
        % perl GA2GAQ.pl -i gene_association.txt

    Generate file specifying output names
        % perl GA2GAQ.pl -i gene_association.txt -o gaq_detail.txt -s gaq_summary.txt

    Generate a file of GAQ scores designating do not exclude obsolete go ids
        % perl GA2GAQ.pl -i gene_association.txt -b

    Print Help message
        % perl GA2GAQ.pl -h


This tool builds a hash of the depth of each go id by reading the gene_ontology obo file.
This tool builds a hash of evidence code ranks by reading the evidence_rank.txt file.       
This tool parses the input gene association text file. 

It will ignore any lines that begin with !  these should be comment lines

Script Validation Checks:
--------------------------
  If any line/column fails a check it should print in the error file.

   1) Each line contains at least 15 columns. 
   2) The GO ID column is in the format of GO:NNNNNNN.  
   3) If the default option of exclude obsolete GO IDs is not overridden
      (by using the -b option) then the program consider these annotations in the 
      GAQ calculation and will print the line in the error file.  

      If you choose not to exclude obsolete GO IDs then the GO ID should appear in 
      the GAQ calculation, but it will also print in the error file so that you know 
      the GO ID is obsolete and a depth of 1 was used for the record.
   4) If the script is unable to provide a GO Term Name then it will print 
         go id not found in gene ontology file.

Output:
-------

There can be up to three output files from the script.

   1) Detail output file (default name GAQ_Detail.txt) 
	The Detail file will contain a line for each annotation considered in the GAQ score.
	If obsolete GO IDs are excluded from the GAQ calculation then they will not be listed 
	in the detail file.

	The format of the output is:
	Gene_product_ID, GO_ID, GO_name, DB:Reference, With(or)From, Aspect, Evidence, GAQ_score, Annotation_date

	A tab is used to delimit the fields. 

   2) Summary output file (default name GAQ_Summary.txt) 
	The Summary file will contain a line for each gene product and the sum of the gene products         GAQ scores.
 
        The Summary file will also contain a Total (Overall) GAQ score, the number of non-redundant
 	gene products and a Mean GAQ score.
              
	The format of the output is:
	Gene_product_ID, GAQ_Scores

        Then Summary and the 3 overall total lines

	A tab is used to delimit the fields. 

   3) Error output file (default name is GAQ_Error.txt)
      The error file will list each line that contained a warning or error.
      Examples of warnings/errors:
           GO ID is obsolete will not be in output (if exclude option set)
           GO Depth not found 
           GO ID is obsolete will be in output but depth will be set to 1.
     
III -  Computational requirements for GOanna
---------------------------------------------
This tool should not be computationally intensive.  


IV  -  Sample Input Records
---------------------------------------------
The input records must be tab-delimited text records.

Example of gene association formatted file*:
-------------------------------------------------------------------------------------
UniProtKB	Q84J55	GF14A		GO:0019904	GOA:interpro|GO_REF:0000002	IEA	InterPro:IPR000308	F	14-3-3-like protein GF14-A	GF14A|Os08g0480800|OJ1113_A10.40|LOC_Os08g37490|OSJNBb0092C08.10|14331_ORYSJ	protein	39947	20090112	UniProtKB
UniProtKB	Q84J55	GF14A		GO:0000067	GOA:interpro|GO_REF:0000002	IEA	InterPro:IPR000308	F	14-3-3-like protein GF14-A	GF14A|Os08g0480800|OJ1113_A10.40|LOC_Os08g37490|OSJNBb0092C08.10|14331_ORYSJ	protein	39947	20090112	UniProtKB
UniProtKB	Q84J55	GF14A		GO:0000174	GOA:interpro|GO_REF:0000002	IEA	InterPro:IPR000308	F	14-3-3-like protein GF14-A	GF14A|Os08g0480800|OJ1113_A10.40|LOC_Os08g37490|OSJNBb0092C08.10|14331_ORYSJ	protein	39947	20090112	UniProtKB

Example of GORetriever formatted file* (same format as gene association file but with GO Term Name added as 6th column):
-------------------------------------------------------------------------------------------------------------------------
UniProtKB	Q84J55	GF14A		GO:0019904	protein domain specific binding	GOA:interpro|GO_REF:0000002	IEA	InterPro:IPR000308	F	14-3-3-like protein GF14-A	GF14A|Os08g0480800|OJ1113_A10.40|LOC_Os08g37490|OSJNBb0092C08.10|14331_ORYSJ	protein	39947	20090112	UniProtKB
UniProtKB	Q84J55	GF14A		GO:0000067	protein domain specific binding	GOA:interpro|GO_REF:0000002	IEA	InterPro:IPR000308	F	14-3-3-like protein GF14-A	GF14A|Os08g0480800|OJ1113_A10.40|LOC_Os08g37490|OSJNBb0092C08.10|14331_ORYSJ	protein	39947	20090112	UniProtKB
UniProtKB	Q84J55	GF14A		GO:0000174	protein domain specific binding	GOA:interpro|GO_REF:0000002	IEA	InterPro:IPR000308	F	14-3-3-like protein GF14-A	GF14A|Os08g0480800|OJ1113_A10.40|LOC_Os08g37490|OSJNBb0092C08.10|14331_ORYSJ	protein	39947	20090112	UniProtKB

	

V   -  Sample Output Records
--------------------------------------------

GAQ_Detail.txt: 
Gene_product_ID GO_ID   GO_name DB:Reference    With(or)From    Aspect  Evidence        GAQ_score       Annotation_date
Q84J55  GO:0019904      protein domain specific binding GOA:interpro|GO_REF:0000002     InterPro:IPR000308      F       IEA     8       20090112
Q06967  GO:0005634      nucleus GOA:spkw|GO_REF:0000004 SP_KW:KW-0539   C       IEA     10      20090112
Q06967  GO:0005634      nucleus GOA:spsl|GO_REF:0000023 SP_SL:SL-0191   C       IEA     10      20090112
Q06967  GO:0005737      cytoplasm       GOA:spkw|GO_REF:0000004 SP_KW:KW-0963   C       IEA     12      20090112
Q06967  GO:0005737      cytoplasm       GOA:spsl|GO_REF:0000023 SP_SL:SL-0086   C       IEA     12      20090112
Q06967  GO:0019904      protein domain specific binding GOA:interpro|GO_REF:0000002     InterPro:IPR000308      F       IEA     8       20090112
Q6Q6X0  GO:0005737      cytoplasm       GOA:spkw|GO_REF:0000004 SP_KW:KW-0963   C       IEA     12      20090112
Q6Q6X0  GO:0005737      cytoplasm       GOA:spsl|GO_REF:0000023 SP_SL:SL-0086   C       IEA     12      20090112
Q6Q6X0  GO:0019904      protein domain specific binding GOA:interpro|GO_REF:0000002     InterPro:IPR000308      F       IEA     8       20090112


GAQ_Summary.txt:
Gene_product_ID GAQ_scores
Q84J55 8
Q06967 52
Q6Q6X0 32

SUMMARY
Total GAQ score: 92
Number of non-redundant gene products: 3
Mean GAQ score: 30.6666666666667


GAQ_Error.txt:
Q84J55  GO:0000067      protein domain specific binding GOA:interpro|GO_REF:0000002     InterPro:IPR000308      F       IEA     20090112 GO ID is obsolete record will not be in output.
Q84J55  GO:0000174      protein domain specific binding GOA:interpro|GO_REF:0000002     InterPro:IPR000308      F       IEA     20090112 GO ID is obsolete record will not be in output.
Q06967  GO:0004863      diacylglycerol-activated phospholipid-dependent protein kinase C inhibitor activity     PMID:8405471            F       IDA     20060831 GO Depth not found.

VI   -  Software use disclaimer
--------------------------------------------

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

