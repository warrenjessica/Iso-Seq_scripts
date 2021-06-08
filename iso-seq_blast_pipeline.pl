#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Bio::SearchIO; 
use sloan;

my $usage = 
"\nUsage: perl $0 [options/arguments]
   
   This script extracts sequences from multiple Iso-Seq clustered sequence
   output files.
   
   REQUIRED ARGUMENTS
   
   Query Fasta File
         --query
         Fasta file containing query sequences to search against Iso-Seq
         databases.
   
   Database Fasata Files
         --db_files
         Comma-separated list of fasta file names containing clustered Iso-Seq
         output. Matching blast databases must also exist.
         
   Output Name
         --output
         Base name for all output files (additional extensions will be added)
   
   
   OPTIONAL ARGUMENTS
 
   Minimum BLAST ID
         --min_id [default: 0]      
         Minimum proportion of identical alignment positions in blast hit
         (expressed as decimal value between 0 and 1).

   Minimum BLAST Length Coverage
         --min_len [default: 0]      
         Minimum proportion of query length coverage in blast hit (expressed
         as decimal value between 0 and 1).

   Minimum Read Number
         --min_reads [default: 0]      
         Minimum number of CCSs that were merged into the Iso-Seq cluster in
         the blast hit, as parsed from coverage term in fasta header line.

   Number of Threads
         --num_threads [default: 1]      
         Number of threads to be used in BLAST search.
                 
";

our $QUERY;
our $DB_FILES;
our $OUTPUT;
our $MIN_ID = 0;
our $MIN_LEN = 0;
our $MIN_READS = 0;
our $NUM_THREADS = 1;

##Print out start time and command line call of Perl script
print "\n" . (localtime) . "\nRunning $0 with the following command:\n", qx/ps -o args $$/, "\n";


GetOptions(
    'query=s'  => \$QUERY,
    'db_files=s'  => \$DB_FILES,
    'output=s'  => \$OUTPUT,
    'min_id=f'  => \$MIN_ID,    
    'min_len=f'  => \$MIN_LEN,    
    'min_reads=i'  => \$MIN_READS,
    'num_threads=i'  => \$NUM_THREADS
);


$DB_FILES or die ("\n$usage\n\nERROR: Comma-delimited list of database files must be provided with --db_files\.\n\n");
$QUERY or die ("\n$usage\n\nERROR: Query file must be provided with --query\.\n\n");
$OUTPUT or die ("\n$usage\n\nERROR: Output basename must be provided with --output\.\n\n");


my @dbs = split (/\,/, $DB_FILES);
my @db_shortnames;
my %fasta_string_hash;

foreach (@dbs){
	
	my @split_path = split (/\//, $_);
	my $file_name = $split_path[-1];
	my @split_file = split (/\./, $file_name);
	my $species = $split_file[0];
	push (@db_shortnames, $species);
	
	system ("tblastn -query $QUERY -db $_ -evalue 0.001 -num_threads $NUM_THREADS -out $OUTPUT\.$species\.blast.txt");
	
	my %fasta = fasta2hash($_);

	my $SearchIO_obj = new Bio::SearchIO(-format => 'blast', -file   => "$OUTPUT\.$species\.blast.txt");

	while( my $result_obj = $SearchIO_obj->next_result ) {
		my $query_name = $result_obj->query_name;
		my $query_len = $result_obj->query_length;
		while ( my $hit_obj = $result_obj->next_hit ) {
			my $hit_name = $hit_obj->name;
			my $hit_desc = $hit_obj->description;
			my $hsp_obj = $hit_obj->next_hsp;
			my $align_len = $hsp_obj->length('query');
			$align_len / $query_len >= $MIN_LEN or next;
			$hsp_obj->frac_identical >= $MIN_ID or next;
			my $coverage;
			if ($hit_desc =~ /coverage\=(\d+)\;/){
				$coverage = $1;
			}else{
				die ("\nERROR: cannot parse coverage from $hit_desc\n\n");
			}
			$coverage >= $MIN_READS or next;
			my @split_transcript = split (/\//, $hit_name);
			my $transcript = $split_transcript[1];
			my $key = $hit_name . " " . $hit_desc;
			$fasta_string_hash{$query_name} .= ">$species\_$transcript\_Cov:$coverage\n$fasta{$key}\n";
		}
	}	
}
	
	
foreach my $seq (sort keys %fasta_string_hash){
	my $FHO = open_output ("$OUTPUT\.$seq\.fas");
	
	print $FHO $fasta_string_hash{$seq};
	
	close $FHO;
	
	system ("mafft $OUTPUT\.$seq\.fas > $OUTPUT\.$seq\.aligned.fas");
	
}	



