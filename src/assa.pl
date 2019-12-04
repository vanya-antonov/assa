#!/usr/bin/perl --

use strict;
use warnings;

my $SVN_STR  = '$Id$';
my $ASSA_VER = '1.00';

###
# Ivan Antonov (antonov1986@gmail.com)
#

$| = 1; # Turn off buffering

use FindBin;
use lib $FindBin::RealBin;

use Data::Dumper;
use Getopt::Long;
use File::Spec;
use File::Path qw(remove_tree);
use Cwd qw(abs_path);
use IPC::Open2;

use ASSA::ASSA;
use ASSA::Lib qw(compute_gc_content compute_ali_identity revcomp sum);

###
# CONSTANTS
my $Q_MASKED_FN  = undef;
my $NUM_THREADS  = 1;
my $MIN_SCORE    = 36;
my $LASTAL_A     = 12;   # -a: gap existence cost
my $LASTAL_B     = 6;    # -b: gap extension cost
my $MIN_SITE_LEN = 0;
my $MAX_SITE_LEN = undef;
my $MAX_W        = undef;
my $ALL_SITES    = undef;
my $TMP_DIR      = '.';
my $SILENT       = 0;

###
# Parse input data
GetOptions(
	'q_masked=s'      => \$Q_MASKED_FN,
	'num_threads=i'   => \$NUM_THREADS,
	'min_score=i'     => \$MIN_SCORE,
	'lastal_a=i'      => \$LASTAL_A,
	'lastal_b=i'      => \$LASTAL_B,
	'min_site_len=i'  => \$MIN_SITE_LEN,
	'max_site_len=i'  => \$MAX_SITE_LEN,
	'max_w=i'         => \$MAX_W,
	'all_sites=s'     => \$ALL_SITES,
	'tmp_dir=s'       => \$TMP_DIR,
	'silent'          => \$SILENT,
) || die usage();

die usage() if @ARGV!=2;

###
my $START_TIME = time;

run(
	q_fn            => $ARGV[0],
	h_fn            => $ARGV[1],
	q_masked_fn     => $Q_MASKED_FN,
	num_threads     => $NUM_THREADS,
	lastal_score    => $MIN_SCORE,
	lastal_a        => $LASTAL_A,
	lastal_b        => $LASTAL_B,
	min_site_len    => $MIN_SITE_LEN,
	max_site_len    => $MAX_SITE_LEN,
	rnaup_max_w     => $MAX_W,
	all_sites       => $ALL_SITES,
	tmp_root        => $TMP_DIR,
	verbose         => $SILENT ? 0 : 1,
);

warn "\nElapsed time: ".(time-$START_TIME)." sec\n" if !$SILENT;
###

###
# SUBROUTINES
sub run
{
	my %opts = @_;
	
	$opts{tmp_dir} = create_tmp_dir(root_dir => $opts{tmp_root}, prefix => '__assa_');
	
	# Read all the sequences
	my $assa = ASSA::ASSA->new(%opts);
	
	# Run LASTAL and identify putative sites
	$assa->run_lastal(min_site_len => $opts{min_site_len}, max_site_len => $opts{max_site_len});
	
	# Based on the LASTAL sites prepare regions and run RNAup
	$assa->run_RNAup();
	
	$assa->compute_pvalues();
	
	output_site_table  ($assa, $opts{all_sites}) if  $opts{all_sites};
	output_interactions($assa, %opts);
	
	remove_tree($opts{tmp_dir}) or warn "Couldn't remove tmp dir $opts{tmp_dir}";
}

sub output_interactions
{
	my($assa, %opts) = @_;
	
	my @head  = qw(
		name1   name2   len1   len2   gc1   gc2   numSites   sumEnergy   Pvalue   Padj
	);
	print "#\t".join("\t", @head)."\n";
	
	my $i = 1;
	foreach my $trx_pair ( sort { $a->{PVALUE_ADJ} <=> $b->{PVALUE_ADJ} } @{$assa->get_trx_pairs_aoh} )
	{
		my $q_trx = $assa->{q_trx_hash}{$trx_pair->{Q_TRX_ID}};
		my $h_trx = $assa->{h_trx_hash}{$trx_pair->{H_TRX_ID}};
		my $h = {
			name1     => $trx_pair->{Q_TRX_ID},
			name2     => $trx_pair->{H_TRX_ID},
			len1      => length($q_trx->{SEQ}),
			len2      => length($h_trx->{SEQ}),
			gc1       => 100*compute_gc_content($q_trx->{SEQ}),
			gc2       => 100*compute_gc_content($h_trx->{SEQ}),
			numSites  => scalar(keys %{$trx_pair->{ALL_SITES}}),
			sumEnergy => $trx_pair->{SUM_ENERGY},
			Pvalue    => sprintf('%.2e', $trx_pair->{PVALUE}),
			Padj      => sprintf('%.2e', $trx_pair->{PVALUE_ADJ}),
		};
		print join("\t", $i++, map { defined $h->{$_} ? $h->{$_} : '-' } @head)."\n";
	}
}

sub output_site_table
{
	my($assa, $out_fn) = @_;
	
	my @head  = qw(
		name1        name2
		siteStart1   siteEnd1   siteStart2   siteEnd2   siteLen   siteCmpl   siteScore
	    regStart1    regLen1    regStart2    regLen2    regGC     deltaG
	);
	
	my $reg_hash = $assa->get_region_hash();
	my $site_aoh = [sort {$a->{DELTA_G} <=> $b->{DELTA_G}} @{$assa->get_all_sites_aoh}];
	
	open(OUT, '>', $out_fn) or die "Can't open file to write: $!";
	print OUT "#\t".join("\t", @head)."\n";
	for(my $i = 0; $i < scalar(@$site_aoh); $i++)
	{
		my $site   = $site_aoh->[$i];
		my $region = $reg_hash->{$site->{SITE_ID}};
		
		my %h = (
			name1      => $region->{Q_TRX_ID},
			name2      => $region->{H_TRX_ID},
			siteStart1 => $site->{Q_START},
			siteEnd1   => $site->{Q_END},
			siteStart2 => $site->{H_START},
			siteEnd2   => $site->{H_END},
			siteLen    => $site->{SITE_LEN},
			siteCmpl   => $site->{SITE_CMPL},
			siteScore  => $site->{SCORE},
			regStart1  => $region->{Q_START},
			regLen1    => length($region->{Q_SEQ}),
			regStart2  => $region->{H_START},
			regLen2    => length($region->{H_SEQ}),
			regGC      => compute_gc_content($region->{Q_SEQ}.$region->{H_SEQ}),
			deltaG     => $site->{DELTA_G},
		);
		
		print OUT join("\t", $i+1, map { defined $h{$_} ? $h{$_} : '-' } @head)."\n";
		
		my $rnaup_seqs = $site->{RNAUP}{ALL_LINES}[1];
		$rnaup_seqs =~ tr/Uu/Tt/;
		print OUT "# $site->{RNAUP}{DOT_BRACKET}\n# $rnaup_seqs\n";
	}
	close OUT;
}

sub create_tmp_dir
{
	my(%opts) = @_;
	$opts{root_dir} ||= '.';
	$opts{prefix}   ||= '';
	
	my $dir;
	while(1)
	{
		$dir = $opts{root_dir}.'/'.$opts{prefix}.int(rand(1000000));
		last unless -d $dir;
	}
	
	mkdir($dir, 0777);
	
	return abs_path($dir);
}

sub usage
{
	my($msg) = @_;
	$msg = $msg ? $msg."\n" : '';
	
	# $Id$
	my($revision, $date) = $SVN_STR =~ /assa.*?\.pl\s+(\S+)\s+(\S+)/;
	
	my $script = File::Spec->splitpath($0);
	
	return"$msg
ASSA, version $ASSA_VER, revision $revision ($date)

USAGE:
    $script   [OPTIONS]   <QUERIES.fasta>   <TARGETS.fasta>   >   <INTERACTIONS.txt>

OPTIONS:
    --num_threads <INT>  --  number of threads to use ($NUM_THREADS)
    --tmp_dir     <DIR>  --  where the temporary folder will be created ($TMP_DIR)
    --q_masked <FN.fna>  --  masked versions of the sequences from <QUERIES.fna> for LASTAL search
    --min_score   <NUM>  --  [LASTAL] minimum score for gapped alignments ($MIN_SCORE)
    --lastal_a    <NUM>  --  [LASTAL] -a: gap existence cost ($LASTAL_A)
    --lastal_b    <NUM>  --  [LASTAL] -b: gap extension cost ($LASTAL_B)
    --min_site_len  <N>  --  [LASTAL] minimal length of the LASTAL alignments ($MIN_SITE_LEN)
    --max_site_len  <N>  --  [LASTAL] maximal length of the LASTAL alignments (infinity)
    --max_w       <INT>  --  [RNAup] maximum value for the -w option of the RNAup (see ASSA_MODEL)
    --all_sites  <FILE>  --  save information about all the regions analyzed by RNAup
    --silent

";
}

