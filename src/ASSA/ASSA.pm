package ASSA::ASSA;

use strict;
use warnings;

# $Id$

###
# Ivan Antonov (antonov1986@gmail.com)
#
###
# 
# !!! ALL COORDINATES ARE 1-BASED !!!
#
# $self
#   |
#   |
#   |--->{unique_id}     --  next unique ID
#   |--->{lastdb_path}   --  path to the BLASTn db
#   |--->{tmp_dir}
#   |--->{num_threads}
#   |--->{verbose}
#   |
#   |
#   |--->{q_trx_hash}
#   |--->{h_trx_hash}
#   |          |
#   |          |--->{$trx_id}
#   |                   |
#   |                   |--->{TRX_ID}
#   |                   |--->{FULL_NAME}
#   |                   |--->{SEQ}
#   |                   |--->{SEQ_MASKED}        --   optional
#   |                   |--->{LEN_MASKED_WO_N}   --   optional
#   |
#   |
#   |---{trx_pairs}
#   |       |
#   |       |--->{$q_trx_id}
#   |                 |
#   |                 |--->{$h_trx_id}
#   |                           |
#   |                           |--->{Q_TRX_ID}
#   |                           |--->{H_TRX_ID}
#   |                           |--->{SUM_ENERGY}
#   |                           |--->{PVALUE}
#   |                           |--->{PVALUE_ADJ}
#   |                           |
#   |                           |--->{ALL_SITES}
#   |                                    |
#   |                                    |--->{$site_id}
#   |                                             |
#   |                                             |--->{SITE_ID}          --  unique site ID
#   |                                             |--->{Q_START}
#   |                                             |--->{Q_END}
#   |                                             |--->{Q_LEN}
#   |                                             |--->{H_START}
#   |                                             |--->{H_END}
#   |                                             |--->{H_LEN}
#   |                                             |--->{SITE_LEN}
#   |                                             |--->{SITE_CMPL}
#   |                                             |--->{SCORE}
#   |                                             |--->{DELTA_G}
#   |                                             |--->{RNAUP}            -- hash
#   |
#   |
#   |
#   |
#

use Data::Dumper;
use Carp qw(confess); 
use IPC::Open2;
use POSIX qw(:sys_wait_h);
use Cwd qw(abs_path);

#use Math::CDF;
use Statistics::R;

use ASSA::Lib qw(fasta2hash intersect_digits sum min max revcomp compute_gc_content read_file_to_str vector_product compute_ali_identity unique);

#use Exporter;
#use vars qw(@ISA @EXPORT_OK);
#@ISA       = qw(Exporter);
#@EXPORT_OK = qw(
#	compute_ddG_pvalue
#	compute_dpl_len_pvalue
#);

###
# CONSTANTS


###
# CONSTRUCTOR

sub new
{
	my $class = shift;
	my(%o) = @_;
	
	die "q_fn option is required: ". Dumper(\%o) if !$o{q_fn};
	die "h_fn option is required: ". Dumper(\%o) if !$o{h_fn};
	die "tmp_dir is required: "    . Dumper(\%o) if !$o{tmp_dir} || !-d $o{tmp_dir};
	die "LASTAL score is required: ".Dumper(\%o) if !$o{lastal_score};
	
	# $assa_dir is the ASSA lib dir, e.g. '/Users/antonov/_my/Programming/Perl/lib/ASSA'
	(my $assa_dir = $INC{"ASSA/ASSA.pm"}) =~ s/.ASSA.pm$//;
	
	my $model_fn = $o{model_fn} || "$assa_dir/MODEL.txt";
	
	my $self = bless {
		model        => eval( read_file_to_str($model_fn).";" ),
		q_trx_hash   => fasta2hash($o{q_fn}),
		h_trx_hash   => fasta2hash($o{h_fn}),
		unique_id    => 1,
		num_threads  => $o{num_threads} || 1,
		lastal_score => $o{lastal_score},
		lastal_mtx   => "$assa_dir/LASTAL_RNA_RNA.txt",  # -p: match/mismatch score matrix
		lastal_a     => $o{lastal_a},                              # -a: gap existence cost
		lastal_b     => $o{lastal_b},                              # -b: gap extension cost
		lastdb_path  => "$o{tmp_dir}/target_lastdb",
		tmp_dir      => $o{tmp_dir},
		verbose      => $o{verbose},
	}, $class;
	
	$self->{model}{max_w} = $o{rnaup_max_w} if defined $o{rnaup_max_w};
	
	add_masked_seqs_to_trx_hash($self->{q_trx_hash}, $o{q_masked_fn}) if $o{q_masked_fn};
	
	# verification
	foreach my $trx ( values(%{$self->{q_trx_hash}}) )
	{
		if( $trx->{SEQ} =~ /[^ACGT]/i )
		{
			warn "Sequence '$trx->{TRX_ID}' contains non-ACGT letters and will be ignored!";
			delete $self->{q_trx_hash}{$trx->{TRX_ID}};
		}
	}
	foreach my $trx ( values(%{$self->{h_trx_hash}}) )
	{
		if( $trx->{SEQ} =~ /[^ACGT]/i )
		{
			warn "Sequence '$trx->{TRX_ID}' contains non-ACGT letters and will be ignored!";
			delete $self->{h_trx_hash}{$trx->{TRX_ID}};
		}
	}
	
	return $self;
}

sub add_masked_seqs_to_trx_hash
{
	my($trx_hash, $masked_fn) = @_;
	
	my $masked_hash = fasta2hash( $masked_fn );
	foreach my $trx_id ( keys %$masked_hash )
	{
		warn "Masked file '$masked_fn' contains unknown trx '$trx_id'" and next if !$trx_hash->{$trx_id};
		warn "Masked sequence '$trx_id' has different length!" and next if length($masked_hash->{$trx_id}{SEQ}) != length($trx_hash->{$trx_id}{SEQ});
		$trx_hash->{$trx_id}{SEQ_MASKED} = $masked_hash->{$trx_id}{SEQ};
		
		(my $seq_wo_n = $masked_hash->{$trx_id}{SEQ}) =~ s/[^ACGT]//i;
		$trx_hash->{$trx_id}{LEN_MASKED_WO_N} = length($seq_wo_n);
	}
}


###
# PUBLIC METHODS

sub compute_pvalues
{
	my $self = shift;
	
	my $R = Statistics::R->new();
	$R->startR;
	
	foreach my $pair ( @{$self->get_trx_pairs_aoh} )
	{
		my $q_trx = $self->{q_trx_hash}{ $pair->{Q_TRX_ID} };
		my $h_trx = $self->{h_trx_hash}{ $pair->{H_TRX_ID} };
		my $q_len = $q_trx->{SEQ_MASKED} ? $q_trx->{LEN_MASKED_WO_N} : length($q_trx->{SEQ});
		my $h_len = $h_trx->{SEQ_MASKED} ? $h_trx->{LEN_MASKED_WO_N} : length($h_trx->{SEQ});
		my $L1    = ($q_len < $h_len) ? $q_len : $h_len;
		my $L2    = ($q_len < $h_len) ? $h_len : $q_len;
		my $q_gc  = 100*compute_gc_content($q_trx->{SEQ});
		my $h_gc  = 100*compute_gc_content($h_trx->{SEQ});
		
		$pair->{SUM_ENERGY} = sum( map { $_->{DELTA_G} } values %{$pair->{ALL_SITES}} );
		$pair->{PVALUE}     = $self->_compute_pvalue_for_features($R,
			X        => -1*$pair->{SUM_ENERGY},
			LOG_L1   => log($L1),
			LOG_L2   => log($L2),
			GC_AVER  => ($q_gc + $h_gc)/2,
			GC_DIFF  => abs($q_gc - $h_gc),
			SCORE    => $self->{lastal_score},
		);
	}
	
	$R->stopR();
	
	# Compute the FDR adjusted pvalues: https://www.youtube.com/watch?v=K8LQSvtjcEo
	my @all_pairs = sort { $a->{PVALUE} <=> $b->{PVALUE} } @{$self->get_trx_pairs_aoh};
	for(my $i = $#all_pairs; $i >= 0; $i--)
	{
		my $pair = $all_pairs[$i];
		if( $i == $#all_pairs )    # The last (weakest) prediction
		{
			$pair->{PVALUE_ADJ} = $pair->{PVALUE};
		}
		else
		{
			my $padj = $all_pairs[$i]{PVALUE}*scalar(@all_pairs)/($i+1);
			$pair->{PVALUE_ADJ} = min($padj, $all_pairs[$i+1]{PVALUE_ADJ});
		}
	}
}

sub _compute_pvalue_for_features
{
	my $self = shift;
	my($R, %f) = @_;
	
	return 1 if $f{X} == 0;
	
	my @f_names = qw(LOG_L1   LOG_L2   GC_AVER   GC_DIFF   SCORE);
	my %range   = %{$self->{model}{feature_range}};
	foreach my $name ( @f_names )
	{
		if( $f{$name} < $range{$name}[0] || $range{$name}[1] < $f{$name} )
		{
			warn "WARNING: the value of the feature $name = $f{$name} is outside of the training range [$range{$name}[0], $range{$name}[1]]!";
		}
	}
	
	my $theta    = $self->{model}{theta};
	my @f_vals   = map { $f{$_} } @f_names;	
	my $logit_P  = vector_product([1, @f_vals], $theta->{logit_P} );
	my $log_mean = vector_product([1, @f_vals], $theta->{log_mean});
	my $log_var  = vector_product([1, @f_vals], $theta->{log_var} );
	
	my $P_non_zero  = 1/(1 + exp(-1*$logit_P));
	my $Gamma_mean  = exp( $log_mean );
	my $Gamma_var   = exp( $log_var );
	my $Gamma_shape = $Gamma_mean**2 / $Gamma_var;
	my $Gamma_rate  = $Gamma_mean    / $Gamma_var;
	
	$R->send( "cat($P_non_zero * pgamma($f{X}, shape=$Gamma_shape, rate=$Gamma_rate, lower.tail=FALSE))" );
	return $R->read;
}

sub run_lastal
{
	my $self = shift;
	my(%opts) = @_;
	
	$self->{trx_pairs} = {};
	print STDERR "[".localtime()."] Running LASTAL:\n" if $self->{verbose};
	foreach my $q_trx ( values %{$self->{q_trx_hash}} )
	{
		# Usually we have 1 (or few) query and many targets. LASTAL is optimized to align many reads on 1 genome.
		# So, let's create database from queries!
		my $q_seq = $q_trx->{SEQ_MASKED} || $q_trx->{SEQ};
		open(PIPE, "| lastdb $self->{lastdb_path}") or die "Can't open LASTDB pipe: $!";
		print PIPE ">$q_trx->{TRX_ID}\n$q_seq\n";
		close PIPE;
		
		my @all_targets = values %{$self->{h_trx_hash}};
		if( $self->{num_threads} == 1 )
		{
			print STDERR "[".localtime()."] $q_trx->{TRX_ID} vs ".scalar(@all_targets)." targets ...      " if $self->{verbose};
			my $hits = $self->_get_lastal_hits_for_trx_arr( \@all_targets );
			$self->_create_trx_pairs_from_lastal_hits($hits, %opts);
		}
		else
		{
			require Parallel::ForkManager;
			
			# https://perlmaven.com/speed-up-calculation-by-running-in-parallel
			my $pm = Parallel::ForkManager->new( $self->{num_threads} );
			$pm->run_on_finish( sub {
				my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hits) = @_;
				$self->_create_trx_pairs_from_lastal_hits($hits, %opts);
			});
			
			my($num_t, $i) = (scalar(@all_targets), 1);
			foreach my $t_trx ( @all_targets )
			{
				print STDERR "\r[".localtime()."] $q_trx->{TRX_ID} vs $t_trx->{TRX_ID}: ".($i++)."/$num_t           " if $self->{verbose};
				my $pid  = $pm->start and next;
				my $hits = $self->_get_lastal_hits_for_trx_arr( [$t_trx] );
				$pm->finish(0, $hits);
			}
			$pm->wait_all_children;
		}
		print STDERR "\n" if $self->{verbose};
	}
	
	# Make sure that the LASTAL sites do not overlap
	$self->_discard_overlapping_sites();
	
	my $trx_pairs_aoh = $self->get_trx_pairs_aoh();
	my $total_pairs   = scalar(@$trx_pairs_aoh);
	my $total_sites   = sum( map { scalar(keys %{$_->{ALL_SITES}}) } @$trx_pairs_aoh );
	print STDERR "[".localtime()."] Total number of transcript pairs after LASTAL search = $total_pairs ($total_sites sites)\n" if $self->{verbose};
}

sub run_RNAup
{
	my $self = shift;
	
	my @site_aoh = sort { $b->{SITE_LEN} <=> $a->{SITE_LEN} } @{$self->get_all_sites_aoh};
	return if scalar(@site_aoh) == 0;
	
	# Generate larger job's because RNAup usually takes <1 sec to  process one site,
	# but Parallel::ForkManager checks the jobs only one time per second
	my $reg_hash = $self->get_region_hash();
	my($job_arr, $site_i) = ([], 0);
	while( $site_i < scalar(@site_aoh) )
	{
		for(my $job_i = 0; $job_i < $self->{num_threads}; $job_i++)
		{
			my $site   = $site_aoh[ $site_i ];
			my $region = $reg_hash->{ $site->{SITE_ID} };
			die "Can't get region for site id '$site->{SITE_ID}'" if !$region;
			
			push @{$job_arr->[$job_i]}, {site => $site, region => $region};
			
			$site_i++;
			last if $site_i >= scalar(@site_aoh);
		}
	}
	
	if( $self->{num_threads} == 1 )
	{
		my $data = $self->_compute_energy_for_site_and_region_arr( $job_arr->[0] );
		$self->_save_site_energies( $data );
	}
	else
	{
		# https://perlmaven.com/speed-up-calculation-by-running-in-parallel
		require Parallel::ForkManager;
		
		my $pm = Parallel::ForkManager->new( $self->{num_threads} );
		$pm->run_on_finish( sub {
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data) = @_;
			$self->_save_site_energies( $data );
		});
		
		foreach my $job ( @$job_arr )
		{
			my $pid  = $pm->start and next;
			my $data = $self->_compute_energy_for_site_and_region_arr( $job );
			$pm->finish(0, $data);
		}
		$pm->wait_all_children;
	}
	print STDERR "\n" if $self->{verbose};
}

sub get_region_hash
{
	my $self = shift;
	
	my $reg_hash = {};
	foreach my $trx_pair ( @{$self->get_trx_pairs_aoh} )
	{
		foreach my $site ( values %{$trx_pair->{ALL_SITES}} )
		{
			my $q_trx  = $self->{q_trx_hash}{$trx_pair->{Q_TRX_ID}};
			my $h_trx  = $self->{h_trx_hash}{$trx_pair->{H_TRX_ID}};
			my %region = ();
			$region{SITE_ID}   = $site->{SITE_ID};
			$region{Q_TRX_ID}  = $trx_pair->{Q_TRX_ID};
			$region{H_TRX_ID}  = $trx_pair->{H_TRX_ID};
			
			if( $site->{SITE_LEN} > $self->{model}{max_w} )    # Do not add flanks to long sites  (e.g. length > 50)
			{
				$region{Q_START}   = $site->{Q_START};
				$region{Q_END}     = $site->{Q_END};
				$region{H_START}   = $site->{H_START};
				$region{H_END}     = $site->{H_END};
			}
			else
			{
				my $q_len    = length( $q_trx->{SEQ} );
				my $h_len    = length( $h_trx->{SEQ} );
				my $flank    = $site->{SITE_LEN};
				$region{Q_START}   = max(1,      $site->{Q_START} - $flank);
				$region{Q_END}     = min($q_len, $site->{Q_END}   + $flank);
				$region{H_START}   = max(1,      $site->{H_START} - $flank);
				$region{H_END}     = min($h_len, $site->{H_END}   + $flank);
			}
			
			$region{Q_SEQ} = substr($q_trx->{SEQ}, $region{Q_START}-1, $region{Q_END}-$region{Q_START}+1);
			$region{H_SEQ} = substr($h_trx->{SEQ}, $region{H_START}-1, $region{H_END}-$region{H_START}+1);
			
			$reg_hash->{$site->{SITE_ID}} = \%region;
		}
	}
	
	return $reg_hash;
}

sub get_all_sites_aoh
{
	my $self = shift;
	
	my $sites_aoh = [];
	foreach my $trx_pair ( @{$self->get_trx_pairs_aoh} )
	{
		push @$sites_aoh, values %{$trx_pair->{ALL_SITES}};
	}
	
	return [sort { $a->{SITE_ID} <=> $b->{SITE_ID} } @$sites_aoh];
}

sub get_trx_pairs_aoh
{
	my $self = shift;
	
	my($pairs_hash, $pairs_aoh) = ($self->{trx_pairs}, []);
	foreach my $q_trx_id ( keys %$pairs_hash )
	{
		foreach my $h_trx_id ( keys %{$pairs_hash->{$q_trx_id}} )
		{
			push @$pairs_aoh, $pairs_hash->{$q_trx_id}{$h_trx_id};
		}
	}
	
	return $pairs_aoh;
}

###
# PRIVATE METHODS

sub _discard_overlapping_sites
{
	my $self = shift;
	
	foreach my $trx_pair ( @{$self->get_trx_pairs_aoh} )
	{
		my $no_overlap = [];
		foreach my $site ( sort { $b->{SCORE} <=> $a->{SCORE} } values %{$trx_pair->{ALL_SITES}} )
		{
			push @$no_overlap, $site if !has_double_overlap_with_other_sites($site, $no_overlap);
		}
		$trx_pair->{ALL_SITES} = { map { $_->{SITE_ID} => $_ } @$no_overlap };
	}
}

sub _save_site_energies
{
	my $self = shift;
	my($all_data) = @_;
	foreach my $d ( @$all_data )
	{
		my $q_trx_id = $d->{region}{Q_TRX_ID};
		my $h_trx_id = $d->{region}{H_TRX_ID};
		my $site_id  = $d->{site}{SITE_ID};
		my $delta_G  = $d->{delta_G};
		if( $d->{site}{SITE_LEN} > $self->{model}{max_w} )
		{
			# Approximate the energy of a long site by increasing its energy
			# according to how much it is longer than the max_w RNAup parameter
			$delta_G *= $d->{site}{SITE_LEN} / $self->{model}{max_w}
		}
		$self->{trx_pairs}{$q_trx_id}{$h_trx_id}{ALL_SITES}{$site_id}{DELTA_G} = $delta_G;
		$self->{trx_pairs}{$q_trx_id}{$h_trx_id}{ALL_SITES}{$site_id}{RNAUP}   = $d->{rnaup};
	}
}

sub _compute_energy_for_site_and_region_arr
{
	my $self = shift;
	my($site_region_aoh) = @_;
	
	my($max_w, $data) = ($self->{model}{max_w}, []);
	for(my $i = 0; $i < scalar(@$site_region_aoh); $i++)
	{
		my $percent = sprintf('%.1f', 100*$i/scalar(@$site_region_aoh));
		print STDERR "\r[".localtime()."] Computing site energies: $percent%            " if $self->{verbose};
		
		my $site   = $site_region_aoh->[$i]{site};
		my $region = $site_region_aoh->[$i]{region};
		
		my $prm_w = $site->{SITE_LEN} < $max_w ? $site->{SITE_LEN} : $max_w;
		my $rnaup = _run_RNAup_for_region($region, $prm_w);
		
		my $delta_G = 0;
		if( intersect_digits($site->{Q_START}, $site->{Q_END}, $rnaup->{Q_START}, $rnaup->{Q_END})
			&&
			intersect_digits($site->{H_START}, $site->{H_END}, $rnaup->{H_START}, $rnaup->{H_END}) )
		{
			# RNAup duplex matches the LASTAL site
			$delta_G = $rnaup->{DELTA_G};
		}
		
		push @$data, {delta_G => $delta_G, rnaup => $rnaup, site => $site, region => $region};
	}

	return $data;
}

sub _create_trx_pairs_from_lastal_hits
{
	my $self = shift;
	my($lastal_hits, %opts) = @_;
	
	foreach my $q_trx_id ( keys %{$self->{q_trx_hash}} )
	{
		foreach my $h_trx_id ( keys %{$self->{h_trx_hash}} )
		{
			$self->{trx_pairs}{$q_trx_id}{$h_trx_id} = {
				Q_TRX_ID  => $q_trx_id,   H_TRX_ID  => $h_trx_id,   ALL_SITES => {},
			} if !exists $self->{trx_pairs}{$q_trx_id}{$h_trx_id};
		}
	}
	
	foreach my $lastal ( @$lastal_hits )
	{
		my $site_len = length( $lastal->{q_alignment} );
		next if $opts{min_site_len} && $site_len < $opts{min_site_len};
		next if $opts{max_site_len} && $site_len > $opts{max_site_len};
		
		my $q_trx_id = $lastal->{q_name};
		my $h_trx_id = $lastal->{h_name};
		warn "Unknown seqname in the LASTAL output: '$q_trx_id'" and next if !$self->{q_trx_hash}{$q_trx_id};
		warn "Unknown seqname in the LASTAL output: '$h_trx_id'" and next if !$self->{h_trx_hash}{$h_trx_id};
		
		my $h_len   = length($self->{h_trx_hash}{$h_trx_id}{SEQ});
		my $h_start = $h_len - $lastal->{h_alnSize} - $lastal->{h_start};
		
		my $q_seq_chunk  = substr($self->{q_trx_hash}{$q_trx_id}{SEQ}, $lastal->{q_start}, $lastal->{q_alnSize});
		(my $q_ali_chunk = $lastal->{q_alignment}) =~ s/\-//g;
		warn "Wrong coordinates in LASTAL output!" and next if uc($q_seq_chunk) ne uc($q_ali_chunk);
		
		my $h_seq_chunk = substr($self->{h_trx_hash}{$h_trx_id}{SEQ}, $h_start, $lastal->{h_alnSize});
		(my $h_ali_chunk = revcomp($lastal->{h_alignment})) =~ s/\-//g;
		warn "Wrong coordinates in LASTAL output:\n$h_seq_chunk\n$h_ali_chunk\n" and next if uc($h_seq_chunk) ne uc($h_ali_chunk);
		
		my $site_id = $self->{unique_id}++;
		$self->{trx_pairs}{$q_trx_id}{$h_trx_id}{ALL_SITES}{$site_id} = {
			SITE_ID  => $site_id,
			Q_START  => $lastal->{q_start} + 1,
			Q_END    => $lastal->{q_start} + $lastal->{q_alnSize},
			Q_LEN    => $lastal->{q_alnSize},
			H_START  => $h_start + 1,
			H_END    => $h_start + $lastal->{h_alnSize},
			H_LEN    => $lastal->{h_alnSize},
			SITE_LEN => $site_len,
			SITE_CMPL=> compute_ali_identity($lastal->{q_alignment}, $lastal->{h_alignment}),
			SCORE    => $lastal->{score},
		};
	}
}

sub _get_lastal_hits_for_trx_arr
{
	my $self = shift;
	my($trx_arr) = @_;
	
	my $params = "-s 0 -m 9999999 -p $self->{lastal_mtx} -a $self->{lastal_a} -b $self->{lastal_b} -e $self->{lastal_score}";
	my $pid = open2(*CHLD_OUT, *CHLD_IN, "lastal $params $self->{lastdb_path}");
	print CHLD_IN ">$_->{TRX_ID}\n$_->{SEQ}\n" foreach @$trx_arr;
	close CHLD_IN;
	
	my @head = qw(name start alnSize strand seqSize alignment);
	my($all_hits, $n_lines) = ([], 0);
	while( my $line = <CHLD_OUT> )
	{
		$n_lines++;
		next if $line =~ /^#/;
		if( $line =~ /^a / )
		{
			# a score=47 EG2=3.8e+09 E=0.02
			# s KU881769.1  2506 17 + 2882 TGCATCCTGGACCCCAG
			# s NM_002229.2  436 17 - 1832 TGCTCCCTGGACCCCAG
			my %info = map { /(\S+)\=(\S+)/; $1 => $2; } grep { /\=/ } split /\s+/, $line;
			
			$line = <CHLD_OUT>;
			warn "Wrong LASTAL format: '$line'" and next unless $line =~ s/^s //;
			my @q_vals = split /\s+/, $line;
			$info{"q_$head[$_]"} = $q_vals[$_] for 0..$#head;
			
			$line = <CHLD_OUT>;
			warn "Wrong LASTAL format: '$line'" and next unless $line =~ s/^s //;
			my @h_vals = split /\s+/, $line;
			$info{"h_$head[$_]"} = $h_vals[$_] for 0..$#head;
			
			push @$all_hits, \%info;
		}
	}
	close CHLD_OUT;
	
	waitpid( $pid, 0 );   # to avoid zombie processes
	
	warn "WARNING: the LASTAL output was empty!" if $n_lines == 0;	
	
	return $all_hits;
}

###
# SUBROUTINES

sub _run_RNAup_for_region
{
	my($region, $param_w) = @_;
	
	my $pid = open2(*CHLD_OUT, *CHLD_IN, "RNAup -o -b -w $param_w 2> /dev/null");
	print CHLD_IN "$region->{Q_SEQ}&$region->{H_SEQ}";
	close CHLD_IN;
	
	my @all_lines = ();
	while (<CHLD_OUT>)
	{
		s/[\n\r]//g;
		push @all_lines, $_;
	}
	close CHLD_OUT;
	waitpid( $pid, 0 );   # to avoid zombie processes
	
	warn "Something is wrong witn RNAup output: ".Dumper(\@all_lines) and return(undef) if scalar(@all_lines) != 2;
	
	# $lines[0] = '((((..((((((((&)))))))).))))  18,31  :   1,13  (-10.68 = -26.08 + 6.12 + 9.28)'
	# $lines[1] = 'CACCCCCGCCCGGC&GCCGGGCGCGGUG'
	if( $all_lines[0] =~ /^(\S+)\s+(\d+)\,(\d+)\s+\:\s+(\d+)\,(\d+)\s*\((\S+)/ )
	{
		my $rnaup = {DOT_BRACKET => $1, DELTA_G => $6, ALL_LINES => \@all_lines};
		if( length($region->{Q_SEQ}) >= length($region->{H_SEQ}) )
		{
			$rnaup->{Q_START} = $region->{Q_START} - 1 + $2;
			$rnaup->{Q_END}   = $region->{Q_START} - 1 + $3;
			$rnaup->{H_START} = $region->{H_START} - 1 + $4;
			$rnaup->{H_END}   = $region->{H_START} - 1 + $5;
			($rnaup->{Q_SEQ}, $rnaup->{H_SEQ}) = split /\&/, $all_lines[1];
		}
		else
		{
			# The longer sequence ALWAYS goes first
			$rnaup->{Q_START} = $region->{Q_START} - 1 + $4;
			$rnaup->{Q_END}   = $region->{Q_START} - 1 + $5;
			$rnaup->{H_START} = $region->{H_START} - 1 + $2;
			$rnaup->{H_END}   = $region->{H_START} - 1 + $3;
			($rnaup->{H_SEQ}, $rnaup->{Q_SEQ}) = split /\&/, $all_lines[1];
		}
		$rnaup->{Q_SEQ} =~ tr/Uu/Tt/;
		$rnaup->{H_SEQ} =~ tr/Uu/Tt/;
		
		warn "Can't find RNAup seq in the region:\n'$rnaup->{Q_SEQ}'\n'$region->{Q_SEQ}'" and return(undef) if $region->{Q_SEQ} !~ /$rnaup->{Q_SEQ}/;
		warn "Can't find RNAup seq in the region:\n'$rnaup->{H_SEQ}'\n'$region->{H_SEQ}'" and return(undef) if $region->{H_SEQ} !~ /$rnaup->{H_SEQ}/;
		
		return $rnaup;
	}
	else
	{
		warn "Can't parse RNAup output: '$all_lines[0]'";
		return undef;
	}
}

sub has_double_overlap_with_other_sites
{
	my($site, $site_arr) = @_;
	
	# check overlaps on the query AND on the hit
	foreach my $s ( @$site_arr )
	{
		if( $s->{Q_TRX_ID} && $site->{Q_TRX_ID} && $s->{H_TRX_ID} && $site->{H_TRX_ID} )
		{
			next if $s->{Q_TRX_ID} ne $site->{Q_TRX_ID} || $s->{H_TRX_ID} ne $site->{H_TRX_ID};
		}
		
		if( intersect_digits($site->{Q_START},$site->{Q_END},$s->{Q_START},$s->{Q_END})
		    &&
		    intersect_digits($site->{H_START},$site->{H_END},$s->{H_START},$s->{H_END}))
		{
			return $s;
		}
	}
	
	return undef;
}


1;

