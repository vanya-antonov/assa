package ASSA::Lib;

use strict;
use warnings;

# $Id$

use Data::Dumper;
use Exporter;
use vars qw(@ISA @EXPORT_OK);
@ISA       = qw(Exporter);
@EXPORT_OK = qw(
	fasta2hash
	compute_gc_content
	compute_ali_identity
	max
	min
	read_file_to_str
	sum
	vector_product

	ah2a
	array2hash
	find_longest_orfs
	intersect_digits
	log10
	read_table_wo_header
	read_table_wo_header_from_fh
	revcomp
	unique
	write_to_file
);


###
# SUBROUTINES

sub fasta2hash
{
    my($fn) = shift;
	
	open(F, '<', $fn) or die "Can't open '$fn': $!";
	
	my($seqs_hash, $cur_seq) = ({}, {});
	while( my $line = <F> )
	{
		$line =~ s/[\n\r]//g;
		$line =~ s/^\s*(.*?)\s*$/$1/;   # trim spaces on both sides
		
		if($line =~ /^>/)
		{
			if( $cur_seq->{SEQ} )
			{
				warn "ID '$cur_seq->{TRX_ID}' is duplicated!" if $seqs_hash->{$cur_seq->{TRX_ID}};
				$seqs_hash->{$cur_seq->{TRX_ID}} = $cur_seq;
			}
			# We use the idea of the TRX_ID because LASTAL uses it anyway
			my($full_name, $trx_id) = $line =~ /^>((\S+).*)/;
			$cur_seq = {FULL_NAME => $full_name, TRX_ID => $trx_id, SEQ => ''};
		}
		else
		{
			$cur_seq->{SEQ} .= $line;
		}
	}
	if( $cur_seq->{SEQ} )
	{
		warn "TRX_ID '$cur_seq->{TRX_ID}' is duplicated!" if $seqs_hash->{$cur_seq->{TRX_ID}};
		$seqs_hash->{$cur_seq->{TRX_ID}} = $cur_seq;
	}
	
	close F;
	
	return $seqs_hash;
}

sub revcomp
{
	my($seq) = @_;
	$seq = reverse $seq;
	$seq =~ tr/ACGTacgt/TGCAtgca/;
	return $seq;
}

sub compute_gc_content
{
	my($seq) = @_;
	my $num_gc = $seq =~ tr/CGcg/CGcg/;
	return sprintf('%.5f', $num_gc/length($seq));
	
}

sub compute_ali_identity
{
	my($seq1, $seq2) = @_;
	
	warn "Sequences have different lengths!" and return(undef) if length($seq1) != length($seq2);
	
	my @lets1 = split //, uc($seq1);
	my @lets2 = split //, uc($seq2);
	
	my $num_matches = 0;
	foreach my $i ( 0..$#lets1 )
	{
		$num_matches++ if $lets1[$i] eq $lets2[$i];
	}
	
	return sprintf('%.5f', $num_matches/length($seq1))
}

sub intersect_digits
{
	#################################
	# Example:
	# 
	# |            INPUT            |          RETURNS            |
	# |-----------------------------|-----------------------------|
	# |                             |                             |
	# |  s1=3      e1=14            |                             |
	# |     |----------|            |        6 (=14-9+1)          |
	# |           |---------|       |                             |
	# |        s2=9     e2=19       |                             |
	# |                             |                             |
	# =============================================================
	#
	my $s1 = shift;   # start1
	my $e1 = shift;   # end1
	my $s2 = shift;   # start2
	my $e2 = shift;   # end2

	confess("Wrong params: s1='$s1', e1='$e1', s2='$s2', e2='$e2'") if !defined $s1 || !defined $e1 || !defined $s2 || !defined $e2;

	if($s1<=$s2 && $e1>=$e2)
	{
		# 
		# |----------------|
		#     |--------|
		# 
		return $e2-$s2+1;
	}
	if($s1>=$s2 && $e1<=$e2)
	{
		# 
		#     |--------|
		# |----------------|
		# 
		return $e1-$s1+1;
	}
	elsif($s1<=$s2 && $e1>=$s2 && $e1<=$e2)
	{
		# 
		# |------------|
		#       |-------------|
		# 
		return $e1-$s2+1;
	}
	elsif($s1>=$s2 && $s1<=$e2 && $e1>=$e2)
	{
		# 
		#       |-------------|
		# |------------|
		# 
		return $e2-$s1+1;
	}
	else
	{
		return 0;
	}
}

sub min
{
	my @arr = @_;
	return undef if @arr == 0;
	my $min = $arr[0];
	for( 1 .. $#arr )
	{
		$min = $arr[$_] if $arr[$_] < $min;
	}
	return $min;
}

sub max
{
	my @arr = @_;
	return undef if @arr == 0;
	my $max = $arr[0];
	for( 1 .. $#arr ){
		$max = $arr[$_] if $arr[$_] > $max;
	}
	return $max;
}

sub sum
{
	my $res = 0;
	$res += $_ foreach @_;
	return $res;
}

sub vector_product
{
	my($v1, $v2) = @_;
	
	die "The two vectors are of different sizes!" if scalar(@$v1) != scalar(@$v2);
	
	my $product = 0;
	for(my $i=0; $i < scalar(@$v1); $i++)
	{
		$product += $v1->[$i] * $v2->[$i];
	}
	
	return $product;
}

sub read_file_to_str
{
	my($fn) = @_;
	
	open(F, $fn eq '-' ? '<&STDIN' : $fn) or die "Can't open file '$fn': $!";
	my @lines = <F>;
	close F;
	
	return join('', @lines);
}













sub read_table_wo_header
{
	my($table_fn, %opts) = @_;
	
	open(my $fh, $table_fn eq '-' ? '<&STDIN' : $table_fn) or die "Can't open file '$table_fn': $!";
	my $tbl = read_table_wo_header_from_fh($fh,%opts);
	close $fh;
	
	return $tbl;
}

###
# %opts
#   |
#   |--->{delim_re}
#   |--->{strict_t}
#   |--->{skip_head}
#   |--->{max_lines} -- max number of lines to read (excluding skip_head if present)
#
sub read_table_wo_header_from_fh
{
	my($fh, %opts) = @_;
	my $delim_re = defined $opts{delim_re} ? $opts{delim_re} : qr/\t/;
	my $strict_t = defined $opts{strict_t} ? $opts{strict_t} : 1;
	if( $opts{skip_head} )
	{
		$_ = <$fh> for 1..$opts{skip_head};
	}
	
	my($num_cols, $tbl) = (undef, []);
	while ( my $str = <$fh> )
	{
		next if $str =~ /^\s*$/;
		$str =~ s/[\n\r]//g;
		$str =~ s/^$delim_re//;   # remove delimiter from the beginning to avoid empty values
		my @vals = split $delim_re, $str;
		$num_cols = scalar @vals if !defined $num_cols;
		die "Wrong number of columns: '$str'" if $strict_t && @vals != $num_cols;
		push @$tbl, \@vals;
		
		last if $opts{max_lines} && scalar(@$tbl) >= $opts{max_lines};
	}
	
	return $tbl;
}

sub write_to_file
{
	my($txt, $out_fn, %opts) = @_;
	
	open(my $fh, '>', $out_fn) or die("Can't open file '$out_fn' to write: $!");
	binmode $fh if $opts{binary};
	print $fh $txt;
	close $fh;
	
#	chmod($opts{mode}, $out_fn) if defined $opts{mode};
	`chmod $opts{mode} $out_fn` if defined $opts{mode};
}

sub ah2a
{
	my($key, $aoh) = @_;
	return [ map { $_->{$key} } @$aoh ];
}

sub unique
{
	return keys %{array2hash(@_)};
}

sub array2hash
{
	my %hash = ();
	$hash{$_}++ foreach @_;
	return \%hash;
}

###
# Arguments:
#     $seq  -- sequence to find ORFs in
#     %opts -- options
#       |
#       |--->{min_len}  --  minimum lenght of the ORF
#       |--->{verbose}  --  boolean
# 
# Returns:
#      $orfs      --  reference to array of hashes
#      $orfs[i]
#        |
#        |--->{start}  -- the ORF at the very beginning of $seq has start = 0
#        |--->{len}
# 
sub find_longest_orfs
{
	my($seq, %opts) = @_;
	
	my $full_len = length($seq);   # Length of the initial full sequence
	my $offset   = 0;              # Offset from the start of initial sequence
	my $orfs     = []; 
	while( $seq =~ /ATG((?:[ACGT]{3})+?(?:TAG|TAA|TGA))/i )
	{
		$seq = $1.$';   # In the next iteration we will search the same seq but without the start codon
		$offset += length($`)+3;
		
		push @$orfs, {
			seq   => $&,
			start => $offset-3,
			len   => length($&)
		};
		
		printf STDERR ("\r\t%d ORFs found in %.1f %% of sequence", scalar(@$orfs), 100-100*length($seq)/$full_len) if $opts{verbose};
	}
	print STDERR "\n" if $opts{verbose};
	
	@$orfs = grep { $_->{len} > $opts{min_len} } @$orfs if $opts{min_len};
	
	return [sort { $b->{len} <=> $a->{len} } @$orfs];
}

sub log10
{
	return log($_[0])/log(10);
}

1;

