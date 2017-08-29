
=head1 NAME

Bio::Pfam::Scan::PfamScan

=cut

package Bio::Pfam::Scan::PfamScan;

=head1 SYNOPSIS

  my $ps = Bio::Pfam::Scan::PfamScan->new(
             -cut_off => $hmmscan_cut_off,
             -dir => $dir,
             -clan_overlap => $clan_overlap,
             -fasta => $fasta,
             -align => $align,
             -as => $as
           );

  $ps->search;
  $ps->write_results;

=head1 DESCRIPTION

$Id: PfamScan.pm,v 1.4 2010-01-12 09:41:42 jm14 Exp $

=cut

use strict;
use warnings;

use Bio::Pfam::HMM::HMMResultsIO;
use Bio::Pfam::Active_site::as_search;
use Bio::SimpleAlign;
use Bio::Pfam::Scan::Seq;

use Carp;
use IPC::Run qw( start finish );

#-------------------------------------------------------------------------------
#- constructor -----------------------------------------------------------------
#-------------------------------------------------------------------------------

=head1 METHODS

=head2 new

The only constructor for the object. Accepts a set of arguments that specify
the parameters for the search:

=over

=item -cut_off

=item -dir

=item -clan_overlap

=item -fasta

=item -sequence

=item -align

=item -hmm

=item -as

=back

=cut

sub new {
  my ( $class, @args ) = @_;

  my $self = {};
  bless $self, $class;

  # To avoid hard coding the location for the binary, we assume it will be on the path.....
  $self->{_HMMSCAN} = 'hmmscan';

  # handle arguments, if we were given any here
  $self->_process_args(@args) if @args;

  return $self;
}

#-------------------------------------------------------------------------------
#- public methods --------------------------------------------------------------
#-------------------------------------------------------------------------------

=head2 search

The main method on the object. Performs a C<hmmscan> search using the supplied
sequence and the specified HMM library.

=cut

sub search {
  my ( $self, @args ) = @_;

  # handle the arguments, if we were handed any here
  $self->_process_args(@args) if @args;

  # set up the output header
  $self->_build_header;

  croak qq(FATAL: no sequence given; set the search parameters before calling "search")
    unless defined $self->{_sequence};

  my ( %AllResults, $pfamB, $firstResult );

  foreach my $hmmlib ( @{ $self->{_hmmlib} } ) {

    my ( @hmmscan_cut_off, $seq_evalue, $dom_evalue );
    if ( $hmmlib !~ /Pfam\-B/ ) {
      @hmmscan_cut_off = @{ $self->{_hmmscan_cutoff} };
    }
    else {
      $pfamB      = 1;
      $seq_evalue = 0.001;
      $dom_evalue = 0.001;

      # It's a pfamB search so use some default cut off values
      push @hmmscan_cut_off, '-E', $seq_evalue, '--domE', $dom_evalue;
    }

    push @{ $self->{_header} },
      "#     cpu number specified: " . $self->{_cpu} . "\n"
      if ( $hmmlib !~ /Pfam\-B/ and $self->{_cpu} );

    push @{ $self->{_header} },
      "#        searching against: "
      . $self->{_dir}
      . "/$hmmlib, with cut off "
      . join( " ", @hmmscan_cut_off ) . "\n";
    my @params;
    if ( $self->{_cpu} ) {
      @params = (
        'hmmsearch', '--notextw', '--cpu', $self->{_cpu}, @hmmscan_cut_off,
        $self->{_dir} . '/' . $hmmlib,
        $self->{_fasta}
      );
    }
    else {
      @params = (
        'hmmsearch', '--notextw', @hmmscan_cut_off, $self->{_dir} . '/' . $hmmlib,
        $self->{_fasta}
      );

    }

    print STDERR "PfamScan::search: hmmscan command: |@params|\n"
      if $ENV{DEBUG};
    print STDERR 'PfamScan::search: sequence: |' . $self->{_sequence} . "|\n"
      if $ENV{DEBUG};

    my $run = start \@params, '<pipe', \*IN, '>pipe', \*OUT, '2>pipe', \*ERR
      or croak qq(FATAL: error running hmmscan; IPC::Run returned '$?');

    # print IN $self->{_sequence}; ;
    close IN;

    $self->{_hmmresultIO} = Bio::Pfam::HMM::HMMResultsIO->new;
    $self->{_all_results} = $self->{_hmmresultIO}->parseMultiHMMER3( \*OUT );
    close OUT;

    $self->{_all_results} = $self->_convert_results_search_to_scan($self->{_all_results});

    my $err;
    while (<ERR>) {
      $err .= $_;
    }
    close ERR;

    finish $run
      or croak qq|FATAL: error running hmmscan ($err); ipc returned '$?'|;

    unless ( $hmmlib =~ /Pfam\-B/ ) {

      if ( $self->{_clan_overlap} ) {
        push( @{ $self->{_header} }, "#    resolve clan overlaps: off\n" );
      }
      else {
        push( @{ $self->{_header} }, "#    resolve clan overlaps: on\n" );
        $self->_resolve_clan_overlap;
      }

      if ( $self->{_as} ) {
        push( @{ $self->{_header} }, "#     predict active sites: on\n" );
        $self->_pred_act_sites;
      }
      else {
        push( @{ $self->{_header} }, "#     predict active sites: off\n" );
      }

      if ( $self->{_translate} ) {
        push @{ $self->{_header} },  "#   translate DNA sequence: " . $self->{_translate} . "\n";
      }
    }

    # Determine which hits are significant
    foreach my $result ( @{ $self->{_all_results} } ) {
      foreach
        my $unit ( sort { $a->seqFrom <=> $b->seqFrom } @{ $result->units } )
      {

        unless ($pfamB) {

          $unit->sig(0);
          if ( $result->seqs->{ $unit->name }->bits >=
            $self->{_seqGA}->{ $unit->name } )
          {
            if ( $unit->bits >= $self->{_domGA}->{ $unit->name } ) {
              $unit->sig(1);
            }
          }
        }
      }
    }

    if ($firstResult) {
      $AllResults{ $self->{_all_results} } = $self->{_all_results};
    }
    else {
      $firstResult = $self->{_all_results};
    }

  }    # end of "foreach $hmmlib"

  # If more than one search, merge results into one object
  if ( keys %AllResults ) {

    foreach my $AllResult ( keys %AllResults ) {

      foreach my $seq_id ( keys %{ $self->{_seq_hash} } ) {

        my $flag;

        #If seq exists in both, add all units from $AllResult to $firstResult
        foreach my $result ( @{$firstResult} ) {

          if ( $result->seqName eq $seq_id ) {
            $flag = 1;

            foreach my $result2 ( @{ $AllResults{$AllResult} } ) {

              if ( $result2->seqName eq $seq_id ) {
                foreach my $hmmname ( keys %{ $result2->seqs } ) {
                  $result->addHMMSeq( $result2->seqs->{$hmmname} );
                }
                foreach my $unit ( @{ $result2->units } ) {
                  $result->addHMMUnit($unit);
                }
              }
            }
          }
        }

        #If seq doesn't exist in $firstResult, need to add both sequence and units to $firstResult
        unless ($flag) {
          foreach my $result2 ( @{ $AllResults{$AllResult} } ) {
            if ( $result2->seqName eq $seq_id ) {
              push @{$firstResult}, $result2;
            }
          }
        }
      }
    }
    $self->{_all_results} = $firstResult;

  }    # end of "if keys %AllResults"

  push @{ $self->{_header} }, "# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =\n#\n";

  if ( $self->{_as} ) {
    push @{ $self->{_header} }, "# <seq id> <alignment start> <alignment end> <envelope start> <envelope end> <hmm acc> <hmm name> <type> <hmm start> <hmm end> <hmm length> <bit score> <E-value> <significance> <clan> <predicted_active_site_residues>";
  }
  else {
    push @{ $self->{_header} }, "# <seq id> <alignment start> <alignment end> <envelope start> <envelope end> <hmm acc> <hmm name> <type> <hmm start> <hmm end> <hmm length> <bit score> <E-value> <significance> <clan>";
  }

  if ( $self->{_translate} ) {
    push @{ $self->{_header} }, " <strand> <nt start> <nt end>";
  }
  push @{ $self->{_header} }, "\n";
}

#-------------------------------------------------------------------------------

=head2 write_results

Writes the results of the C<hmmscan> search. Takes a single argument, which can
be an open filehandle or a filename. A fatal error is generated if a file of the
given name already exists.

=cut

sub write_results {
  my ( $self, $out, $e_seq, $e_dom, $b_seq, $b_dom ) = @_;

  my $fh;

  if ( ref $out eq 'GLOB' ) {

    # we were handed a filehandle
    $fh = $out;
  }
  elsif ( $out and not ref $out ) {

    # we were handed a filename
    croak qq(FATAL: output file "$out" already exists) if -f $out;

    open( FH, ">$out" )
      or croak qq(FATAL: Can\'t write to your output file "$out": $!);
    $fh = \*FH;
  }
  else {

    # neither filehandle nor filename, default to STDOUT
    $fh = \*STDOUT;
  }

  if ( $self->{_header} ) {
    my $header = join '', @{ $self->{_header} };
    print $fh "$header\n";
  }

  foreach my $result ( @{ $self->{_all_results} } ) {
    $self->{_hmmresultIO}
      ->write_ascii_out( $result, $fh, $self, $e_seq, $e_dom, $b_seq, $b_dom );
  }
  close $fh;
}

#-------------------------------------------------------------------------------

=head2 results

Returns the search results.

=cut

sub results {
  my ( $self, $e_value ) = @_;

  unless ( defined $self->{_all_results} ) {
    carp "WARNING: call search() before trying to retrieve results";
    return;
  }

  my @search_results = ();

  foreach my $hmm_result ( @{ $self->{_all_results} } ) {
    push @search_results, @{ $hmm_result->results( $self, $e_value ) };
  }

  return \@search_results;
}

#-------------------------------------------------------------------------------
#- private methods -------------------------------------------------------------
#-------------------------------------------------------------------------------

=head1 PRIVATE METHODS

=head2 _process_args

Handles the input arguments.

=cut

sub _process_args {
  my ( $self, @args ) = @_;

  # accept both a hash and a hash ref
  my $args = ( ref $args[0] eq 'HASH' ) ? shift @args : {@args};

  # make sure we get a sequence
  if ( $args->{-fasta} and $args->{-sequence} ) {
    croak qq(FATAL: "-fasta" and "-sequence" are mutually exclusive);
  }
  elsif ( $args->{-fasta} ) {
    croak qq(FATAL: fasta file "$args->{-fasta}" doesn\'t exist)
      unless -s $args->{-fasta};
  }
  elsif ( $args->{-sequence} ) {
    croak qq(FATAL: no sequence given)
      unless length( $args->{-sequence} );
  }
  else {
    croak qq(FATAL: must specify either "-fasta" or "-sequence");
  }

  # check the cut off
  if ( ( $args->{-e_seq} and ( $args->{-b_seq} || $args->{-b_dom} ) )
    or ( $args->{-b_seq} and ( $args->{-e_seq} || $args->{-e_dom} ) )
    or ( $args->{-b_dom} and $args->{-e_dom} ) )
  {
    croak qq(FATAL: can\'t use e value and bit score threshold together);
  }

  $self->{_hmmscan_cutoff} = ();
  if ( $args->{-e_seq} ) {
    croak qq(FATAL: the E-value sequence cut-off "$args->{-e_seq}" must be a positive non-zero number)
      unless $args->{-e_seq} > 0;

    push @{ $self->{_hmmscan_cutoff} }, '-E', $args->{-e_seq};
  }

  if ( $args->{-e_dom} ) {
    croak q(FATAL: if you supply "-e_dom" you must also supply "-e_seq")
      unless $args->{-e_seq};

    croak qq(FATAL: the E-value domain cut-off "$args->{-e_dom}" must be positive non-zero number)
      unless $args->{-e_dom} > 0;

    push @{ $self->{_hmmscan_cutoff} }, '--domE', $args->{-e_dom};
  }

  if ( $args->{-b_seq} ) {
    push @{ $self->{_hmmscan_cutoff} }, '-T', $args->{-b_seq};
  }

  if ( $args->{-b_dom} ) {
    croak q(FATAL: if you supply "-b_dom" you must also supply "-b_seq")
      unless $args->{-b_seq};

    push @{ $self->{_hmmscan_cutoff} }, '--domT', $args->{-b_dom};
  }

  unless ( $self->{_hmmscan_cutoff} ) {
    push @{ $self->{_hmmscan_cutoff} }, '--cut_ga';
  }

  # make sure we have a valid directory for the HMM data files
  croak qq(FATAL: directory "$args->{-dir}" does not exist)
    unless -d $args->{-dir};

  # populate the object
  $self->{_cut_off}      = $args->{-cut_off};
  $self->{_dir}          = $args->{-dir};
  $self->{_clan_overlap} = $args->{-clan_overlap};
  $self->{_fasta}        = $args->{-fasta};
  $self->{_align}        = $args->{-align};
  $self->{_as}           = $args->{-as};
  $self->{_sequence}     = $args->{-sequence};
  $self->{_cpu}          = $args->{-cpu};
  $self->{_translate}    = $args->{-translate};

  $self->{_hmmlib} = [];
  if ( $args->{-hmmlib} ) {
    if ( ref $args->{-hmmlib} eq 'ARRAY' ) {
      push @{ $self->{_hmmlib} }, @{ $args->{-hmmlib} };
    }
    else {
      push @{ $self->{_hmmlib} }, $args->{-hmmlib};
    }
  }
  else {
    push @{ $self->{_hmmlib} }, "Pfam-A.hmm";
  }

  # Now check that the library exists in the data dir!
  foreach my $hmmlib ( @{ $self->{_hmmlib} } ) {

    croak qq(FATAL: can't find $hmmlib and/or $hmmlib binaries in "$args->{-dir}")
      unless (
      -s $self->{_dir},
      "/$hmmlib"
      and -s $self->{_dir} . "/$hmmlib.h3f"
      and -s $self->{_dir} . "/$hmmlib.h3i"
      and -s $self->{_dir} . "/$hmmlib.h3m"
      and -s $self->{_dir} . "/$hmmlib.h3p"
      and -s $self->{_dir} . "/$hmmlib.dat"
      );

    # read the necessary data, if it's not been read already
    $self->_read_pfam_data;
  }

  $self->{_max_seqname} = 0;

  # if there's nothing in "_sequence" try to load a fasta file
  $self->_read_fasta
    unless $self->{_sequence};

  # check again for a sequence. If we don't have one at this point, bail with
  # an error
  croak qq(FATAL: no sequence given)
    unless $self->{_sequence};

  # read fasta file, store maximum sequence name and store sequences for active
  # sites prediction
  $self->_parse_sequence
    unless $self->{_max_seqname};

  if ( $self->{_as} ) {
    $self->_parse_act_site_data
      unless $self->{_read_read_act_site_data};
  }

  if ( $self->{_translate} ) {
    $self->_translate_fasta;
  }

  # see if a version number was specified
  $self->{_version} = $args->{version};

}

#-------------------------------------------------------------------------------

=head2 _build_header

Adds version to the header object

=cut

sub _build_header {
  my ( $self, $version ) = @_;

  unshift @{ $self->{_header} },
    '#      query sequence file: ' . $self->{_fasta} . "\n";

  unshift @{ $self->{_header} }, <<EOF_license;
# Copyright (c) 2009 Genome Research Ltd
# Freely distributed under the GNU 
# General Public License
#
# Authors: Jaina Mistry (jm14\@sanger.ac.uk), John Tate (jt6\@sanger.ac.uk), 
#          Rob Finn (rdf\@sanger.ac.uk)
#
# This is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>. 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
EOF_license

  my $v =
    ( defined $self->{_version} )
    ? "version $version, "
    : '';

  unshift @{ $self->{_header} },
    "# pfam_scan.pl, $v run at " . scalar(localtime) . "\n#\n";
    
    
}

#-------------------------------------------------------------------------------

=head2 _read_fasta

Reads a sequence from the fasta-format file that was specified in the
parameters.

=cut

sub _read_fasta {
  my $self = shift;

  open( FASTA, $self->{_fasta} )
    or croak qq(FATAL: Couldn't open fasta file "$self->{_fasta}" $!\n);
  my @rows = <FASTA>;
  close FASTA;

  $self->{_sequence_rows} = \@rows;

  $self->{_sequence} = join '', @rows;
}

#-------------------------------------------------------------------------------

=head2 _resolve_clan_overlap

Resolves overlaps between clans.

=cut

sub _resolve_clan_overlap {
  my $self = shift;

  my @no_clan_overlap = ();
  foreach my $result ( @{ $self->{_all_results} } ) {
    my $new =
      $result->remove_overlaps_by_clan( $self->{_clanmap}, $self->{_nested} );

    push @no_clan_overlap, $new;
  }

  $self->{_all_results} = \@no_clan_overlap;
}

#-------------------------------------------------------------------------------

=head2 _pred_act_sites

Predicts active sites. Takes no arguments. Populates the "act_site" field on
each results object.

=cut

sub _pred_act_sites {
  my $self = shift;

  # print STDERR "predicting active sites...\n";

  my $hmm_file = $self->{_dir} . '/Pfam-A.hmm';

RESULT: foreach my $result ( @{ $self->{_all_results} } ) {

    # print STDERR "result: |" . $result->seqName . "|\n";

  UNIT: foreach my $unit ( @{ $result->units } ) {

      # print STDERR "family: |" . $unit->name . "|\n";

      next UNIT
        unless ( $self->{_act_site_data}->{ $unit->name }->{'alignment'} );

      my $seq_region = substr(
        $self->{_seq_hash}->{ $result->seqName },
        $unit->seqFrom - 1,
        $unit->seqTo - $unit->seqFrom + 1
      );

      my $seq_se = $unit->seqFrom . '-' . $unit->seqTo;

      # print STDERR "seq_id:     |" . $result->seqName . "|\n";
      # print STDERR "seq_se:     |" . $seq_se . "|\n";
      # print STDERR "seq_region: |" . $seq_region . "|\n";
      # print STDERR "family:     |" . $unit->name . "|\n";
      # print STDERR "hmm_file:   |" . $hmm_file . "|\n";
      # print STDERR "dir:        |" . $self->{_dir} . "|\n";

      $unit->{act_site} = Bio::Pfam::Active_site::as_search::find_as(
        $self->{_act_site_data}->{ $unit->name }->{'alignment'},
        $self->{_act_site_data}->{ $unit->name }->{'residues'},
        $result->seqName,
        $seq_se,
        $seq_region,
        $unit->name,
        $hmm_file
      );
    }
  }
}

#-------------------------------------------------------------------------------

=head2 _read_pfam_data

Reads the Pfam data file ("Pfam-A.scan.dat") and populates the C<accmap>,
C<nested> and C<clanmap> hashes on the object.

=cut

sub _read_pfam_data {
  my $self = shift;

  #print STDERR "reading " . $self->{_hmmlib} . ".dat\n" if($ENV{DEBUG});
  $self->{_accmap}    = {};
  $self->{_nested}    = {};
  $self->{_clanmap}   = {};
  $self->{_desc}      = {};
  $self->{_seqGA}     = {};
  $self->{_domGA}     = {};
  $self->{_type}      = {};
  $self->{_model_len} = {};

  foreach my $hmmlib ( @{ $self->{_hmmlib} } ) {
    my $scandat = $self->{_dir} . '/' . $hmmlib . '.dat';
    open( SCANDAT, $scandat )
      or croak qq(FATAL: Couldn't open "$scandat" data file: $!);
    my $id;
    while (<SCANDAT>) {
      if (m/^\#=GF ID\s+(\S+)/) {
        $id = $1;
      }
      elsif (m/^\#=GF\s+AC\s+(\S+)/) {
        $self->{_accmap}->{$id} = $1;
      }
      elsif (m/^\#=GF\s+DE\s+(.+)/) {
        $self->{_desc}->{$id} = $1;
      }
      elsif (m/^\#=GF\s+GA\s+(\S+)\;\s+(\S+)\;/) {
        $self->{_seqGA}->{$id} = $1;
        $self->{_domGA}->{$id} = $2;
      }
      elsif (m/^\#=GF\s+TP\s+(\S+)/) {
        $self->{_type}->{$id} = $1;
      }
      elsif (m/^\#=GF\s+ML\s+(\d+)/) {
        $self->{_model_len}->{$id} = $1;
      }
      elsif (/^\#=GF\s+NE\s+(\S+)/) {
        $self->{_nested}->{$id}->{$1} = 1;
        $self->{_nested}->{$1}->{$id} = 1;
      }
      elsif (/^\#=GF\s+CL\s+(\S+)/) {
        $self->{_clanmap}->{$id} = $1;
      }
    }

    close SCANDAT;

    # set a flag to show that we've read the data files already
    $self->{ '_read_' . $hmmlib } = 1;
  }

}

#-------------------------------------------------------------------------------

=head2 _convert_results_search_to_scan

Converts the search format to the scan format

=cut

sub _convert_results_search_to_scan {
  my $self  = shift;
  my $search_results =  shift;

  my $scan_results = {};

  foreach my $search_result (@{$search_results}) {

    foreach my $seq_id (keys %{$search_result->seqs}) {

      my $this_seq_obj = $search_result->seqs->{$seq_id};

      if (! defined($scan_results->{$seq_id})) {
        my $new_scan_result = Bio::Pfam::HMM::HMMResults->new;
        $new_scan_result->seqName($this_seq_obj->name);
        $new_scan_result->description($this_seq_obj->desc);
        $new_scan_result->program($search_result->program);

        $scan_results->{$seq_id} = $new_scan_result;
      }

      my $this_scan_result = $scan_results->{$seq_id};

      $this_scan_result->addHMMSeq(
        Bio::Pfam::HMM::HMMSequence->new({
          evalue     => $this_seq_obj->evalue,
          bits       => $this_seq_obj->bits,
          bias       => $this_seq_obj->bias,
          exp        => $this_seq_obj->exp,
          numberHits => $this_seq_obj->numberHits,
          name       => $search_result->seedName,
          desc       => $search_result->description
        })
      );

      foreach my $search_unit (@{$this_seq_obj->hmmUnits}) {
        $this_scan_result->addHMMUnit(
          Bio::Pfam::HMM::HMMUnit->new({
            name      => $search_result->seedName,
            domain    => $search_unit->domain,
            hmmalign  => $search_unit->hmmalign,
            bits      => $search_unit->bits,
            bias      => $search_unit->bias,
            domEvalue => $search_unit->domEvalue,
            evalue    => $search_unit->evalue,
            hmmFrom   => $search_unit->hmmFrom,
            hmmTo     => $search_unit->hmmTo,
            seqFrom   => $search_unit->seqFrom,
            seqTo     => $search_unit->seqTo,
            envFrom   => $search_unit->envFrom,
            envTo     => $search_unit->envTo,
            aliAcc    => $search_unit->aliAcc
          })
        );
      }

      $this_scan_result->eof(1);

    }
  }

  my @ordered_keys = sort {$a cmp $b} keys(%{$scan_results});
  my @values = @{$scan_results}{@ordered_keys};

  return \@values;
}

#-------------------------------------------------------------------------------

=head2 _read_act_site_data

Reads the Pfam active site data file ("active_site.dat") and populates
the C<act_site_data> hashes on the object.

=cut

sub _parse_act_site_data {
  my $self   = shift;
  my $as_dat = $self->{_dir} . '/active_site.dat';

  $self->{_act_site_data} = {};

  open( AS, $as_dat )
    or croak qq(FATAL: Couldn\'t open "$as_dat" data file: $!);

  my ( $fam_id, $aln );

  while (<AS>) {
    if (/^ID\s+(\S+)/) {
      $fam_id = $1;
      $aln    = new Bio::SimpleAlign;
    }
    elsif (/^AL\s+(\S+)\/(\d+)\-(\d+)\s+(\S+)/) {
      my ( $seq_id, $st, $en, $seq ) = ( $1, $2, $3, $4 );

      $aln->add_seq(
        Bio::Pfam::Scan::Seq->new(
          '-seq'   => $seq,
          '-id'    => $seq_id,
          '-start' => $st,
          '-end'   => $en,
          '-type'  => 'aligned'
        )
      );
    }
    elsif (/^RE\s+(\S+)\s+(\d+)/) {
      my ( $seq_id, $res ) = ( $1, $2 );
      push(
        @{ $self->{_act_site_data}->{$fam_id}->{'residues'}->{$seq_id} },
        $res
      );

    }
    elsif (/^\/\//) {

      $self->{_act_site_data}->{$fam_id}->{'alignment'} = $aln;

      $fam_id = "";
      $aln    = "";

    }
    else {
      warn "Ignoring line:\n[$_]";
    }
  }
  close AS;
  $self->{_read_read_act_site_data} = 1;
}

#-------------------------------------------------------------------------------

=head2 _parse_sequence

This method is used to parse the sequence and hash it on sequence
identifier. It also stores the length of the longest sequence id

=cut

sub _parse_sequence {
  my $self = shift;

  my $seq_hash = {};
  my $seq_id;
  foreach ( @{ $self->{_sequence_rows} } ) {

    next if m/^\s*$/;    #Ignore blank lines

    if (m/^>(\S+)/) {
      $seq_id = $1;

      if ( exists( $seq_hash->{$seq_id} ) ) {
        croak "FATAL: Sequence identifiers must be unique. Your fasta file contains two sequences with the same id ($seq_id)";
      }

      #Store the max length of seq name, use this later when printing in ascii
      $self->{_max_seqname} = length($seq_id)
        if ( !$self->{_max_seqname}
        or length($seq_id) > $self->{_max_seqname} );
    }
    else {
      croak "FATAL: Unrecognised format of fasta file. Each sequence must have a header line in the format '>identifier  <optional description>'"
        unless $seq_id;
      chomp;
      $seq_hash->{$seq_id} .= $_;
    }
  }

  $self->{_seq_hash} = $seq_hash;
}

#-------------------------------------------------------------------------------

=head2 _translate_fasta

Uses the HMMER v2.3.2 progam "translate" to perform a six-frame translation of
the input sequence. Checks the parameter "-translate".

Accepted arguments are "all" and "orf", where "all" means (from the "translate"
help text) "translate in full, with stops; no individual ORFs" and "orf" means
"report only ORFs greater than minlen" where minlen is set to the default of
20.

=cut

sub _translate_fasta {
  my ($self) = @_;
  my $translatedFasta = $self->{_fasta} . ".translated";

  my @params = ( 'translate', '-q', );
  if ( $self->{_translate} eq 'all' ) {
    push( @params, '-a' );
  }
  elsif ( $self->{_translate} eq 'orf' ) {
    push( @params, '-l', '20' );
  }
  else {
    croak qq(Unexpected parameter '$self->{_translate}');
  }
  push( @params, '-o', $translatedFasta, $self->{_fasta} );

  print STDERR "PfamScan::translate_fasta: translate command: |@params|\n"
    if $ENV{DEBUG};

  my $run = start \@params, '<pipe', \*IN, '>pipe', \*OUT, '2>pipe', \*ERR
    or croak qq(FATAL: error running translate; IPC::Run returned '$?');

  close IN;
  close OUT;

  my $err;
  while (<ERR>) {
    $err .= $_;
  }
  close ERR;

  finish $run
    or croak qq|FATAL: error running translate ($err); ipc returned '$?'|;
  open( F, "<", $translatedFasta )
    or croak qw(Could not open $translatedFasta '$!');
  if ( $self->{_translate} eq 'orf' ) {
    while (<F>) {
      if (/^>\s?(\S+).*nt (\d+)\.+(\d+)/) {
        $self->{_orf}->{$1}->{start}  = $2;
        $self->{_orf}->{$1}->{end}    = $3;
        $self->{_orf}->{$1}->{strand} = ( $2 < $3 ) ? '+' : '-';
      }
    }
  }
  else {
    my $currentSeq;
    my $currentFrame;
    my $currentLen = 0;
    my $maxEnd = 0;
    while (<F>) {
      chomp;
      if (/^>\s?(\S+\:)(\d+)/) {
        if ( $currentLen > 0 ) {
          my $seqName = $currentSeq . $currentFrame;
          if ( $currentFrame < 3 ) {
            my $start = 1 + $currentFrame;
            my $end   = $start + $currentLen - 1;
            $self->{_orf}->{$seqName}->{strand} = '+';
            $self->{_orf}->{$seqName}->{start}  = $start;
            $self->{_orf}->{$seqName}->{end}    = $end;
            $maxEnd = $end if ( $end > $maxEnd );
          }
          else {
            my $start = $maxEnd - ( $currentFrame - 3 );
            my $end = $start - $currentLen + 1;
            $self->{_orf}->{$seqName}->{strand} = '-';
            $self->{_orf}->{$seqName}->{start}  = $start;
            $self->{_orf}->{$seqName}->{end}    = $end;
          }
        }
        $currentLen   = 0;
        $currentSeq   = $1;
        $currentFrame = $2;
      }
      else {
        $currentLen += length($_) * 3;
      }
    }
    my $seqName = $currentSeq . $currentFrame;
    if ( $currentFrame < 3 ) {
      my $start = 1 + $currentFrame;
      my $end   = $start + $currentLen - 1;
      $self->{_orf}->{$seqName}->{strand} = '+';
      $self->{_orf}->{$seqName}->{start}  = $start;
      $self->{_orf}->{$seqName}->{end}    = $end;
      $maxEnd = $end if ( $end > $maxEnd );
    }
    else {
      my $start = $maxEnd - ( $currentFrame - 3 );
      my $end = $start - $currentLen + 1;
      $self->{_orf}->{$seqName}->{strand} = '-';
      $self->{_orf}->{$seqName}->{start}  = $start;
      $self->{_orf}->{$seqName}->{end}    = $end;
    }
  }
  $self->{_fasta} = $translatedFasta;
}
#-------------------------------------------------------------------------------

=head1 COPYRIGHT

Copyright (c) 2009: Genome Research Ltd.

Authors: Jaina Mistry (jm14@sanger.ac.uk), John Tate (jt6@sanger.ac.uk), Rob Finn (finnr@janelia.hhmi.org)

This is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
or see the on-line version at http://www.gnu.org/copyleft/gpl.txt

=cut

  1;

