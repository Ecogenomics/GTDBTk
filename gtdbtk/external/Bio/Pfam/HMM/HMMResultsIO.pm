# HMMResultsIO.pm
#
# Author:        rdf
# Maintainer:    $Id: HMMResultsIO.pm,v 1.2 2009-12-01 15:42:20 jt6 Exp $
# Version:       $Revision: 1.2 $
# Created:       Nov 16, 2008
# Last Modified: $Date: 2009-12-01 15:42:20 $

=head1 NAME

Template - a short description of the class

=cut

package Bio::Pfam::HMM::HMMResultsIO;

=head1 DESCRIPTION

A more detailed description of what this class does and how it does it.

$Id: HMMResultsIO.pm,v 1.2 2009-12-01 15:42:20 jt6 Exp $

=head1 COPYRIGHT

File: HMMResultsIO.pm

Copyright (c) 2007: Genome Research Ltd.

Authors: Rob Finn (rdf@sanger.ac.uk), John Tate (jt6@sanger.ac.uk)

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

use strict;
use warnings;
use Moose;
use Carp;

#All the things we need to objectfy the search results
use Bio::Pfam::HMM::HMMResults;
use Bio::Pfam::HMM::HMMSequence;
use Bio::Pfam::HMM::HMMUnit;

#-------------------------------------------------------------------------------

=head1 ATTRIBUTES



=cut

has 'align' => (
  isa     => 'Int',
  is      => 'rw',
  default => 0
);

has 'outfile' => (
  isa     => 'Str',
  is      => 'rw',
  default => 'OUTPUT'
);

has 'pfamout' => (
  isa     => 'Str',
  is      => 'rw',
  default => 'PFAMOUT'
);

has 'scores' => (
  isa     => 'Str',
  is      => 'rw',
  default => 'scores'
);

#-------------------------------------------------------------------------------

=head1 METHODS

=head2 parseHMMER3

  Title    : parseHMMER 
  Usage    : $hmmResIO->parseHMMSearch( filename )
  Function : Parse the output from a HMMER3 search results 
  Args     : Filename containing the search 
  Returns  : A Bio::Pfam::HMM::HMMResults object
  
=cut

sub parseHMMER3 {
  my ( $self, $filename ) = @_;
  my $fh;
 
  if(ref($filename) eq 'GLOB'){
    $fh = $filename;
  }else{
    open( $fh, $filename ) or confess "Could not open $filename:[$!]\n";
  }
  
#  open( $fh, $filename ) or confess "Could not open $filename:[$!]\n";
  my $hmmRes = Bio::Pfam::HMM::HMMResults->new;
  $self->_readHeader( $fh, $hmmRes );
  $self->_readSeqHits( $fh, $hmmRes );
  $self->_readUnitHits( $fh, $hmmRes );
  $self->_readFooter($fh, $hmmRes);
  return ($hmmRes);
}



sub parseMultiHMMER3 {
  my ( $self, $filename ) = @_;
  my $fh;
  
  if(ref($filename) eq 'GLOB'){
    $fh = $filename;
  }elsif( ref($filename) and $filename->isa('IO::File') ) {
      $fh = $filename;
  }else{
    open( $fh, $filename ) or confess "Could not open $filename:[$!]\n";
  }  
  
  my @hmmResAll;
  my $program;
  while(!eof($fh)){
    my $hmmRes = Bio::Pfam::HMM::HMMResults->new;
    my $eof = $self->_readHeader( $fh, $hmmRes ); 
    last if($eof);
    push(@hmmResAll, $hmmRes);
    if($hmmRes->program) {
	$program = $hmmRes->program;
    }
    else {
	$hmmRes->program($program);
    }
    $self->_readSeqHits( $fh, $hmmRes );
    $self->_readUnitHits( $fh, $hmmRes );
    $self->_readFooter($fh, $hmmRes);
  }
  return (\@hmmResAll);
}

sub parseSplitHMMER3 {
  my($self, $files ) = @_;
  
  my $hmmRes = Bio::Pfam::HMM::HMMResults->new;
  
  foreach my $filename (@{$files}){
    my ($fh);
    open( $fh, $filename ) or confess "Could not open $filename:[$!]\n";
    $self->_readHeader( $fh, $hmmRes );
    $self->_readSeqHits( $fh, $hmmRes );
    $self->_readUnitHits( $fh, $hmmRes );
    $self->_readFooter($fh, $hmmRes);
  }
  
  return ( $hmmRes );
    
}


#-------------------------------------------------------------------------------

=head2 convertHMMSearch 

  Title    : convertHMMSearch
  Usage    : $hmmResIO->convertHMMSearch('SEARCHFILE') 
  Function : This wraps up a couple of methods to convert the more complex hmmsearch 
           : results in to nice clean format that we Pfam-ers are used to. 
  Args     : The filename of the hmmsearch output file
  Returns  : Nothing
  
=cut

sub convertHMMSearch {
  my ( $self, $filename ) = @_;

  unless ($filename) {
    confess "No filename passed in to convertHMMSearch\n";
  }
  unless ( -s $filename ) {
    confess "$filename does not exists\n";
  }

  #Now parse in the raw HMM output and write out the results as a PFAMOUT
  my $hmmRes = $self->parseHMMER3($filename);
  $self->writePFAMOUT($hmmRes);
  return $hmmRes;
}

#-------------------------------------------------------------------------------

=head2 writePFAMOUT 

  Title    : writePFAMOUT
  Usage    : $hmmResIO->writePFAMOUT( $hmmRes ) 
  Function : Writes a Bio::Pfam::HMM:HMMResults object in to a PFAMOUT file. 
  Args     : A Bio::Pfam::HMM:HMMResults
  Returns  : Nothing
  
=cut

sub writePFAMOUT {
  my ( $self, $hmmRes ) = @_;

  unless ($hmmRes) {
    confess "A Bio::Pfam::HMM::HMMResults object was not parsed in\n";
  }
  unless ( $hmmRes->isa("Bio::Pfam::HMM::HMMResults") ) {
    confess("Variable passed in is not a Bio::Pfam::HMM::Results object");
  }

  my $fh;
  open( $fh, ">" . $self->pfamout )
    or confess "Could not open " . $self->pfamout . ":[$!]\n";

  print $fh <<HEAD;
# ===========
# Pfam output
# ===========
#
# Sequence scores
# ---------------
#
# name      description                                   bits      evalue   n   exp  bias  
 
HEAD

  foreach
    my $seq ( sort { $b->bits <=> $a->bits } ( @{ $hmmRes->eachHMMSeq } ) )
  {
    $_ = $seq->desc;
    my ($desc) = /^(.{1,42})/;
    $desc = uc($desc);
    printf $fh (
      "%-10s  %-42s %8.1f  %9s %3d %5.1f %5.1f\n",
      $seq->name,
      $desc,
      $seq->bits,
      $seq->evalue,
      scalar( @{ $seq->hmmUnits } ),
      defined( $seq->exp )  ? $seq->exp  : "-",
      defined( $seq->bias ) ? $seq->bias : "-"
    );
  }

  print $fh <<HEAD;
#
# Domain scores
# -------------
#
# name      env-st  env-en  ali-st  ali-en  hmm-st  hmm-en   bits      evalue    hit   bias
#

HEAD

  foreach my $dom ( sort { $b->bits <=> $a->bits } @{ $hmmRes->units } ) {
    
    
    printf $fh (
      "%-10s  %6d  %6d  %6d  %6d  %6s  %6s %6.1f  %9s %6d %6.1f\n",
      $dom->name,
      $dom->envFrom,
      $dom->envTo,
      $dom->seqFrom,
      $dom->seqTo,
      $dom->hmmFrom,
      $dom->hmmTo,
      $dom->bits,
      $dom->evalue,
      $dom->domain,
      defined( $dom->bias ) ? $dom->bias : "-",

    );
  }
}

#-------------------------------------------------------------------------------

=head2 parsePFAMOUT 

  Title    : parsePFAMOUT
  Usage    : $self->parsePFAMOUT($filename)
  Function : Reads in a PFAMOUT file.  This file contains the minimal amount of information
           : require to constrcut a pfam ALIGN file.   
  Args     : A filename. Normally this is filename
  Returns  : A Bio::Pfam::HMM::HMMResults object
  
=cut

sub parsePFAMOUT {
  my $self     = shift;
  my $filename = shift;

  unless ($filename) {
    confess('No filename or filehandle passed to parsePFAMOUT');
  }

  my $fh;
  if ( ref($filename) eq 'GLOB' ) {
    $fh = $filename;
  }
  else {
    open( $fh, $filename ) or confess "Could not open $filename:[$!]\n";
  }
  my $hmmRes = Bio::Pfam::HMM::HMMResults->new;

  while (<$fh>) {
    /^# Domain scores/ && last;

    #if (/^(\S+)\s+(.*?)\s+(\S+)\s+(\S+)\s+(\d+)\s*$/) {
    if (/^(\S+)\s+(.*?)\s+(\S+)\s+(\S+)\s+(\d+)\s+\S+\s+(\S+)\s*$/) {
      
      $hmmRes->addHMMSeq(
        Bio::Pfam::HMM::HMMSequence->new(
          {
            name       => $1,
            desc       => $2,
            bits       => $3,
            evalue     => $4,
            numberHits => $5,
	    bias       => $6
          }
        )
      );
    }
    elsif (/^#|^\s+$/) {
      next;
    }
    else {
      warn "Did not parse|$_|\n";
    }
  }
  while (<$fh>) {

    #if (/^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s*$/) {
    if (
      /^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)/)
    {
      $hmmRes->addHMMUnit(
        Bio::Pfam::HMM::HMMUnit->new(
          {
            name     => $1,
            envFrom  => $2,
            envTo    => $3,
            seqFrom  => $4,
            seqTo    => $5,
            hmmFrom  => $6,
            hmmTo    => $7,
            bits     => $8,
            evalue   => $9,
	    bias     => $10
          }
        )
      );
    }
    elsif (/^#|^\s+$/) {
      next;
    }
    elsif (/^$/) {
      next;
    }
    else {
      warn "Did not parse: |$_|";
    }
  }
  close($fh);
  return ($hmmRes);
}

#-------------------------------------------------------------------------------

=head2 _readHeader 

  Title    : _readHeader
  Usage    : Private method.  $self->_readHeader(\*FH, $hmmResults)
  Function : Reads the header section from a HMMER3 hmmsearch   
  Args     : The file handle to hmmsearch output, a Bio::Pfam::HMM::HMMResults object
  Returns  : Nothing
  
=cut

#Parse the header part of the output first;
sub _readHeader {
  my ( $self, $hs, $hmmRes ) = @_;

  #Check the $hs is defined and a GLOB

  while (<$hs>) {
    if (/^Scores for complete/) {
      last;
    }
    elsif (/^# query HMM file:\s+(\S+)/) {
      $hmmRes->hmmName($1);
    }
    elsif (/^# target sequence database:\s+(\S+)/) {
      $hmmRes->seqDB($1);
    }
    elsif (/^output directed to file:\s+(\S+)/) {
      $hmmRes->thisFile($1);
    }
    elsif (/^Query:\s+(\S+)\s+\[M\=(\d+)\]/) {
      $hmmRes->seedName($1);
      $hmmRes->hmmLength($2);
    }elsif(/^Query:\s+(\S+)\s+\[L\=(\d+)\]/) {
	$hmmRes->seqName($1); 
      $hmmRes->seqLength($2);
    }elsif (/^sequence E-value threshold: <= (\d+)/) {
      $hmmRes->evalueThr($1);
    }
    elsif (/^# Random generator seed:      (\d+)/) {
      $hmmRes->randSeedNum($1);
    }elsif(/^Description:\s+(.*)/){
      $hmmRes->description($1);  
    }elsif(/^# (phmmer|hmmsearch|hmmscan|jackhmmer)/){
      $hmmRes->program($1);
    }elsif (/(^#)|(^$)/) {
      next;
    }elsif(/^Accession/){
      next; 
    } elsif(/^\[ok\]/) {
      return(1);
    } else {
      die "Failed to parse hmmsearch results |$_| in header section\n";
    }
  }
}

#-------------------------------------------------------------------------------

=head2 _readSeqHits 

  Title    : _readSeqHits
  Usage    : Private method.  $self->_readSeqHits(\*FH, $hmmResults)
  Function : Reads the sequence hits from a HMMER3 hmmsearch   
  Args     : The file handle to hmmsearch output, a Bio::Pfam::HMM::HMMResults object
  Returns  : Nothing
  
=cut

sub _readSeqHits {
  my ( $self, $hs, $hmmRes ) = @_;
  while (<$hs>) {

#Match a line like this
# E-value  score  bias    E-value  score  bias    exp  N  Sequence Description
#    ------- ------ -----    ------- ------ -----   ---- --  -------- -----------
#      4e-83  285.8  10.0    5.3e-83  285.5   7.0    1.1  1  Q14SN3.1 Q14SN3_9HEPC Polyprotein (Fragment).
    if (/^Domain annotation for each [sequence|model]/) { # This is the format for HMMER3b3
      last;
    }
    elsif (/^Domain and alignment annotation for each [sequence|model]/) {  #This is the format for HMMER3b2 - can be removed later
	last;
    }
    elsif (/^\s+(E-value|---)/) {
      next;
    }
    elsif (/^$/) {
      next;
    }
    else {
      next if(/No hits detected that satisfy reporting thresholds/);

      #Assume that we have a sequence match
      my @sMatch = split( /\s+/, $_ );
      unless ( scalar(@sMatch) >= 10 ) {
        die "Expected at least 10 pieces of data: $_;\n";
      }
      my $desc;
      if ( scalar(@sMatch) >= 11 ) {
        $desc = join( " ", @sMatch[ 10 .. $#sMatch ] );
      }

      $hmmRes->addHMMSeq(
        Bio::Pfam::HMM::HMMSequence->new(
          {
            evalue     => $sMatch[1],
            bits       => $sMatch[2],
            bias       => $sMatch[3],
            exp        => $sMatch[7],
            numberHits => $sMatch[8],
            name       => $sMatch[9],
            desc       => defined($desc) ? $desc : "-",
          }
        )
      );
      
      next;
    }
    die "Failed to parse $_ in sequence section\n";
  }

}

#------------------------------------------------------------------------------

=head2 _readUnitHits 

  Title    : _readUnitHits
  Usage    : Private method.  $self->_readUnitHits(\*FH, $hmmResults)
  Function : Reads the unit (domain) hits from a HMMER3 hmmsearch   
  Args     : The file handle to hmmsearch output, a Bio::Pfam::HMM::HMMResults object
  Returns  : Nothing
  
=cut

no warnings 'recursion';

sub _readUnitHits {
  my ( $self, $hs, $hmmRes ) = @_;
  
  if($hmmRes->eof){
    return; 
  }

#Parse the domain hits section
#>> P37935.1  MAAY4_SCHCO Mating-type protein A-alpha Y4.
#     # bit score    bias    E-value ind Evalue hmm from   hmm to    ali from   ali to    env from   env to    ali-acc
#   --- --------- ------- ---------- ---------- -------- --------    -------- --------    -------- --------    -------
#     1     244.0     0.5    9.5e-76    1.7e-70        1      146 [.        1      145 [.        1      146 [.    0.99
#
#  Alignments for each domain:
#  == domain 1    score: 244.0 bits;  conditional E-value: 9.5e-76
#      SEED   1 medrlallkaisasakdlvalaasrGaksipspvkttavkfdplptPdldalrtrlkeaklPakaiksalsayekaCarWrsdleeafdktaksvsPanlhllealrirlyteqvekWlvqvlevaerWkaemekqrahiaatmgp 146
#               m+++la+l++isa+akd++ala+srGa+++ +p++tt+++fd+l++P+ld++rtrl+ea+lP+kaik++lsaye+aCarW++dleeafd+ta+s+sP+n+++l++lr+rly+eqv+kWl++vl+v+erWkaemekqrahi+atmgp
#  P37935.1   1 MAELLACLQSISAHAKDMMALARSRGATGS-RPTPTTLPHFDELLPPNLDFVRTRLQEARLPPKAIKGTLSAYESACARWKHDLEEAFDRTAHSISPHNFQRLAQLRTRLYVEQVQKWLYEVLQVPERWKAEMEKQRAHINATMGP 145
#               899***************************.******************************************************************************************************************8 PP

  while (<$hs>) {
    if (/^Internal/) {
      last;
    }
    elsif (/\>\>\s+(\S+)/) {
      my $seqId = $1;    
      $self->_readUnitData( $seqId, $hs, $hmmRes );
      if($hmmRes->eof){
        return; 
      }
    }
  }
}

sub _readUnitData {
  my ( $self, $id, $hs, $hmmRes ) = @_;
  
  if($hmmRes->eof){
    return; 
  }
  my $hmmName = $hmmRes->seedName();

  my $seqName = $hmmRes->seqName;

# bit score    bias    E-value ind Evalue hmm from   hmm to    ali from   ali to    env from   env to    ali-acc
#   --- --------- ------- ---------- ---------- -------- --------    -------- --------    -------- --------    -------
#     1     244.0     0.5    9.5e-76    1.7e-70        1      146 [.        1      145 [.        1      146 [.    0.99
#
#  Alignments for each domain:

  my @units;
  my $align   = 1;
  my $recurse = 0;
  my $eof = 0;
  my ($nextSeqId);
  while (<$hs>) {
    if (/^[(\/\/|Internal)]/ ) {
      $align   = 0;
      $recurse = 0;
      $eof = 1;
      last;
    }
    elsif (/^\>\>\s+(\S+)/) {
      $nextSeqId = $1;
      $align     = 0;
      $recurse   = 1;
      last;
    }
    elsif (/^\s+Alignments for each domain:/) {
      $align   = 1;
      $recurse = 0;
      last;
    }
    elsif (/^\s+(#\s+score|---)/){

      #Two human readable lines
      next;
    }
    elsif (/^$/) {

      #blank line
      next;
    }
    elsif (/^\s+\d+\s+/) {
      my @dMatch = split( /\s+/, $_ );
      unless ( scalar(@dMatch) == 17 ) {
        die "Expected 16 elements of data: $_\n";
      }

      push(
        @units,
        Bio::Pfam::HMM::HMMUnit->new(
          {
            name      => $id,
            domain    => $dMatch[1],
            bits      => $dMatch[3],
            bias      => $dMatch[4],
            domEvalue => $dMatch[5],
            evalue    => $dMatch[6],
            hmmFrom   => $dMatch[7],
            hmmTo     => $dMatch[8],
            seqFrom   => $dMatch[10],
            seqTo     => $dMatch[11],
            envFrom   => $dMatch[13],
            envTo     => $dMatch[14],
            aliAcc    => $dMatch[16]
          }
        )
      );
      
      next;
    }
    elsif(/^\s+\[No individual domains/) {
	$align=0;
	next;
    }
    else {
      confess("Did not parse line: $_");
    }
  }

#  == domain 1    score: 244.0 bits;  conditional E-value: 9.5e-76
#      SEED   1 medrlallkaisasakdlvalaasrGaksipspvkttavkfdplptPdldalrtrlkeaklPakaiksalsayekaCarWrsdleeafdktaksvsPanlhllealrirlyteqvekWlvqvlevaerWkaemekqrahiaatmgp 146
#               m+++la+l++isa+akd++ala+srGa+++ +p++tt+++fd+l++P+ld++rtrl+ea+lP+kaik++lsaye+aCarW++dleeafd+ta+s+sP+n+++l++lr+rly+eqv+kWl++vl+v+erWkaemekqrahi+atmgp
#  P37935.1   1 MAELLACLQSISAHAKDMMALARSRGATGS-RPTPTTLPHFDELLPPNLDFVRTRLQEARLPPKAIKGTLSAYESACARWKHDLEEAFDRTAHSISPHNFQRLAQLRTRLYVEQVQKWLYEVLQVPERWKAEMEKQRAHINATMGP 145
#               899***************************.******************************************************************************************************************8 PP
#
# OR....
#
#  == domain 1    score: 27.6 bits;  conditional E-value: 7.4e-10
#   PF00018  17 LsfkkGdvitvleksee.eWwkaelkdg.keGlvPsnYvep 55 
#               L++++Gd+++++++++e++Ww++++++++++G++P+n+v+p
#  P15498.4 617 LRLNPGDIVELTKAEAEqNWWEGRNTSTnEIGWFPCNRVKP 657
#               7899**********9999*******************9987 PP


  if ($align) {
    my ($pattern1, $pattern2);

    if($hmmName and $hmmRes->program eq 'hmmsearch'){
      $pattern1 = qr/^\s+$hmmName\s+\d+\s+(\S+)\s+\d+/;
      $id =~ s/(\W)/\\$1/g; # escape any non-word character
      # $id =~ s/\|/\\|/g;  #Escape '|', '[' and ']' characters
      # $id =~ s/\[/\\[/g;
      # $id =~ s/\]/\\]/g;
      $pattern2 = qr/^\s+$id\s+\d+\s+(\S+)\s+\d+/;
    }elsif($seqName and $hmmRes->program eq 'hmmscan'){
      my $tmpSeqName = $seqName;
      $tmpSeqName =~ s/(\W)/\\$1/g; # escape any non-word character
      # $tmpSeqName =~ s/\|/\\|/g; #Escape '|', '[' and ']' characters
      # $tmpSeqName =~ s/\[/\\[/g;  
      # $tmpSeqName =~ s/\]/\\]/g;
      $pattern1 = qr/^\s+$id\s+\d+\s+(\S+)\s+\d+/;
      $pattern2 = qr/^\s+$tmpSeqName\s+\d+\s+(\S+)\s+\d+/;
    }elsif($seqName and ($hmmRes->program eq 'phmmer' or $hmmRes->program eq 'jackhmmer') ){
      $seqName =~ s/(\W)/\\$1/g; # escape any non-word character
      # $seqName =~ s/\|/\|/g; #Escape '|', '[' and ']' characters
      # $seqName =~ s/\[/\\[/g;
      # $seqName =~ s/\]/\\]/g;
      $pattern1 = qr/^\s+$seqName\s+\d+\s+(\S+)\s+\d+/;
      $pattern2 = qr/^\s+$id\s+\d+\s+(\S+)\s+\d+/;
    }
  

    $recurse = 0;
    my $matchNo;
    my $hmmlen = 0;
    while (<$hs>) {
      if (/$pattern1/) {
        $units[ $matchNo - 1 ]->hmmalign->{hmm} .= $1;
        $hmmlen = length($1);
      }
      elsif (/$pattern2/) {
        $units[ $matchNo - 1 ]->hmmalign->{seq} .= $1;
      }
      elsif (/^\s+([x\.]+)\s+RF$/) {
        my $rf = $1;
        $units[ $matchNo - 1 ]->hmmalign->{rf} .= $rf;
      }
      elsif (/^\s+([0-9\*\.]+)\s+PP$/) {
        my $pp = $1;
        $units[ $matchNo - 1 ]->hmmalign->{pp} .= $pp;
      }elsif (/^\s+(\S+)\s+CS$/) {
        my $cs = $1;
        $units[ $matchNo - 1 ]->hmmalign->{cs} .= $cs;
      }elsif (/^\s+==\s+domain\s+(\d+)/) {
        $matchNo = $1;
      }
      elsif (/^\s+(.*)\s+$/) {
      	# $1 is *not* the match - this fails if there are prepended
      	# or appended spaces
        # $units[ $matchNo - 1 ]->hmmalign->{match} .= $1;
        # Let's get a right substring based on the HMM length
        chomp;
        my $m1 = substr($_,-$hmmlen);
        $units[ $matchNo - 1 ]->hmmalign->{match} .= $m1;
      }elsif (/^$/) {
        next;
      }
      elsif (/^[(\/\/|Internal)]/) {
        $align   = 0;
        $recurse = 0;
        $eof = 1;
        last;
      }
      elsif (/^\>\>\s+(\S+)/) {
        $nextSeqId = $1;
        $recurse   = 1;
        last;
      }
     
      else {
        confess("Did not parse |$_| in units");
      }
    }
  }

  foreach my $u (@units) {
    $hmmRes->addHMMUnit($u);
  }
      
  $hmmRes->eof($eof);

  if ($recurse and $nextSeqId) {
    $self->_readUnitData( $nextSeqId, $hs, $hmmRes );
  }
  return;
}
use warnings 'recursion';

#-------------------------------------------------------------------------------

=head2 parseHMMER2 

  Title    : parseHMMER2
  Usage    : $self->parseHMMER2(\*FH )
  Function : This is a minimal parser for reading in the output of HMMER2 hmmsearch   
  Args     : The file handle to hmmsearch output
  Returns  : A Bio::Pfam::HMM::HMMResults object
  
=cut

sub parseHMMER2 {
  my $self = shift;
  my $file = shift;

  my $hmmRes = Bio::Pfam::HMM::HMMResults->new;

  my %seqh;
  my $count = 0;

  while (<$file>) {
    /^Scores for complete sequences/ && last;
  }

  while (<$file>) {
    /^Parsed for domains/ && last;
    if ( my ( $id, $de, $sc, $ev, $hits ) =
      /^(\S+)\s+(.*?)\s+(\S+)\s+(\S+)\s+(\d+)\s*$/ )
    {
      $hmmRes->addHMMSeq(
        Bio::Pfam::HMM::HMMSequence->new(
          {
            bits       => $sc,
            evalue     => $ev,
            name       => $id,
            desc       => $de,
            numberHits => $hits
          }
        )
      );

      $seqh{$id} = $sc;
    }
  }

  while (<$file>) {
    /^Histogram of all scores/ && last;
    if ( my ( $id, $sqfrom, $sqto, $hmmf, $hmmt, $sc, $ev ) =
      /^(\S+)\s+\S+\s+(\d+)\s+(\d+).+?(\d+)\s+(\d+)\s+\S+\s+(\S+)\s+(\S+)\s*$/ )
    {
      $hmmRes->addHMMUnit(
        Bio::Pfam::HMM::HMMUnit->new(
          {
            name    => $id,
            seqFrom => $sqfrom,
            seqTo   => $sqto,
            hmmFrom => $hmmf,
            hmmTo   => $hmmt,
            bits    => $sc,
            evalue  => $ev
          }
        )
      );

    }
  }

  return $hmmRes;
}

#-------------------------------------------------------------------------------

=head2 parseHMMER1

  Title    : parseHMMER1
  Usage    : $self->parseHMMER1(\*FH )
  Function : This is a minimal parser for reading in the output of HMMER1 hmmsearch.
           : There are a few hacks to get round some of them requirements     
  Args     : The file handle to hmmsearch output
  Returns  : A Bio::Pfam::HMM::HMMResults object
  
=cut

sub parseHMMER1 {
  my $self = shift;
  my $file = shift;

  my $hmmRes = Bio::Pfam::HMM::HMMResults->new;

  my %seqh;
  my $count = 0;

  while (<$file>) {
    if ( my ( $bits, $s, $e, $id, $de ) =
/^(-?\d+\.?\d*)\s+\(bits\)\s+f:\s+(\d+)\s+t:\s+(\d+)\s+Target:\s+(\S+)\s+(.*)/
      )
    {
      if ( $id =~ /(\S+)\/(\d+)-(\d+)/ ) {
        $id = $1;
        $s  = $2 + $s - 1;
        $e  = $2 + $e - 1;
      }

      if ( !$hmmRes->seqs->{$id} ) {
        $hmmRes->addHMMSeq(
          Bio::Pfam::HMM::HMMSequence->new(
            {
              bits       => $bits,
              evalue     => 1,
              name       => $id,
              desc       => $de,
              numberHits => 1
            }
          )
        );
      }
      $hmmRes->addHMMUnit(
        Bio::Pfam::HMM::HMMUnit->new(
          {
            name    => $id,
            seqFrom => $s,
            seqTo   => $e,
            hmmFrom => "1",
            hmmTo   => "1",
            bits    => $bits,
            evalue  => "1"
          }
        )
      );
      if ( $bits > $hmmRes->seqs->{$id}->bits ) {
        $hmmRes->seqs->{$id}->bits($bits);
      }
    }
  }
  return $hmmRes;
}

#-------------------------------------------------------------------------------

=head2 writeScoresFile

  Title    : writeScoresFile
  Usage    : $hmmResIO->writeScoresFile( $hmmRes)
  Function : Writes a scores file for a Bio::Pfam::HMM::HMMResults object.
  Args     : Bio::Pfam::HMM::HMMResults
  Returns  : Nothing
  
=cut

sub writeScoresFile {
  my ( $self, $hmmRes ) = @_;

  unless ($hmmRes) {
    confess "A Bio::Pfam::HMM::HMMResults object was not parsed in\n";
  }
  unless ( $hmmRes->isa("Bio::Pfam::HMM::HMMResults") ) {
    confess("Variable passed in is not a Bio::Pfam::HMM::Results object");
  }

  my $fh;
  open( $fh, ">" . $self->scores )
    or confess "Could not open " . $self->scores . ":[$!]\n";

  my ( $lowSeq, $lowDom, $highSeq, $highDom );
  $lowSeq  = $lowDom  =  999999.99;
  $highSeq = $highDom = -999999.99;
  unless ( defined $hmmRes->domThr and defined $hmmRes->seqThr ) {
    warn "No threshold set, setting to 25.0 bits\n";
    $hmmRes->domThr("25.0");
    $hmmRes->seqThr("25.0");
  }

  my @sigUnits;

  foreach my $seqId ( keys %{ $hmmRes->seqs } ) {

    #Does this sequence score above or equal to the sequence threshold?
    if ( $hmmRes->seqs->{$seqId}->bits >= $hmmRes->seqThr ) {

      #Is this the lowest sequence thresh
      if ( $hmmRes->seqs->{$seqId}->bits < $lowSeq ) {
        $lowSeq = $hmmRes->seqs->{$seqId}->bits;
      }

#For each of the regions found on the sequence, look to see if the match is great
#than the domain threshold.  If it is, is it lower than we we have seen previously
      foreach my $unit ( @{ $hmmRes->seqs->{$seqId}->hmmUnits } ) {
        if ( $unit->bits >= $hmmRes->domThr ) {
          push( @sigUnits, $unit );
          if ( $unit->bits < $lowDom ) {
            $lowDom = $unit->bits();
          }
        }
        elsif ( $unit->bits > $highDom ) {
          $highDom = $unit->bits;
        }
      }
    }
    else {

      #Is this the highest sequence thres below the cut-off
      if ( $hmmRes->seqs->{$seqId}->bits > $highSeq ) {
        $highSeq = $hmmRes->seqs->{$seqId}->bits;
      }

#For each of the regions found on the sequence, look to see if the match is great
#than the domain threshold.  If it is, is it lower than we we have seen previously
      foreach my $unit ( @{ $hmmRes->seqs->{$seqId}->hmmUnits } ) {
        if ( $unit->bits < $hmmRes->domThr && $unit->bits > $highDom ) {
          $highDom = $unit->bits;
        }
      }
    }
  }

  $hmmRes->domTC($lowDom);
  $hmmRes->seqTC($lowSeq);
  $hmmRes->domNC($highDom);
  $hmmRes->seqNC($highSeq);

  #Print the domains to the scores file
  foreach my $u ( sort { $b->bits <=> $a->bits } @sigUnits ) {
    print $fh
      sprintf( "%.1f %s/%s-%s %s-%s\n", $u->bits, $u->name, $u->envFrom, $u->envTo, $u->seqFrom, $u->seqTo );
  }
  close($fh);

}

#-------------------------------------------------------------------------------

#TODO - write _readAlign

=head2 _readAlign

  Title    :
  Usage    :  
  Function :
  Args     :
  Returns  :
  
=cut

sub _readAlign {
  my ( $self, $fh, $hmmRes ) = @_;

}

#Parse the alignment section
#if($pp){

#}else{
#  while(<HS>){
#    last if(/^\/\//)
#  }
#}



sub _readFooter {
  my($self, $fh, $hmmRes ) = @_;
 
  # We are going to parse something like this! 
    
  #  Internal pipeline statistics summary:
#-------------------------------------
#Query sequence(s):                         1  (360 residues)
#Target model(s):                           7  (836 nodes)
#Passed MSV filter:                         2  (0.285714); expected 0.1 (0.02)
#Passed Vit filter:                         1  (0.142857); expected 0.0 (0.001)
#Passed Fwd filter:                         1  (0.142857); expected 0.0 (1e-05)
#Initial search space (Z):                  7  [actual number of targets]
#Domain search space  (domZ):               1  [number of targets reported over threshold]
## CPU time: 0.00u 0.00s 00:00:00.00 Elapsed: 00:00:00
## Mc/sec: inf
#//  

  while(<$fh>){
    if(/\/\//){
      last;
    }
  }
}


#Parse the internal summary section
#Internal statistics summary:
#----------------------------
#Query HMM(s):      1  (0 nodes)
#Target sequences:  5323441  (0 residues)
#Passed MSV filter: 116519  (-37389918065567040729448769671768824784852036328367855636063687997915136.000; expected 19991592792512146725679052970637918208.000)
#Passed Vit filter: 7579  (-0.0000; expected -35982214160587876085407389642471051723332987952235753317595472501307733302049608744822636544.0000)
#Passed Fwd filter: 1687  (8.3e+165; expected -7.5e-266)
#Mc/sec:            828.85
# CPU time: 115.36u 4.45s 00:01:59.81 Elapsed: 00:03:01

#sub writeHMMSearch {
#  my ( $self, $hmmRes ) = @_;
#  my $fh;
#  open($fh, ">".$self->outfile."\n");
#
#  $self->_writeHeader($fh, $hmmRes);
#  $self->_writeSeqHits( $fh, $hmmRes);
#  $self->_writeDomHits( $fh, $hmmRes);
#  $self->_writeAlign( $fh, $hmmRes) if($self->align);
#  $self->_writeInternalSummary( $fh, $hmmRes);
#}
#sub mergeHMMSearch {
#  my ( $self, $filenames ) = @_;
#}


sub write_ascii_out {

    my ($self, $HMMResults, $fh, $scanData, $e_seq, $e_dom, $b_seq, $b_dom) = @_;


    $scanData->{_max_seqname} = 20 unless($scanData->{_max_seqname} or $scanData->{_max_seqname} < 1);
    
    my $ga;

    if($e_seq or $e_dom) {
	$e_seq = $e_dom unless($e_seq);
	$e_dom = "10" unless($e_dom);
    } 
    elsif($b_seq or $b_dom) {
	$b_seq = $b_dom unless($b_seq);
	$b_dom = "0" unless($b_dom);
    }
    else {
	$ga = 1;
    }


    foreach my $unit ( sort { $a->seqFrom <=> $b->seqFrom } @{ $HMMResults->units } ) {    

        if($unit->name =~ /Pfam\-B/) {

	    next unless($HMMResults->seqs->{$unit->name}->evalue <= "0.001" and $unit->evalue <= "0.001");


	    printf $fh "%-".$scanData->{_max_seqname}."s %6d %6d %6d %6d %-11s %-16s %7s %5d %5d %5d %8s %9s %3s %-8s\n",
	    $HMMResults->seqName,
	    $unit->seqFrom,
	    $unit->seqTo,
	    $unit->envFrom,
	    $unit->envTo,
	    $scanData->{_accmap}->{ $unit->name },
	    $unit->name,
	    "Pfam-B",
	    $unit->hmmFrom,
	    $unit->hmmTo,
	    $scanData->{_model_len}->{ $unit->name },
	    $unit->bits,
	    $unit->evalue,
	    "NA",
	    "NA";


	}
	else {

            #Filter results based on thresholds
	    if($ga) {
		next unless($unit->sig);
	    }
	    if($e_seq) {
		next unless($HMMResults->seqs->{$unit->name}->evalue <= $e_seq and $unit->evalue <= $e_dom);
	    }
	    if($b_seq) {
		
		next unless($HMMResults->seqs->{$unit->name}->bits >= $b_seq and $unit->bits >= $b_dom);
	    }
	    
	    my $clan = $scanData->{_clanmap}->{ $unit->name } || "No_clan";
	    
	    
	    printf $fh "%-".$scanData->{_max_seqname}."s %6d %6d %6d %6d %-11s %-16s %7s %5d %5d %5d %8s %9s %3d %-8s ",
	    $HMMResults->seqName,
	    $unit->seqFrom,
	    $unit->seqTo,
	    $unit->envFrom,
	    $unit->envTo,
	    $scanData->{_accmap}->{ $unit->name },
	    $unit->name,
	    $scanData->{_type}->{ $unit->name },
	    $unit->hmmFrom,
	    $unit->hmmTo,
	    $scanData->{_model_len}->{ $unit->name },
	    $unit->bits,
	    $unit->evalue,
	    $unit->sig, 
	    $clan;
	
	    
	    if($unit->{'act_site'}) {
		local $" = ",";
		print $fh "predicted_active_site[@{$unit->{'act_site'}}]";
	    }
	
	    if($scanData->{_translate}){
	      my $strand = '?';
	      my $start = '-';
	      my $end   = '-';
	      if(exists($scanData->{_orf}->{$HMMResults->seqName})){
	       $strand = $scanData->{_orf}->{$HMMResults->seqName}->{strand};  
	       if($strand eq '+'){
	         $start = $scanData->{_orf}->{$HMMResults->seqName}->{start} + ($unit->envFrom * 3) - 3;
	         $end = $scanData->{_orf}->{$HMMResults->seqName}->{start} + ($unit->envTo * 3) - 3;
	       }elsif($strand eq '-'){
	         $start = $scanData->{_orf}->{$HMMResults->seqName}->{start} - ($unit->envFrom * 3) + 3;
           $end = $scanData->{_orf}->{$HMMResults->seqName}->{start} - ($unit->envTo * 3) + 3;
	       }
	      }
	      print $fh "$strand $start $end";
	    }
	
	    print $fh "\n";
	}

	if($scanData->{_align}){
	    print $fh sprintf( "%-10s %s\n", "#HMM",   $unit->hmmalign->{hmm} );
	    print $fh sprintf( "%-10s %s\n", "#MATCH", $unit->hmmalign->{match} );
	    print $fh sprintf( "%-10s %s\n", "#PP",   $unit->hmmalign->{pp});
	    print $fh sprintf( "%-10s %s\n", "#SEQ",   $unit->hmmalign->{seq});
	    print $fh sprintf( "%-10s %s\n", "#CS",   $unit->hmmalign->{cs}) if($unit->hmmalign->{cs});
	}
	
    }
    
}

1;
