# Bio::Pfam::HMM::HMMResults.pm
#
# Author:        finnr
# Maintainer:    $Id: HMMResults.pm,v 1.3 2009-12-15 14:38:08 jt6 Exp $
# Version:       $Revision: 1.3 $
# Created:       Nov 19, 2008
# Last Modified: $Date: 2009-12-15 14:38:08 $

=head1 NAME

Bio::Pfam::HMM::HMMResults - A object to represents the results from hmmsearch

=cut

package Bio::Pfam::HMM::HMMResults;

=head1 DESCRIPTION

A more detailed description of what this class does and how it does it.

$Id: HMMResults.pm,v 1.3 2009-12-15 14:38:08 jt6 Exp $

=head1 COPYRIGHT

File: Bio::Pfam::HMM::HMMResults.pm

Copyright (c) 2007: Genome Research Ltd.

Authors: Rob Finn (rdf@sanger.ac.uk)

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
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 
=cut

use strict;
use warnings;

use Moose;
use Moose::Util::TypeConstraints;
use Bio::Pfam::HMM::HMMSequence;
use Bio::Pfam::HMM::HMMUnit;

#
#-------------------------------------------------------------------------------
# Attributes

has 'hmmerVersion' => (
  isa       => 'Str', 
  is        => 'rw', 
);

has 'hmmName' => (
  isa => 'Str',
  is  => 'rw'
);

has 'seqDB' => (
  isa => 'Str',
  is  => 'rw'
);

has hmmLength => (
  isa => 'Int',
  is  => 'rw'
);

has 'thisFile' => (
  isa => 'Str',
  is  => 'rw'
);

has seedName => (
  isa => 'Str',
  is  => 'rw'
);

has 'seqs' => (
  isa => 'HashRef',
  is  => 'rw',
  default => sub { {} },
);

has 'units' => (
  isa => 'ArrayRef',
  is  => 'rw',
  default => sub { [] },
);

has 'domThr' => (
  isa => 'Num',
  is   => 'rw',
  default => '25.0'
);

has 'seqThr' => (
  isa => 'Num',
  is  => 'rw',
  default => '25.0'
);

has 'evalueThr' => (
  isa => 'Num',
  is  => 'rw'
);

has 'domTC' => (
  isa => 'Num',
  is  => 'rw'
);

has 'seqTC' => (
  isa => 'Num',
  is  => 'rw'
);

has 'domNC' => (
  isa => 'Num',
  is  => 'rw'
);

has 'seqNC' => (
  isa => 'Num',
  is  => 'rw'
);

has 'randSeedNum' => (
  isa => 'Int',
  is  => 'rw'
);

has 'description' => (
  isa => 'Str',
  is  => 'rw'
);

has 'seqName' => (
  isa => 'Str',
  is  => 'rw'
);

has 'seqLength' => (
  isa => 'Int',
  is  => 'rw'
);


has 'eof' => (
  isa => 'Int',
  is  => 'rw',
  default => 0
);

has 'program' => (
  isa => 'Str',
  is  => 'rw'
);

=head1 METHODS

=head2 addHMMSeq 

  Title    : addHMMSeq
  Usage    : $hmmRes->addHMMSeq( $hmmSeqObj )
  Function : Adds a Bio::Pfam::HMM::HMMSequence object to the results object 
  Args     : A Bio::Pfam::HMM::HMMSequence object
  Returns  : nothing
  
=cut

sub addHMMSeq {
  my( $self, $hmmSeq ) = @_;

  unless($hmmSeq->isa('Bio::Pfam::HMM::HMMSequence')){
    die 'Trying to add a non Bio::Pfam::HMM::HMMSequence object';  
  }

  if($self->seqs){
    if($self->seqs->{$hmmSeq->name}){
      die "Trying to add the same sequence twice"; 
    }
  }
  
  $self->seqs->{$hmmSeq->name} = $hmmSeq;
}


=head2 eachHMMSeq 

  Title    : eachHMMSeq
  Usage    : my @seqs = $hmmRes->eachHMMSeq
  Function : Returns an array reference containing the references to all of the
           : Bio::Pfam::HMM::HMMSequence objects stored in the HMMResults object.
  Args     : None
  Returns  : Array reference
  
=cut

sub eachHMMSeq {
  my ($self ) = @_;
  my @seqs;
  my $seqRefs = $self->seqs;
  foreach my $n (keys %{ $seqRefs }){
    push(@seqs, $seqRefs->{$n});
  } 
  return(\@seqs);
}

#-------------------------------------------------------------------------------

=head2 addHMMUnit 

  Title    : addHMMUnit
  Usage    : $hmmRes 
  Function : Adds HMM units (the actual region hit) to the HMMSequence in the object
           : and for convenience to the the results sets. All we store are duplicates
           : of the references. 
  Args     : A Bio::Pfam::HMM:HMMUnit
  Returns  : Nothing
  
=cut

sub addHMMUnit {
  my ($self, $hmmUnit) = @_;
  
  unless($hmmUnit->isa('Bio::Pfam::HMM::HMMUnit')){
    die "Trying to add an non-Bio::Pfam::HMM::HMMUnit\n";
  }
  
  if($self->seqs){
    if($self->seqs->{$hmmUnit->name}){
      $self->seqs->{$hmmUnit->name}->addHMMUnit($hmmUnit);
    }else{
      warn "Could not add hmmUnit as the sequence has not been added\n";
    }
  }
  
  #More conveinence we store the point to the hmmunit in an array
  push(@{$self->units},$hmmUnit); 
}


#-------------------------------------------------------------------------------

=head2  domainBitsCutoffFromEvalue

  Title    : domainBitsCutoffFromEvalue
  Usage    : $hmmRes->domainBitsCutoffFromEvalue(0.01)
  Function : From the supplied evalue, it scans through all of the evalues in the results
           : and calulates the bits score.  
  Args     : An evalue. 
  Returns  : A bits score. If no evalue is specified, returns nothing
  
=cut

sub domainBitsCutoffFromEvalue {
  my ($self, $eval) = @_;
  my ($dom,$prev,@doms,$cutoff,$sep,$seen);
  
  unless(defined ($eval) ){
    warn "No evalue specified\n";
    return; 
  }


  $seen = 0;
  foreach $_ ( sort { $b->bits <=> $a->bits } @{$self->units}, @{$self->eachHMMSeq} ) {
    if( $_->evalue > $eval ) {
	    $seen = 1;
	    $dom = $_;
	    last;
	  } 
	  $prev = $_;
  }
    
  if( ! defined $prev || $seen == 0) {
	  carp("Evalue is either above or below the list...");
	  return undef;
  }

  $sep = $prev->bits - $dom->bits ;
    
  if( $sep < 1 ) {
	    return $prev->bits();
  }
  if( $dom->bits < 25 && $prev->bits > 25 ) {
	  return 25;
  }

  return $dom->bits + sprintf("%.1f",$sep/2);
}


#-------------------------------------------------------------------------------

=head2 lowestTrue

  Title    :
  Usage    :  
  Function :
  Args     :
  Returns  :
  
=cut

sub lowestTrue {
  my $self = shift;
  
  unless($self->domTC && $self->seqTC) {
    unless($self->domThr and $self->seqThr){
      die "Could not define TC as I am missing a threshold\n";
    }
    #Set it wildly high!
    my ($lowSeq, $lowDom);
    $lowSeq = $lowDom = 999999.99;  
    
    foreach my $seqId (keys %{$self->seqs} ){ 
      if($self->seqs->{$seqId}->bits >= $self->seqThr){
        #Is this the lowest sequence thresh
        if($self->seqs->{$seqId}->bits < $lowSeq){
          $lowSeq = $self->seqs->{$seqId}->bits; 
        }
        #For each of the regions found on the sequence, look to see if the match is great
        #than the domain threshold.  If it is, is it lower than we we have seen previously    
        foreach my $unit (@{ $self->seqs->{$seqId}->hmmUnits } ){
          if( $unit->bits() >= $self->domThr && $unit->bits() < $lowDom ) {
		        $lowDom  = $unit->bits;
          }
        } 
      }
      
    }
    $self->domTC($lowDom);
    $self->seqTC($lowSeq);
  }
  return($self->seqTC, $self->domTC);
}

#-------------------------------------------------------------------------------

=head2 highestNoise

  Title    :
  Usage    :  
  Function :
  Args     :
  Returns  :
  
=cut

sub highestNoise {
  my $self = shift;
  
  #See if it is already set
  unless($self->domNC && $self->seqNC) {
    unless($self->domThr and $self->seqThr){
      die "Could not define TC as I am missing a threshold\n";
    }
    
    #Set it wildly low
    my ($highSeq, $highDom);
    $highSeq = $highDom = -999999.99;  
    
    foreach my $seqId (keys %{$self->seqs} ){
      
      if($self->seqs->{$seqId}->bits < $self->seqThr){
        #Is this the highest sequence thres below the cut-off
        if($self->seqs->{$seqId}->bits > $highSeq){
          $highSeq = $self->seqs->{$seqId}->bits; 
        }
      }
      
      #For each of the regions found on the sequence, look to see if the match is great
      #than the domain threshold.  If it is, is it lower than we we have seen previously    
      foreach my $unit (@{ $self->seqs->{$seqId}->hmmUnits } ){
        if( $unit->bits < $self->domThr && $unit->bits > $highDom ) {
		      $highDom  = $unit->bits;
        }
      }
    }
    $self->domNC($highDom);
    $self->seqNC($highSeq); 
  }

  return($self->seqNC, $self->domNC);
}


sub applyEdits {
  my ($self, $edits) = @_;
  foreach my $e (@$edits){
    #{ seq => $1, oldFrom => $2, oldTo => $3, newFrom => $5, newTo => $6 }  
    if($self->seqs->{$e->{seq}}){
      my $matched = 0;
      foreach my $u (@{ $self->seqs->{ $e->{seq} }->hmmUnits }){
          if($u->envFrom == $e->{oldFrom} and 
             $u->envTo == $e->{oldTo}) {
            $matched = 1;
            if(defined $e->{newFrom} and $e->{newTo}){
              
              #Modify the start end positions
              $u->envFrom($e->{newFrom});
              $u->envTo($e->{newTo});
              
              #Check that the ali-positions are still okay
              if($u->seqFrom < $e->{newFrom}){
                $u->seqFrom($e->{newFrom});  
              }
              if($u->seqTo > $e->{newTo}){
                $u->seqTo($e->{newTo});  
              }
            }else{
              #Set the score so low it will never get in the align
              $u->bits(-999999.99);
            }
            last;    
          }  
      }
      unless($matched){
        warn $e->{seq}."/".$e->{oldFrom}."-".$e->{oldTo}." does not appear in the list of hmm units - bad ED line\n";
      }
    }else{
      warn $e->{seq}." does not appear in the list of hmm units - bad ED line\n";  
    }
  }
  
}

sub remove_overlaps_by_clan {

    my ($self, $clanmap, $nested) = @_;

    my $new = Bio::Pfam::HMM::HMMResults->new;
    $new->seqName($self->seqName);
   
    foreach my $unit ( sort { $a->evalue <=> $b->evalue }  @{ $self->units } ) {

        #check if it overlaps before adding
	my $o;
	
	foreach my $u ( @{ $new->units } ) {
	    
	    if( exists($clanmap->{$unit->name}) and exists($clanmap->{$u->name}) and ($clanmap->{$unit->name} eq $clanmap->{$u->name}) ) {
		if( overlap( $unit, $u ) ) {
		    if(exists($$nested{$unit->name}{$u->name})) {
			next;
		    }
		    else {
			$o=1;
			last;
		    }
		}
		
	    }
	}
	unless($o) {
	    if(! $new->seqs->{$unit->name}) {
		
		$new->addHMMSeq( Bio::Pfam::HMM::HMMSequence->new( {  name       => $self->seqs->{$unit->name}->name,
								      desc       => $self->seqs->{$unit->name}->desc,
								      bits       => $self->seqs->{$unit->name}->bits,
								      evalue     => $self->seqs->{$unit->name}->evalue,
								      numberHits => $self->seqs->{$unit->name}->numberHits}) );
	    
	    }
	    $new->addHMMUnit($unit);
	}

    }
    return $new;
}



sub overlap {
    # does unit1 overlap with unit2?
    my $unit1 = shift;
    my $unit2 = shift;
    my( $u1, $u2 ) = sort { $a->seqFrom <=> $b->seqFrom } ( $unit1, $unit2 );


    if( $u2->seqFrom <= $u1->seqTo ) {
        return 1;
    } 

    return 0;
}


sub results {
  my ( $self, $pfamScanData, $e_value ) = @_;

  my @results = ();
  foreach my $unit ( sort { $a->seqFrom <=> $b->seqFrom } @{ $self->units } ) {    

    my $pfamB = $unit->name =~ /^Pfam-B/;

    #Filter results based on thresholds
    if ( $unit->name =~ /^Pfam-B/ ) {
      next unless ( $self->seqs->{$unit->name}->evalue <= 0.001 and $unit->evalue <= 0.001 );
      $pfamB = 1;
    }
    else {	  
      if ( $e_value ) {
        next unless ( $self->seqs->{$unit->name}->evalue <= $e_value and $unit->evalue <= $e_value ) ;
      } 
      else {
       next unless $unit->sig;
      }
    }

    push @results, {
      seq          => { from => $unit->seqFrom,
                        to   => $unit->seqTo,
                        name => $self->seqName },
      env          => { from => $unit->envFrom,
                        to   => $unit->envTo },

      hmm          => { from => $unit->hmmFrom,
                        to   => $unit->hmmTo },

      model_length => $pfamScanData->{_model_len}->{ $unit->name },
      bits         => $unit->bits,
      evalue       => $unit->evalue,
      acc          => $pfamScanData->{_accmap}->{ $unit->name },
      name         => $unit->name,
      desc         => $pfamScanData->{_desc}->{ $unit->name },
      type         => $pfamB ? undef : $pfamScanData->{_type}->{ $unit->name },
      clan         => $pfamB ? undef : 
                         $pfamScanData->{_clanmap}->{ $unit->name } || 'No_clan',

      act_site     => $pfamB ? undef : $unit->{act_site},
      sig          => $pfamB ? "NA" : $unit->sig,
      align        => [ sprintf( '#HMM       %s', $unit->hmmalign->{hmm} ),
                        sprintf( '#MATCH     %s', $unit->hmmalign->{match} ),
                        sprintf( '#PP        %s', $unit->hmmalign->{pp} ),
                        sprintf( '#SEQ       %s', $unit->hmmalign->{seq} ) ]
    };
  }

  return \@results;
}

1;
