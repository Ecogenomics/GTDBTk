# HMM.pm
#
# Author:        finnr
# Maintainer:    $Id: HMMIO.pm,v 1.3 2010-01-12 17:00:26 jm14 Exp $
# Version:       $Revision: 1.3 $
# Created:       Nov 24, 2008
# Last Modified: $Date: 2010-01-12 17:00:26 $
=head1 NAME

Template - a short description of the class

=cut

package Bio::Pfam::HMM::HMMIO;

=head1 DESCRIPTION

A more detailed description of what this class does and how it does it.

$Id: HMMIO.pm,v 1.3 2010-01-12 17:00:26 jm14 Exp $

=head1 COPYRIGHT

File: HMM.pm

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
use Moose::Util::TypeConstraints;
use Carp;
use Bio::Pfam::HMM::HMM;

#-------------------------------------------------------------------------------

=head1 METHODS

=cut
sub readHMM {
  my ($this, $hmm) = @_;
  
  unless($hmm){
    confess("No HMM passed in!"); 
  }
  chomp($hmm);
  
  my @input;
  if(ref($hmm) eq 'GLOB'){
    @input = <$hmm>;
  }elsif($hmm !~ /\n/ and -e $hmm and -s $hmm){
    #Assume that we have a filename and try and open it;
    open(HMM, $hmm) || confess("Could not open $hmm:[$!]");
    @input = <HMM>;
  }else{
    @input = split(/\n/, $hmm); 
  }
  
  
  
  #Parse the header section!  
  #HMMER3/f [3.1b1 | May 2013]
  #NAME  SEED
  #ACC   PF000001.1
  #DESC  A description
  #LENG  55
  #ALPH  amino
  #RF    no
  #MM    no
  #CONS  yes
  #CS    no
  #MAP   yes
  #DATE  Fri Nov 21 09:58:16 2008
  #COM   [1] /Users/finnr/Work/Software/hmmer-3.0.20081101/bin/hmmbuild -o hmmbuild.log HMM SEED
  #NSEQ  279
  #EFFN  4.966292
  #STATS LOCAL MSV      -11.4716  0.69948
  #STATS LOCAL VITERBI  -12.3713  0.69948
  #STATS LOCAL FORWARD   -5.5807  0.69948

  #To add GA, TC, NC, CKSUM, DESC
  my($objHash);
  my $i =0; 
  foreach ( @input ){  
    if(my ($version) = $_ =~ /(HMMER3.*)/){
      $objHash->{version} = $version;
    }elsif(my ($acc) = $_ =~ /^ACC\s+(PF\d+\.\d+)$/){
      $objHash->{accession} = $acc;
    }elsif(/NAME\s+(\S+)/){ 
      $objHash->{name} =  $1 ;
    }elsif(/DESC\s+(.*)/){ 
      $objHash->{description} =   $1 ;
    }elsif(my ($length) = $_ =~ /^LENG\s+(\d+)/){
      $objHash->{length} = $length;
    }elsif( my ($alpha) = $_ =~ /^ALPH\s+(\S+)/){
      $objHash->{alpha} = $alpha;
    }elsif( my ($rf) = $_ =~ /^RF\s+(no|yes)/){
      $objHash->{rf} = ($rf eq "no") ? 0 : 1; 
    }elsif( my ($mm) = $_ =~ /^MM\s+(no|yes)/){
      $objHash->{mm} = ($mm eq "no") ? 0 : 1; 
    }elsif( my ($cons) = $_ =~ /^CONS\s+(no|yes)/){
      $objHash->{cons} = ($cons eq "no") ? 0 : 1; 
    }elsif(my ($cs) = $_ =~ /^CS\s+(no|yes)/ ){
      $objHash->{cs} =  ($cs eq "no") ? 0 : 1; 
    }elsif(my ($map) = $_ =~ /^MAP\s+(no|yes)/){
      $objHash->{map} = ($map eq "no") ? 0 : 1; 
    }elsif(my ($date) = $_ =~ /^DATE\s+(.*)/){
      $objHash->{date} =  $date; 
    }elsif(my ($sm) = $_ =~ /^SM\s+(.*)/){
      $objHash->{searchMethod} =  $sm; 
    
    }elsif(my ($options, $hmmName, $alignName) = $input[$i] =~ /^BM.*hmmbuild(.*)? (\S+) (\S+)$/){
      $objHash->{buildLine} = { cmd     => 'hmmbuild', 
                                options => $options,
                                name    => $hmmName,
                                align   => $alignName } ;
    }elsif(my($noSeqs) = $_ =~ /^NSEQ\s+(\d+)/){
      $objHash->{nSeq} = $noSeqs;
    }elsif( my($effn) = $_ =~ /^EFFN\s+(\d+\.\d+)/){
       #EFFN  4.966292
      $objHash->{effn} =  $effn ;
    }elsif( my ( $cksum ) = $_ =~ /^CKSUM (\d+)/){
      $objHash->{cksum} = $cksum ;
    }elsif(/GA\s+(\S+)\s+(\S+)\;/){ 
      $objHash->{seqGA} = $1;
      $objHash->{domGA} = $2;
    }elsif(/TC\s+(\S+)\s+(\S+)\;/){ 
      $objHash->{seqTC} = $1;
      $objHash->{domTC} = $2;
    }elsif(/NC\s+(\S+)\s+(\S+)\;/){ 
      $objHash->{seqNC} = $1;
      $objHash->{domNC} = $2;
    }elsif( my ($msv_mu, $msv_lambda ) = $_ =~ /^STATS LOCAL MSV\s+(\S+)\s+(0\.\d+)/){
	$objHash->{msvStats} = { mu => $msv_mu, lambda => $msv_lambda};
    }elsif( my ($viterbi_mu, $viterbi_lambda ) = $_ =~ /^STATS LOCAL VITERBI\s+(\S+)\s+(0\.\d+)/){
	$objHash->{viterbiStats} = { mu => $viterbi_mu, lambda => $viterbi_lambda };
    }elsif( my ($forward_tau, $forward_lambda ) = $_ =~ /^STATS LOCAL FORWARD\s+(\S+)\s+(0\.\d+)/){
	$objHash->{forwardStats} = {tau => $forward_tau, lambda => $forward_lambda};
    }elsif( $_ =~ /^HMM\s+A/){
      last;
    }else{
      confess("Got a bad HMM header line:$input[$i]\n"); 
    }
    $i++;
  }

  my $hmmObj = Bio::Pfam::HMM::HMM->new($objHash);
  
  
  #The next two lines are stand lines
  #HMM          A        C        D        E        F        G        H        I        K        L        M        N        P        Q        R        S        T        V        W        Y   
  #          m->m     m->i     m->d     i->m     i->i     d->m     d->d
  $i++;


  #Add the comp line
  for( my $line = 0; $line <=2; $line++){
    $i++;
    my @l = split(/\s+/, $input[$i]);
    my @c;
    if($line == 0 ){
      @c = @l[2..21];
    }elsif( $line == 1){
      @c = @l[1..20]; 
    }elsif($line == 2){
      @c = @l[1..7];  
    }
    $hmmObj->compLines->[$line] = \@c;
  }
  

  
  for(my $pos = 0; $pos < $hmmObj->length; $pos++){
    #There are three lines per position - match emission line, insert emission line, state transition line
    for( my $line = 0; $line <=2; $line++){
      $i++;
      my @l = split(/\s+/, $input[$i]);
      my @e;
      if($line == 0 ){
        @e = @l[2..21];
        if($hmmObj->map){
           $hmmObj->mapPos->[$pos] = $l[22];
        } 
      }elsif( $line == 1){
        @e = @l[1..20]; 
      }elsif($line == 2){
        @e = @l[1..7];  
      }
      $hmmObj->emissionLines->[$pos]->[$line] = \@e;
    }
  }
  
  if($input[$i++] =~ /^\/\/$/){
    confess("Expected file terminator //, but got $input[$i]\n"); 
  }
 
  #No veryifiy that we have COMP line and the the number of emissionlines is equivalent to length 
  unless(scalar( @{ $hmmObj->emissionLines } ) == $hmmObj->length){
    confess("Number of emssionLines does not match the length of the model, got ".scalar( @{ $hmmObj->emissionLines} ).
        " expected ".$hmmObj->length);  
  }
  
  unless($hmmObj->compLines){
    confess("No compLine set on HMM"); 
  }
  
  if($hmmObj->map){
    unless(scalar(@{$hmmObj->mapPos}) == $hmmObj->length ){
      confess("HMM object had map set, but the number of map positions does not match the length of the HMM");  
    }; 
  }
  return $hmmObj;  
}



sub writeHMM {
  my ($this, $hmm, $hmmObj) = @_;
  
  unless($hmm){
    confess("No HMM out file passed in!"); 
  }
  
  unless(ref($hmm) eq 'GLOB'){
    my $hmmFile = $hmm;
    $hmm = undef;
    #Assume that we have a filename and try and open it;
    open($hmm, ">$hmmFile") || confess("Could not open $hmmFile:[$!]");
  }

  print  $hmm $hmmObj->version."\n";
  printf $hmm ("%-5s %s\n", "NAME", $hmmObj->name);
  printf $hmm ("%-5s %s\n", "ACC",  $hmmObj->accession) if($hmmObj->accession);
  printf $hmm ("%-5s %s\n", "DESC", $hmmObj->description) if($hmmObj->description);
  printf $hmm ("%-5s %d\n", "LENG", $hmmObj->length);
  printf $hmm ("%-5s %s\n", "ALPH", $hmmObj->alpha);
  printf $hmm ("%-5s %s\n", "RF", ($hmmObj->rf ? "yes" : "no"));
  printf $hmm ("%-5s %s\n", "MM", ($hmmObj->mm ? "yes" : "no"));
  printf $hmm ("%-5s %s\n", "CONS", ($hmmObj->cons ? "yes" : "no"));
  printf $hmm ("%-5s %s\n", "CS", ($hmmObj->cs ? "yes" : "no"));
  printf $hmm ("%-5s %s\n", "MAP", ($hmmObj->map ? "yes" : "no"));
  printf $hmm ("%-5s %s\n", "DATE", $hmmObj->date);
  printf $hmm ("%-5s %d\n", "NSEQ", $hmmObj->nSeq);
  printf $hmm ("%-5s %f\n", "EFFN", $hmmObj->effn);
  printf $hmm ("%-5s %d\n", "CKSUM", $hmmObj->cksum);
  printf $hmm ("%-5s %.2f %.2f;\n", "GA", $hmmObj->seqGA, $hmmObj->domGA) if(defined($hmmObj->seqGA));
  printf $hmm ("%-5s %.2f %.2f;\n", "TC", $hmmObj->seqTC, $hmmObj->domTC) if(defined($hmmObj->seqTC));
  printf $hmm ("%-5s %.2f %.2f;\n", "NC", $hmmObj->seqNC, $hmmObj->domNC) if(defined($hmmObj->seqNC));
  
  printf $hmm ("%-5s %s %-9s %.4f  %.5f\n", "STATS", "LOCAL", "MSV", $hmmObj->msvStats->{mu},  $hmmObj->msvStats->{lambda});
  printf $hmm ("%-5s %s %-9s %.4f  %.5f\n", "STATS", "LOCAL", "VITERBI", $hmmObj->viterbiStats->{mu},  $hmmObj->viterbiStats->{lambda});
  printf $hmm ("%-5s %s %-9s %.4f  %.5f\n", "STATS", "LOCAL", "FORWARD", $hmmObj->forwardStats->{tau},  $hmmObj->forwardStats->{lambda});
  
  print $hmm <<EOF;
HMM          A        C        D        E        F        G        H        I        K        L        M        N        P        Q        R        S        T        V        W        Y   
            m->m     m->i     m->d     i->m     i->i     d->m     d->d
EOF
  
  printf $hmm ("%7s ", "COMPO");
   foreach my $s (@{$hmmObj->compLines->[0]}){
       printf $hmm ("  %7s", $s);
    }
    print $hmm "\n";
     
    print $hmm (" " x 8);
    foreach my $s (@{$hmmObj->compLines->[1]}){
       printf $hmm ("  %7s", $s);
    }
    print $hmm "\n";
    
    print $hmm (" " x 8);
    foreach my $s (@{$hmmObj->compLines->[2]}){
       printf $hmm ("  %7s", $s);
    }
    print $hmm "\n";
    
  
  my $pos = 1;
  foreach my $el (@{ $hmmObj->emissionLines }){
    printf $hmm ("%7s ", $pos);
    foreach my $s (@{$el->[0]}){
       printf $hmm ("  %7s", $s);
    }
    if($hmmObj->map){
      printf $hmm ("%7s - -\n", $hmmObj->mapPos->[$pos - 1]);  
    }else{
      print $hmm "\n";
    } 
    print $hmm (" " x 8);
    foreach my $s (@{$el->[1]}){
       printf $hmm ("  %7s", $s);
    }
    print $hmm "\n";
    print $hmm (" " x 8);
    foreach my $s (@{$el->[2]}){
       printf $hmm ("  %7s", $s);
    }
    print $hmm "\n"; 
    $pos++; 
  }
  

  print $hmm "//\n";
}



1;

