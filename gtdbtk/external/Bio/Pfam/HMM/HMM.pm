# HMM.pm
#
# Author:        finnr
# Maintainer:    $Id: HMM.pm,v 1.1 2009-10-08 12:27:28 jt6 Exp $
# Version:       $Revision: 1.1 $
# Created:       Nov 24, 2008
# Last Modified: $Date: 2009-10-08 12:27:28 $
=head1 NAME

Template - a short description of the class

=cut

package Bio::Pfam::HMM::HMM;

=head1 DESCRIPTION

A more detailed description of what this class does and how it does it.

$Id: HMM.pm,v 1.1 2009-10-08 12:27:28 jt6 Exp $

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

#-------------------------------------------------------------------------------

=head1 METHODS

=cut


subtype 'hmmVersion'
  => as Str
  => where { $_ =~ m/^HMMER3\/f\s+\[3\.\d+[ab](\d+)?\s+|\s+\[.*\]/ }
  => message { "|$_| does not look like as HMMER3 version" };

subtype 'hmmName'
  => as Str
  => where { $_ =~ m/\S{1,15}/ }
  => message { "|$_| does not look like Pfam name or SEED" };

subtype 'hmmAcc'
  => as Str
  => where { $_ =~ m/PF\d{5}/ }
  => message { "|$_| does not look like Pfam accession" };


subtype 'hmmAlpha'
  => as Str
  => where { $_ eq 'amino' or $_ eq 'nucleic'  }
  => message { "|$_| does not look like a HMMER3 alphabet" };


subtype 'hmmMsvStats'
  => as HashRef
  => where { defined ($_->{mu}) and defined ($_->{lambda}) and ($_->{lambda} <= 0.8) and ($_->{lambda} >= 0.5) }
  => message { "Mu |$_->{mu}| and lambda |$_->{lambda}| must be defined and lambda must be between 0.5 and 0.8" };

  
subtype 'hmmViterbiStats'
  => as HashRef
  => where { defined ($_->{mu}) and defined ($_->{lambda}) and ($_->{lambda} <= 0.8) and ($_->{lambda} >= 0.5) }
  => message { "Mu |$_->{mu}| and lambda |$_->{lambda}| must be defined and lambda must be between 0.5 and 0.8" };

subtype 'hmmForwardStats'
  => as HashRef
  => where { defined ($_->{tau}) and defined ($_->{lambda}) and ($_->{lambda} <= 0.8) and ($_->{lambda} >= 0.5) }
  => message { "Tau |$_->{tau}| and lambda |$_->{lambda}| must be defined and lambda must be between 0.5 and 0.8" };
  
has 'version' => (
  isa      => 'hmmVersion',
  is       => 'rw',
  required => 1
);

has 'name' => (
  isa => 'hmmName',
  is  => 'rw',
  required => 1
);

has 'accession' => (
  isa => 'hmmAcc',
  is  => 'rw'
);

has 'description' => (
  isa => 'Str',
  is  => 'rw'
);

has 'length' => (
  isa  => 'Int',
  is   => 'rw',
  required => 1
);

has 'alpha' => (
  isa  => 'hmmAlpha',
  is   => 'rw',
  required => 1,
);

has 'rf' => (
  isa  => 'Bool',
  is   => 'rw',
  required => 1
);

has 'mm' => (
  isa  => 'Bool',
  is   => 'rw',
);

has 'cons' => (
  isa  => 'Bool',
  is   => 'rw',
);

has 'cs' => (
  isa  => 'Bool',
  is   => 'rw',
  required => 1
);

has 'map' => (
  isa  => 'Bool',
  is   => 'rw',
  required => 1
);

has 'date' => (
  isa  => 'Str',
  is   => 'rw',
  required => 1
);

has 'buildLine' => (
  isa => 'HashRef[Str]',
  is  => 'rw',
  required => 1,
  default => sub { {} },
);

has 'searchMethod' => (
  isa => 'Str',
  is  => 'rw',
);

has 'nSeq' => (
  isa  => 'Int',
  is   => 'rw',
  required => 1
);

has 'msvStats' => (
  isa => 'hmmMsvStats',
  is => 'rw',
  required => 1
);

has 'viterbiStats' => (
  isa => 'hmmViterbiStats',
  is => 'rw',
  required => 1
);

has 'forwardStats' => (
  isa => 'hmmForwardStats',
  is => 'rw',
  required => 1
);


has 'effn' => (
  isa  => 'Num',
  is   => 'rw',
  required => 1
);

has 'cksum' => (
  isa => 'Int',
  is  => 'rw',
  required => 1
);

has 'seqGA' => (
  isa => 'Num',
  is  => 'rw',
);

has 'domGA' => (
  isa => 'Num',
  is  => 'rw',
);

has 'seqTC' => (
  isa => 'Num',
  is  => 'rw',
);

has 'domTC' => (
  isa => 'Num',
  is  => 'rw',
);

has 'seqNC' => (
  isa => 'Num',
  is  => 'rw',
);

has 'domNC' => (
  isa => 'Num',
  is  => 'rw',
);

has 'emissionLines' => (
  isa => 'ArrayRef[ArrayRef]',
  is  => 'rw',
  default => sub { [] },
);

has 'mapPos'=> (
  isa => 'ArrayRef[Int]',
  is  => 'rw',
  default => sub{ [] }
);

has 'compLines' => (
  isa => 'ArrayRef[Str]',
  is  => 'rw',
  default => sub { [] },
);

__PACKAGE__->meta->make_immutable;
1;

