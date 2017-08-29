
package Bio::Pfam::HMM::HMMUnit;

use strict;
use warnings;

use Moose;
use Moose::Util::TypeConstraints;

extends 'Bio::Pfam::HMM::HMMMatch';

subtype 'Domain'
    => as 'Int'
    => where { $_ > 0 };
              
#coerce 'Domain'
#  => from 'Str'
#    => via {
#      my $d;
#      if(/(\d+)\/\d+/){
#        $d = $1;
#      }
#      return $d;
#    };
#    
  
#subtype 'proteinCoos'
#  => as 'Int'
#  => where { $_ > 0 && $_ < 100000 }
#  => message { 'Protein coordinates are expected to be positive and less the 100,000'};


has 'seqEvalue' => (
  isa      => 'Num',
  is       => 'rw',
);

has 'domain' => (
  isa     => 'Domain',
  is       => 'rw'
);

has 'seqFrom' => (
  isa => 'Int',
  is   => 'rw',
  required => 1
);

has 'seqTo' => (
  isa => 'Int',
  is  => 'rw',
  required => 1
);

#has 'indEvalue' => (  
#  isa => 'evalue',
#  is  => 'rw',
#  required => 1,
#);

has 'domEvalue' => (
  isa => 'evalue',
  is  => 'rw',
);

has 'hmmalign' => (
  isa => 'HashRef',
  is  => 'rw',
  default => sub { {} },
);

has 'hmmFrom' => (
  isa => 'Int',
  is   => 'rw',
  required => 1
);

has 'hmmTo' => (
  isa => 'Int',
  is  => 'rw',
  required => 1
);

has 'envFrom' => (
  isa => 'Int',
  is   => 'rw'
);

has 'envTo' => (
  isa => 'Int',
  is  => 'rw'
);

has 'coreFrom' => (
  isa => 'Str',
  is   => 'rw'
);

has 'coreTo' => (
  isa => 'Str',
  is  => 'rw'
);

has 'aliAcc' => (
  isa => 'Num',
  is  => 'rw'
);

has 'sig' => (
  isa => 'Int',
  is  => 'rw'
);


__PACKAGE__->meta->make_immutable;

=head1 COPYRIGHT

Copyright (c) 2007: Genome Research Ltd.

Authors: Rob Finn (rdf@sanger.ac.uk), John Tate (jt6@sanger.ac.uk)

This is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <http://www.gnu.org/licenses/>.

=cut

1;
