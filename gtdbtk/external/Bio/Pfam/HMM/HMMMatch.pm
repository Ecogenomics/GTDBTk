
package Bio::Pfam::HMM::HMMMatch;

use strict;
use warnings;

use Moose;
use Moose::Util::TypeConstraints;


subtype 'evalue'
  => as Str
  => where { $_ =~ m/^(\d+(\.\d+){0,1}e[+|-]\d+|\d+\.\d+|\d+)$/ }
  => message { "$_ does not look like an evalue" };

has 'evalue' => (
  isa       => 'evalue', 
  is        => 'rw',
  required  => 1
);

has 'bits' => (
  isa => 'Str',
  is  => 'rw',
  required => 1
);

has 'name' => (
  isa => 'Str',
  is  => 'rw',
  required => 1
);

has bias => (
  isa => 'Num',
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
