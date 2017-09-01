package Bio::Pfam::Scan::Seq;

use strict;
use warnings;

use Bio::LocatableSeq;
use Bio::Seq::RichSeq;

use base qw(Bio::LocatableSeq Bio::Seq::RichSeq);

sub new {
  my($class, %params ) = @_;
  my( $id, $start, $end, $seq) =
      (
       ($params{'-ID'}          || $params{'-id'}),
       ($params{'-START'}       || $params{'-start'}),
       ($params{'-END'}         || $params{'-end'}),
       ($params{'-SEQ'}         || $params{'-seq'}),
       );

  my $self = $class->SUPER::new( %params );  # this is Bio::Pfam::Root
                      # so we have to set Bio::LocatableSeq fields ourself




  $self->id( $id );
  $self->start( $start );
  $self->end( $end );
  $self->seq( $seq );


  return $self; # success - we hope!
}

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

1
