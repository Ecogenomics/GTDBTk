package Bio::Pfam::Active_site::as_search;

use strict;
use warnings;

use Bio::SeqFeature::Generic;
use Bio::SimpleAlign;
use Bio::Pfam::Scan::Seq;

=head2 find_as

 Title   : find_as
 Usage   : find_as($as_aln, $as_res, $seq_id, $seq_se, $seq_region, $family, $hmm_file)
 Function: finds active sites in a query sequence which 
           has a match to a Pfam active site family

 Returns : An array reference of active site postions
 Args    : Alignment object of active site sequences, hash of arrays containing seq ids => active site positions, 
           start-end sequence in the format "3-50", sequence region, family,  file containing all Pfam models

=cut

sub find_as {
  my ($as_aln, $as_res, $seq_id, $seq_se, $seq_region, $family, $hmm_file) = @_;


  system("hmmfetch $hmm_file $family > /tmp/hmm.$$") and die "FATAL: Problem running [hmmfetch $hmm_file $family > /tmp/hmm.$$]\n";

  $seq_id = "Query_".$seq_id; 

  my $fasta;
  foreach my $seq ($as_aln->each_seq) {
      my $s = $seq->seq;
      $s =~ s/[\-\.]//g; #Remove gaps

      $fasta .= ">" . $seq->id . "/" . $seq->start . "-" . $seq->end . "\n$s\n";
  }
  $fasta .= ">$seq_id/$seq_se\n$seq_region";
  open(SEQ, ">/tmp/seqs.$$") or die "Couldn't open file seqs.$$ $!\n";
  print SEQ $fasta;
  close SEQ;
 

  open(OUT, "hmmalign --outformat Pfam /tmp/hmm.$$ /tmp/seqs.$$ |") or die "Couldn't open fh to hmmalign $!\n";

  my $aln = new Bio::SimpleAlign;
  my ($name, $start, $end, $seq);
  while(<OUT>) {
    if( /^(\S+)\/(\d+)-(\d+)\s+(\S+)\s*/ ) {
	$name = $1;
	$start = $2;
	$end = $3;
	$seq = $4;
	    
	$aln->add_seq(Bio::Pfam::Scan::Seq->new('-seq'=>$seq, '-id'=>$name, '-start'=>$start, '-end'=>$end, '-type'=>'aligned'));
    }
  }
  close OUT;
  

  unlink "/tmp/seqs.$$";
  unlink "/tmp/hmm.$$";
  
  #Locate exp as in fam
  _exp_as($aln, $as_res);

  #Store as patterns
  my $pattern_aln = new Bio::SimpleAlign;
  _pattern_info($aln, $pattern_aln);
  #find pred as
  my $array_ref = _add_pred_as($aln, $pattern_aln, $seq_id);
  return $array_ref;
}

=head2 _exp_as

 Title    : _exp_as
 Usage    : _exp_as($aln,  $hash_of_arrays)
 Function : Adds experimental active site data to alignment object
 Returns  : Nothing, populates the alignment object with active site residue info
 Args     : alignment object

=cut

sub _exp_as {
 
  my ($aln, $as_res) = @_;


  foreach my $seq ($aln->each_seq) {

      foreach my $pos ( @{$as_res->{$seq->id}}) {

        if($pos >= $seq->start and $pos <= $seq->end) { #Feature is in the alignment
                  
             #store column position for seq
             my $col = $aln->column_from_residue_number($seq->id, $pos);              

             #add feature to seq
             my $aa .= uc substr($seq->seq(), $col-1, 1); 

             my $feat = new Bio::SeqFeature::Generic  (  -display_name => 'experimental',
                                                         -primary => $aa,
							 -start => $col);



	     $seq->add_SeqFeature($feat);
	 }

    }
  }
}


=head2 _pattern_info

 Title    : _pattern_info
 Usage    : _pattern_info($aln_object, $aln_object)
 Function : Takes an alignment and extracts active site patterns into a second alignment
 Returns  : Nothing, populates a second alignment object with active site seqences
 Args     : alignment object, empty alignment object

=cut


sub _pattern_info {
    my ($aln, $pattern_aln) = @_;
    my (%pat_col_seq);
  
    foreach my $seq ( $aln->each_seq() ) {  

	next unless($seq->all_SeqFeatures());
           my ($pat, $col);
           foreach my $feat ( sort {$a->start <=> $b->start }  $seq->all_SeqFeatures() ) {            
              $pat .= $feat->primary_tag();   #HEK
              $col .= $feat->start() . " ";    #33 44 55
	   }

           unless(exists($pat_col_seq{"$pat:$col"})) {
	       $pattern_aln->add_seq($seq);
               $pat_col_seq{"$pat:$col"}=1;
	   }

    }
}



=head2 _add_pred_as

 Title    : _add_pred_as
 Usage    : _add_pred_as($aln_object, $aln_object)
 Function : Predicts active sites based on known active site data
 Returns  : array of active site pos
 Args     : alignment, alignment of known active sites

=cut




sub _add_pred_as {
    my ($aln, $pattern_aln, $query_seq_id) = @_;
    my $num_seq=0;
    my ($query_seq, @as_res);

    #locate query seq
    foreach my $seq ( $aln->each_seq() ) {  
	    if($seq->id eq $query_seq_id) {
		$query_seq = $seq;
		last;
	    }
    }
    die "FATAL: Can't locate query sequence [$query_seq_id] in active site alignement\n" unless($query_seq);


    my   $aligns_with = new Bio::SimpleAlign;
    foreach my $seq1 ( $pattern_aln->each_seq() ) {

   
           #See if all active site residues from seq1 exist in query seq
           my $mismatch;
           foreach my $feat ( sort {$a->start <=> $b->start }  $seq1->all_SeqFeatures() ) {

              my $aa1 = $feat->primary_tag();
              my $col = $feat->start();

              my $aa2 = uc substr($query_seq->seq, $col-1, 1);
              unless($aa1 eq $aa2) {
                  $mismatch = 1;
                  last;

              }

           }

           #Store seq1 if all active site residues are present in seq1
           unless($mismatch) {
              $aligns_with->add_seq($seq1);
           }
       }



       $num_seq = $aligns_with->no_sequences();
       return unless($num_seq);
       my (%seq_to_remove, %seq_to_rem);  #two hashes used to collect seq that need removing


        #if query seq matches more than one pattern remove subpatterns and any patterns that overlap

        #first remove sub pat
        if($num_seq>1) {
           foreach my $sequence1 ($aligns_with->each_seq() ) {
              foreach my $sequence2 ($aligns_with->each_seq() ) {

                   next if($sequence1 eq $sequence2);

                   my (%hash1, %hash2, $num_1, $num_2, %smaller, %larger);
                   #collect column positions
                   foreach my $feat1 ($sequence1->all_SeqFeatures() ) {
                       $hash1{$feat1->start} =1;
                       $num_1++;
                   }
                   foreach my $feat2 ($sequence2->all_SeqFeatures() ) {
                       $hash2{$feat2->start} =1;
                       $num_2++;
                   }


                   #see if one is a subpattern of the other
                   my $diff=0;
                   unless($num_1 eq $num_2) {

                       my $remove_seq;

                       if($num_1 > $num_2) {
                           %smaller = %hash2;
                           %larger = %hash1;
                           $remove_seq = $sequence2;

                       }
                       else {
                           %smaller = %hash1;
                           %larger = %hash2;
                           $remove_seq = $sequence1;
                       }


                       foreach my $key (keys %smaller) {
                           $diff = 1 unless(exists($larger{$key}));  #diff is true if it is not a subpattern
                       }


                       $seq_to_rem{$remove_seq}= $remove_seq unless($diff) ;
                       next unless($diff);
                   }
             }

           }
         }

         #Now remove any patterns which need removing
         foreach my $remove (keys %seq_to_rem) {
           $aligns_with->remove_seq($seq_to_rem{$remove});
         }


         unless($num_seq >=1) {
            die "FATAL: All sequences that align with active site sequences have been removed - this should never happen\n";
         }



        $num_seq = $aligns_with->no_sequences();
        #and then any patterns that overlap
        if($num_seq>1) {

           foreach my $sequence1 ($aligns_with->each_seq() ) {

              foreach my $sequence2 ($aligns_with->each_seq() ) {
                   next if($sequence1 eq $sequence2);

                   my ($seq1_st, $seq1_en, $seq2_st, $seq2_en);

                   my (%hash1, %hash2, $num_1, $num_2, %smaller, %larger);

                   #see if patterns overlap - find pattern start ends and collect column positions
                   foreach my $feat1 ($sequence1->all_SeqFeatures() ) {

                       $seq1_st = $feat1->start() if(!$seq1_st or $feat1->start() < $seq1_st);
                       $seq1_en = $feat1->start() if(!$seq1_en or $feat1->start() > $seq1_en);
                   }

                   foreach my $feat2 ($sequence2->all_SeqFeatures() ) {

                       $seq2_st = $feat2->start() if(!$seq2_st or $feat2->start() < $seq2_st);
                       $seq2_en = $feat2->start() if(!$seq2_en or $feat2->start() > $seq2_en);
                   }

                   #then see if patterns overlap - remove sequence with pattern of least identity
                   if(($seq1_st >= $seq2_st and $seq1_st <= $seq2_en) or ($seq2_st >= $seq1_st and $seq2_st <= $seq1_en)) {
                       my $remove = _identity($query_seq, $sequence1, $sequence2);
                       $seq_to_remove{$remove}= $remove;
                   }
             }

           }
         }

         #Now remove any patterns which need removing
         foreach my $remove (keys %seq_to_remove) {
           $aligns_with->remove_seq($seq_to_remove{$remove});
           $num_seq = $aligns_with->no_sequences();
           last if($num_seq eq "1"); #just in case the % identities are identical
         }


         $num_seq = $aligns_with->no_sequences();
         unless($num_seq >=1) {
            die "FATAL: All sequences that align with active site sequences have been removed - this should never happen\n";
         }



           #Add features to seq
           foreach my $sequence ($aligns_with->each_seq() ) {
                foreach my $feat ($sequence->all_SeqFeatures() ) {

                   my $actual_pos = $query_seq->location_from_column($feat->start);
                   $actual_pos = $actual_pos->start();


                   push(@as_res, $actual_pos);

 

               }
           }
           return \@as_res

}


=head2 _identity

 Title    : _identity
 Usage    : _identity($sequence1 , $sequence2, $sequence3)
 Function : Identifies seq with lowest % identity to sequence1
 Returns  : The sequence which has the lowest % id to sequence 1
 Args     : sequence1, sequence2, sequence3.

=cut


sub _identity {
    my $seq1 = shift;
    my @aligns_with = @_;
          my $lower_identity=100;
          my $lower_identity_seq;
          foreach my $s (@aligns_with) {
             my $tmp_aln = new Bio::SimpleAlign;
             $tmp_aln->add_seq($s);
             $tmp_aln->add_seq($seq1);

             my $identity = $tmp_aln->percentage_identity();
             if($identity < $lower_identity) {
                 $lower_identity = $identity;
                 $lower_identity_seq = $s;
             }

          }
          return $lower_identity_seq;
}

=head1 COPYRIGHT

Copyright (c) 2007: Genome Research Ltd.

Authors: Rob Finn (rdf@sanger.ac.uk), John Tate (jt6@sanger.ac.uk),
         Jaina Mistry (jm14@sanger.ac.uk)

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
