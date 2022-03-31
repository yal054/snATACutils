#!/usr/bin/perl -w

use strict;

#####
# Get highest and mean phastCons score for each binding region in BED file
# 10 Sept 2010.  Ian Donaldson.
# 20 Oct 2010 modified to read BED format
#####

# Command line usage

if(@ARGV!=2) {
   die ("USAGE: $0 | Regions BED file | Output\n");
} 

open(BED, "<$ARGV[0]") or die("Could not open input BED file!\n"); 
open(OUTPUT, ">$ARGV[1]") or die("Could not open output file!\n");

# Work thru each line of BED file
while(defined(my $line = <BED>)) {
   chomp($line);

   my ($chr, $start, $end) = '';

   # Skip line if line does not start with 'chr'
   unless($line =~ /^chr/) { next }

   # Remove EOL from each line
   chomp($line);

   # Split each line by TAB
   my @line_bits = split(/\t/, $line);

   $chr = $line_bits[0];
   $start = $line_bits[1];
   $end = $line_bits[2];

   # Run bigWigSummary 
   system("~/apps/kentUtils/bigWigSummary -type=max /projects/ps-renlab/yangli/genome/mm10/mm10.60way.phastCons/mm10.60way.phastCons.bw $chr $start $end 1 > $$.max_tmp 2>&1");
   system("~/apps/kentUtils/bigWigSummary -type=mean /projects/ps-renlab/yangli/genome/mm10/mm10.60way.phastCons/mm10.60way.phastCons.bw $chr $start $end 1 > $$.mean_tmp 2>&1");

   # open temp file
   open(MAX_TEMP, "<$$.max_tmp") or die("Cannot open max temp file!\n");
   open(MEAN_TEMP, "<$$.mean_tmp") or die("Cannot open mean temp file!\n");


   # MAX phastCons
   my $max_temp_output = '';

   while(<MAX_TEMP>) {
      $max_temp_output .= $_;
      chomp($max_temp_output);
   }

   # close max temp file
   close(MAX_TEMP);

   # remove blank lines s/^$//g
   if($max_temp_output =~ /^no/) { 
      $max_temp_output = '';
   }


   # MEAN phastCons
   my $mean_temp_output = '';

   while(<MEAN_TEMP>) {
      $mean_temp_output .= $_;
      chomp($mean_temp_output);
   }

   # close max temp file
   close(MEAN_TEMP);

   # remove blank lines s/^$//g
   if($mean_temp_output =~ /^no/) { 
      $mean_temp_output = '';
   }


   # save temp file variable to $ARGV[1]
   print OUTPUT "$line_bits[0]\t$line_bits[1]\t$line_bits[2]\t$max_temp_output\t$mean_temp_output\n";
}

# remove temporary file to final output
system("rm $$.max_tmp");
system("rm $$.mean_tmp");

# close files
close(BED);
close(OUTPUT);

exit;


