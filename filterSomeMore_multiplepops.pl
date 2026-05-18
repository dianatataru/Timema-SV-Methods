#!/usr/bin/perl


# filter vcf files based on coverage


open(IN,"KeepSNPs.txt") or die "failed initial read\n";
while(<IN>){
	chomp;
	push(@keep,$_);
}
close(IN);


foreach $in (@ARGV){
	open (IN, $in) or die "Could not read the infile = $in\n";
#	$in =~ m/^([a-zA-Z0-9_]+\.vcf)$/ or die "Failed to match the variant file\n";
	open (OUT, "> morefilter_2x_FHA_all_oneref_v2.vcf") or die "Could not write the outfile\n";


	while (<IN>){
		chomp;
		if (m/^\#/){ ## header row, always write
			$flag = 1;
		}
		elsif (m/^\S+/){ ## this is a sequence line, you migh need to edit this reg. expr.
			$flag = shift(@keep);
			if ($flag == 1){
				$cnt++; ## this is a good SNV
			}
		}
		else{
			print "Warning, failed to match the chromosome or scaffold name regular expression for this line\n$_\n";
			$flag = 0;
		}
		if ($flag == 1){
			print OUT "$_\n";
		}
	}
	close (IN);
	close (OUT);

	print "Finished filtering $in\nRetained $cnt variable loci\n";
}
