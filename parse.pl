#! usr/bin/perl

##### Assumptions: 
# 1) The gaps in reference (i.e. deletion is reference or insertion in reads) also gets index
# for example,
# index:        1 2 3 4 5 6 7 8 9 10
# reference:    A T G G - - T A T G
# read(query):  A T - C T T T A T G (cigar: 2M1D1M2I5M)


# Input from command line:
# 1) the sam file containing read information (must be without header)
# 2) The position of interest (the nucleotide position where the nucleotide frequencies will be calculated)
# example: perl parse.pl chr1_small.sam 36932148



#Taking input from command line
@file = `cat $ARGV[0]`;

#set the concerened position in to a variable 
$pos = $ARGV[1];

#setting the initial counts of each nucleotide to 0
$A=0;$T=0;$G=0;$C=0;

# parsing each line of the input file
foreach $line(@file){
	chomp($line);

	# reading each fields of a line of input file into an array (@fields)
        @fields=split("\t",$line);

	# $diff = distance between the start site of each read and the concerened position ($pos) 
	$diff = $pos-$fields[3];

	# reading the cigar string field into $cigar 
	$cigar = $fields[5];

	@typearray = ();

	# The following if condition:
	# 1) discards reads with missing cigar string information
	# 2) discards reads starting after $pos (as they can't overlap with $pos)

	if($cigar != "*" && $diff >= 0){
		
		###### Section 1: subsetting the cigar string information into @match #######
		# for example, 5M3I6M will be: $match[0] = 5M, $match[1] = 3I and $match[2] = 6M
		@match = ($cigar =~ /\d+[M|D|I|S|H|N|P]/g);
		@match= grep { $_ ne '' } @match;
		############################# END of section 1 ##############################
		


		###### section 2: conversion of cigar string into actual string #############
		# for example 2M3D4I ==> MMDDDIIII
		# the converted information is stored in @typearray
		@value=();
		foreach (@match){
			@value = ($_ =~ /(\d+)([a-z])/gi);
			for ($i=0;$i<$value[0];$i++){
				push(@typearray,$value[1]);
			}
		}
		############################# END of section 2 ###############################

		
		
		##### section 3: modification of read sequence according to cigar information stored in @typearray #####
		# A "*" is inserted into the palces of clipping (H and S), deletion (D) and padding (P)
		# However, in contrast to other 3 cases (H,D,P), no shifting of array is done in case of S.
		# example1: ATCCGCATGGATCGTGAC with 3H5M3D13M is converted to --> ***ATCCG***CATGGATCGTGAC
		# example2: ATCCGCATGGATCGTGAC with 3S5M3D13M is converted to --> ***CGCAT***GGATCGTGAC

		@read = split("",$fields[9]);
		#@actualread = @read;	
		foreach($i=0;$i<scalar(@typearray);$i++){
			if($typearray[$i]=~/H|D|P/i){
				splice @read, $i,0,'*';
			}
			elsif($typearray[$i]=~/S/i){
				$read[$i] = "*";
			}
		}
		############################# END of section 4##########################################################
		
		
		##### section 4: calculation of nucleotide frequency from modified read (@read) ########################
		# "*" representing H,S,D,P are not counted
		# $read[$diff] refers to the nucleotide at the position of interest (stored in $pos)
		if(exists$read[$diff]){
			if($read[$diff]=~/A/){
				$A++;
			}
			if($read[$diff]=~/T/){
				$T++;
                        }
			if($read[$diff]=~/G/){
                                $G++;
                        }
			if($read[$diff]=~/C/){
                                $C++;
                        }
		}
		############################# END of section 4 ##########################################################
	}
}


print "Reference allele (T) count = $T, Alternate allele (G) count = $G\n";
