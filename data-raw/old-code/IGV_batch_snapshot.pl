#!/usr/bin/perl -w

##This programs reads in the input file AA-SNV-INPUT.tsv  ---- NOTE: YOU WILL NEED TO ADD a row "0 \t 0 \t 0 \t 0 \t 0 \t END" to the end of the file as the program needs it to correctly operate
##Input file contains SNVs from multiple samples
##Randomly selects 0-100% of the SNVs from each sample
##outputs the selected SNV coordinates in a suitable format for IGV batch snapshots
##
##Note: the folder/file paths are okay but very clunky ... be careful with line 69 and 101 to correctly name the output pdf file.

# For documentaiton on IGV batch snapshots see:
# https://software.broadinstitute.org/software/igv/batch
# https://github.com/igvteam/igv/wiki/Batch-commands


use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use POSIX qw(ceil);

$input_loc = "/home/gmsywss/Arnoud_IGVs/37N-SNV-INPUT.csv";

$working_dir = "/home/gmsywss/Arnoud_IGVs/";
$script_loc = join "", $working_dir, "batch_script_folder/SNV/";
$script_suffix = "_SNV_IGV-batch.txt";

$snapshot_loc = join "", $working_dir, "SNV_snapshot_folder/";
$Tbam_loc = "/gpfs/nextgen1/gmsarno/CCA_Taiwan/37/";
$Nbam_loc = "/gpfs/nextgen1/gmsarno/CCA_Taiwan/37/";
$Tbam_suffix = "T.rmdup.bam";
$Nbam_suffix = "N.rmdup.bam";

$igv_prog = "/home/gmsywss/IGV_2.4.13/igv.sh";
$pdf_loc = join "", $working_dir, "output_pdf_folder/SNV/";

open(INPUT_LOC, $input_loc) or die("Could not find/open the input file.\n"); #just checking if file can be opened

$sample_percent = 1;
$sample_hold = 0;
$first_line = 1;
foreach $line(<INPUT_LOC>) {
 	chomp($line);
	print $line, "\n";
	if ( $first_line == 0 )	{ #from line 2 of input file onwards
		@temp_tab = split( /\t/, $line );
		if ( $temp_tab[5] eq $sample_hold ) { #we are still on the same sample
			push( @list_chr, $temp_tab[0] );
			push( @list_pos, $temp_tab[1] );
			$sample_hold = $temp_tab[5];
		} elsif ( $temp_tab[5] eq "END" ) {
			#sampling the list of variants
			if ( $sample_percent == 1) { #we want ALL entries
				$grab_size = $#list_chr + 1;
				@shuffle_deck = ( 0 .. $#list_chr );
			} else { #we want to sample from list
				@shuffle_deck = shuffle( 0 .. $#list_chr );
				$grab_size = ceil( $sample_percent * ($#list_chr + 1) );
			}
			print "Sample Name:  ", $sample_hold, "\tSNV count:  ", $#list_chr+1, "\t\tSample count:  ", $grab_size, "\n\n";

			#Creating batch script file
			$output_filename = join "", $script_loc, $sample_hold, $script_suffix;
			open ( FILEHANDLE, ">", "$output_filename" ) or die ( "Cannot create IGV batch file.\n" ); #open output file
			print FILEHANDLE "genome hg19\nsnapshotDirectory ", $snapshot_loc, $sample_hold, "_snapshot_folder/\nnew\nload ", $Tbam_loc, $sample_hold, $Tbam_suffix ,"\nload ", $Nbam_loc, $sample_hold, $Nbam_suffix,"\n\nmaxPanelHeight 500\n";
			for ( $i = 0; $i < $grab_size; $i++ ) {
				print FILEHANDLE "goto ", $list_chr[$shuffle_deck[$i]],":", $list_pos[$shuffle_deck[$i]],"\nsort base\nsquish\nsnapshot ", $list_chr[$shuffle_deck[$i]] ,"-", $list_pos[$shuffle_deck[$i]],"_base.png\n";
			}
			print FILEHANDLE "exit\n";
			close(FILEHANDLE);

			#Running IGV to produced snapshots and converting to pdf format
			@igv_call = ($igv_prog, "-b", $output_filename);
			print "\nCommand call: @igv_call\n\n";
			system(@igv_call);
			$png_files = join "", $snapshot_loc, $sample_hold, "_snapshot_folder/*png";
			$pdf_output = join "", $pdf_loc, $sample_hold, "_dinuc_IGV.pdf";
			@conversion_call = ( "convert", $png_files, $pdf_output );
			print "\n\n Command call: @conversion_call\n\n";
			system(@conversion_call);
		} else { #need sample 1% of SNV entries, create output batch file, fill with correct commands, end old output file, reset chr & pos array, load in new array entry for new sample

			#sampling the list of variants
			if ( $sample_percent == 1) { #we want ALL entries
				$grab_size = $#list_chr + 1;
				@shuffle_deck = ( 0 .. $#list_chr );
			} else { #we want to sample from list
				@shuffle_deck = shuffle( 0 .. $#list_chr );
				$grab_size = ceil( $sample_percent * ($#list_chr + 1) );
			}
			print "Sample Name:  ", $sample_hold, "\tSNV count:  ", $#list_chr+1, "\t\tSample count:  ", $grab_size, "\n\n";

			#Creating batch script file
			$output_filename = join "", $script_loc, $sample_hold, $script_suffix;
			open ( FILEHANDLE, ">", "$output_filename" ) or die ( "Cannot create IGV batch file.\n" ); #open output file
			print FILEHANDLE "genome hg19\nsnapshotDirectory ", $snapshot_loc, $sample_hold, "_snapshot_folder/\nnew\nload ", $Tbam_loc, $sample_hold, $Tbam_suffix ,"\nload ", $Tbam_loc, $sample_hold, $Nbam_suffix,"\n\nmaxPanelHeight 500\n";
			for ( $i = 0; $i < $grab_size; $i++ ) {
				print FILEHANDLE "goto ", $list_chr[$shuffle_deck[$i]],":", $list_pos[$shuffle_deck[$i]],"\nsort base\nsquish\nsnapshot ", $list_chr[$shuffle_deck[$i]] ,"-", $list_pos[$shuffle_deck[$i]],"_base.png\n";
			}
			print FILEHANDLE "exit\n";
			close(FILEHANDLE);

			#Running IGV to produced snapshots and converting to pdf format
			@igv_call = ($igv_prog, "-b", $output_filename);
			print "\nCommand call: @igv_call\n\n";
			system(@igv_call);
			$png_files = join "", $snapshot_loc, $sample_hold, "_snapshot_folder/*png";
			$pdf_output = join "", $pdf_loc, $sample_hold, "_dinuc_IGV.pdf";
			@conversion_call = ( "convert", $png_files, $pdf_output );
			print "\n\n Command call: @conversion_call\n\n";
			system(@conversion_call);

			#resetting the chr and pos arrays and start on next sample
			@list_chr = ();
			@list_pos = ();
			@temp_tab = split( /\t/, $line );
			push( @list_chr, $temp_tab[0] );
			push( @list_pos, $temp_tab[1] );
			$sample_hold = $temp_tab[5];
		}
	} elsif ( $first_line == 1 ) { #initialize first line of the input file
		@temp_tab = split( /\t/, $line );
		push( @list_chr, $temp_tab[0] );
		push( @list_pos, $temp_tab[1] );
		$sample_hold = $temp_tab[5];
		$first_line = 0;
		print $sample_hold, "\n";
		print "stuck here???\n";
	}

}
close(INPUT_LOC);
