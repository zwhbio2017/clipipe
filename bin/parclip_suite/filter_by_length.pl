#!/home/CLIPipe_user/bin/miniconda3/envs/clipipe/bin/perl

#
# PARCLIP suite, v1.0
# Pavel Morozov, RNA Biology laboratory, HHMI/Rockefeller University, 2016
# 

# import file with global parameters
require "parameters";

$n=@ARGV;
if($n<3)
	{
	print "\n\nUsage:\n./filter_by_length.pl input_fasta minimal_length maximal_length\n\n\n";
	exit(0);
	}

open(IN,$ARGV[0]) || die();
while(<IN>)
	{
	chomp();
	if(/^\>/)
		{
		$name=$_;
		}
	else
		{
		$l=length($_);
		if($l>=$ARGV[1] && $l<=$ARGV[2])
			{
			print "$name\n$_\n";
			}
		}
	}	
