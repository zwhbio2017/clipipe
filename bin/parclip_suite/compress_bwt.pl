#!/home/CLIPipe_user/bin/miniconda3/envs/clipipe/bin/perl

#
# PARCLIP suite, v1.0
# Pavel Morozov, RNA Biology laboratory, HHMI/Rockefeller University, 2016
# 

open(IN,$ARGV[0]) || die();
while(<IN>)
	{
	if(/^\@/){next;}
	@s=split();
	if($s[1] == 0 || $s[1] == 16)
		{
		print $_;
		}
	}
close(IN);
