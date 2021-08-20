#!/home/CLIPipe_user/bin/miniconda3/envs/clipipe/bin/perl

#
# PARCLIP suite, v1.0
# Pavel Morozov, RNA Biology laboratory, HHMI/Rockefeller University, 2016
# 

#
# Filters out fasta file by substring of nucleotides from 5' end 
#

$min=$ARGV[1];
$max=$ARGV[2];

open(IN,$ARGV[0]) || die();
while(<IN>)
        {
        if(/^\>/)
                {
                chomp();
		$name=$_;
		}
	else
		{
		if(index($_,$ARGV[3])==0 && index($_,"N")==-1)
                        {
                        $l=length($_)-1;
                        if($l>=$min && $l<=$max)
                                {
				$seq{$name}=$_;
                                }
                        }
                }
        }
close(IN);

while(($k,$v)=each(%seq)) 
	{
        print "$k\n$v";
	}
