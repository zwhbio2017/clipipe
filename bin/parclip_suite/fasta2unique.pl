#!/home/CLIPipe_user/bin/miniconda3/envs/clipipe/bin/perl

#
# PARCLIP suite, v1.0
# Pavel Morozov, RNA Biology laboratory, HHMI/Rockefeller University, 2016
# 

# fasta to collapsed (unioque) fasta
# variant for large datasets, collapsing go by subset of reads by starting nucleotides. 


@bin=("AA","AT","AG","AC", "TA","TT","TG","TC", "CA","CT","CG","CC", "GA","GT","GG","GC");
$numbins=16;

$n_out=1;
for($i=0;$i<$numbins;$i++)
	{
	%cnt=();
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
			if(index($_,$bin[$i])==0)
				{
				$cnt{$_}++;
				}
			}
		}
	close(IN);
	
	while(($k,$v)=each(%cnt))
		{
		print ">seq$n_out\|$v\n$k\n";
		$n_out++;
		}
	}
