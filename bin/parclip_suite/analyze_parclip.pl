#!/home/CLIPipe_user/bin/miniconda3/envs/clipipe/bin/perl

#
# PARCLIP suite, v1.0
# Pavel Morozov, RNA Biology laboratory, HHMI/Rockefeller University, 2016
# 

# import file with global parameters
require "parameters";

$n=@ARGV;
if($n<5)
	{
	print "\n\nUsage:\n./analyze_parclip.pl input_fasta prefix minimal_length maximal_length configuration_file output_directory\n\n\n";
	exit(0);
	}

$input=$ARGV[0];
$prefix=$ARGV[1];
$minlen=$ARGV[2];
$maxlen=$ARGV[3];
$config=$ARGV[4];
$outdir=$ARGV[5];

# check if output directory exists and if it can be created
if(!-e $outdir)
        {
        system("mkdir $outdir");
        if(!-e $outdir)
                {
                print "Output directory does not exists and can not be created\n";
                exit(1);
                }
        }

@bin=("AA","AT","AG","AC", "TA","TT","TG","TC", "CA","CT","CG","CC", "GA","GT","GG","GC");
$numbins=16;

$cmd="";

# split to fasta files by letters
for($i=0;$i<$numbins;$i++)
	{
	runcmd("$script_directory/fasta_split_by_letter.pl $input $minlen $maxlen $bin[$i] > $outdir/$prefix\_$bin[$i].fa");
	}

# make a list of fasta files
runcmd("ls $outdir/$prefix\_??.fa > $outdir/$prefix\_fa.lst");

# now run annotate_2 on all collapsed fasta file
for($i=0;$i<$numbins;$i++)
        {
        runcmd("$script_directory/annotate.pl ./ $outdir/$prefix\_$bin[$i].fa $config $outdir/$prefix\_$bin[$i]\_res");
        }

# make a list of anno files
runcmd("ls $outdir/$prefix\_??\_res.anno > $outdir/$prefix\_anno.lst");

# run an_anno_pos.pl on the set of annotations
runcmd("$script_directory/an_anno_lst.pl $outdir/$prefix\_anno.lst $outdir/$prefix\_fa.lst $config $outdir/$prefix $minlen $maxlen");

###############################################
sub runcmd
{
if(length($_[0])>3)
	{
	print "Executing: $_[0]\n";
	system($_[0]);
	}
}
