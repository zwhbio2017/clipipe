#!/home/CLIPipe_user/bin/miniconda3/envs/clipipe/bin/perl

#
# PARCLIP suite, v1.0
# Pavel Morozov, RNA Biology laboratory, HHMI/Rockefeller University, 2016
# 

# import file with global parameters
require "parameters";

# check if there are enough command line parameters
$n=@ARGV;
if($n<2)
	{
	print "prepare_fasta.pl, PAR-CLIP suite v1.0\nUsage:\n\tprepare_fasta.pl input_fastq[.gz] output_prefix output_directory\n\n";
	exit(0);
	}

$output_dir=$ARGV[2];

# check if output directory exists and if it can be created
if(!-e $output_dir)
	{
	system("mkdir $output_dir");
	if(!-e $output_dir)
		{
		print "Output directory does not exists and can not be created\n";
		exit(1);
		}
	}


# check if input fastq file is zipped
$fastqfile=$ARGV[0];
if(index($ARGV[0],".gz")==length($ARGV[0])-3)
	{
	$fastqfile=$ARGV[0]; $fastqfile=~s/\.gz//;
	$cmd="gunzip -c $ARGV[0] > $output_dir/$ARGV[1].fastq";
	print "Unzipping input file:\n$cmd\n";
	system($cmd);
	$fastqfile="$output_dir/$ARGV[1].fastq";
	}

# prepare file names
$fasta=$output_dir."/".$ARGV[1].".fasta";
# $trimmed=$output_dir."/".$ARGV[1].".cfasta";
$trimmed=$output_dir."/".$ARGV[1].".fasta";
$collapsed=$output_dir."/".$ARGV[1].".fau";

# fastq to fasta
$cmd="$script_directory/fastq2fasta.pl $fastqfile > $fasta";
print "Converting FASTQ to FASTA:\n$cmd\n";
system($cmd);

# run cutadapt
if($use_cutadapt eq "Yes" || $use_cutadapt eq "yes" || $use_cutadapt eq "YES")
	{
        $cmd="$cutadapt $fasta $cutadapt_parameters > $trimmed";
	print "Converting FASTQ to FASTA:\n$cmd\n";
	system($cmd);
        }

# collapse
$cmd="$script_directory/fasta2unique.pl $trimmed > $collapsed";
print "Collapsing FASTA:\n$cmd\n";
system($cmd);

# cleanup
print "Cleaning temporary files...";
system("rm $fasta");
system("rm $trimmed");
print "\nDone. Use file $collapsed for further analysis\n\n";
