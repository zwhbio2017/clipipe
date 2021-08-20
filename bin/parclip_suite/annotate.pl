#!/home/CLIPipe_user/bin/miniconda3/envs/clipipe/bin/perl

#
# PARCLIP suite, v1.0
# Pavel Morozov, RNA Biology laboratory, HHMI/Rockefeller University, 2016
# 

# General annotation by groups with coding for all possible substitutions (currently no indels).

# import file with global parameters
require "parameters";

##########################################################
# Read the configuration file and initialization

@db=();		# list of databases
$ndb=0;		# number of databases

$work_directory = $ARGV[0];
$fasta_file = $ARGV[1];
$conf_file = $ARGV[2];
$nick_name = $ARGV[3];

read_configuration($conf_file);

# Types of mismatches ####################################

# 0 for none
$imism{"A"}{"T"}=1; # A>T
$imism{"A"}{"C"}=2; # A>C
$imism{"A"}{"G"}=3; # A>G
$imism{"T"}{"A"}=4; # T>A
$imism{"T"}{"C"}=5; # T>C
$imism{"T"}{"G"}=6; # T>G
$imism{"C"}{"A"}=7; # C>A
$imism{"C"}{"T"}=8; # C>T
$imism{"C"}{"G"}=9; # C>G
$imism{"G"}{"A"}=10; # G>A
$imism{"G"}{"T"}=11; # G>T
$imism{"G"}{"C"}=12; # G>C
# 13 for Ins / for bowtie2
# 14 for Del / for bowtie2

$inv{"A"}="T";
$inv{"T"}="A";
$inv{"C"}="G";
$inv{"G"}="C";

# init variables
# $work_directory\/
$out_bwt="$nick_name.anno_bwt";
# $work_directory\/
$out_bwt_temp="$nick_name.anno_bwt_t";
# $work_directory\/
$out_anno="$nick_name.anno";

##########################################################
if(-e $out_bwt)
        {
        $cmd="rm $out_bwt";
        print "$cmd\n";
        system($cmd);
        }

if(-e $out_anno)
	{
	$cmd="rm $out_anno";
	print "$cmd\n";
	system($cmd);
	}

#########################################################
# run the reads fasta file against all annotation library

for($z=0;$z<$ndb;$z++)
	{
	$cmd="$bowtie --sam -v 2 -p 8 -f -k 20 $db[$z] $fasta_file > $out_bwt_temp";
	print "Executing:\n\t$cmd\n";
	system($cmd);
	system("$script_directory/compress_bwt.pl $out_bwt_temp >> $out_bwt");
	system("rm $out_bwt_temp");
	}

print "Start parsing\n";

# now parsing the bowtie output
open(IN,$out_bwt) || die("Can not open bowtie output $out_bwt\n");
while(<IN>)
	{
	if(/^\@/){next;}
	@s=split();
	if($s[1]!=0 && $s[1]!=16)
		{next;}
	($nam,$cnt)=split(/\|/,$s[0]);

	# categoriy and name
	($ct,$cfnm)=split(/\|/,$s[2]);

	# determine if it is the right strand for this annotation type
	if($strand{$ct} eq "p" && $s[1] == 16){next;}
        if($strand{$ct} eq "m" && $s[1] == 0){next;}
	
	# number of mismatches and type of mismatches
	@w=split(/\:/,$s[11]);

# seq7|1	0	8|SNORD4B_chr17_27050598-27050872	106	255	32M	*	0	0	CCAAATGATGCATATGTTAGCGACCAAAGCCT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	XA:i:1	MD:Z:0G31	NM:i:1
# seq14|1	0	8|SNORD88C_chr19_51305481-51305778	107	255	27M	*	0	0	CTCCCATGATGTCTAGCACTGGGCACT	IIIIIIIIIIIIIIIIIIIIIIIIIII	XA:i:2	MD:Z:13C10T2	NM:i:2

	$mut=$w[2];
	if($mut>0)
		{
		@q=split(/\:/,$s[12]);
		@bytes=split('',$s[9]);
		$kk=$q[2];
		$kk=~s/[0-9]//g;
		@syms=split('',$kk);
		@nums=split(/[A-z]/,$q[2]);

		$pos=0;
		$code=0;

		for($i=0;$i<$mut;$i++)
			{
	                if($s[1] == 16)
        	                {
                	        $first=$inv{$syms[$i]};
                        	$second=$inv{$bytes[$pos+$nums[$i]]};
	                        }
        	        else
                	        {
                        	$first=$syms[$i];
	                        $second=$bytes[$pos+$nums[$i]];
        	                }

			if($i>0)
				{
				$code+=(($i*16)*$imism{$first}{$second});
				}
			else
				{
                                $code+=($imism{$first}{$second});
				}

			$pos+=($nums[$i]+1);
			}
		$mut=$code;
		}

	$rlen{$s[0]}=length($s[9]);

	if(exists $anno{$s[0]})
		{
		$anno{$s[0]}.="$ct,$mut,";
		}
	else
		{$anno{$s[0]}="$ct,$mut,";}
	}
close(IN);

#### Now output (no sorting)
#

open(OUT,">".$out_anno) || die();

# First count stats according to hierarchy
$ge_hits=0;
$tr_hits=0;
while(($k,$v)=each(%anno))
	{
	print OUT "$k\t$rlen{$k}\t$v\n";
	}
close(OUT);

# Zipping the bwt file ################
system("rm $out_bwt.gz");
system("gzip $out_bwt");

#### Subroutines #######################
#

#-----------------------#
sub read_configuration
{
$fil=$_[0];

open(IN,$fil) || die("Can not open configuration file $fil\n");
while(<IN>)
	{
	chomp();
	@s=split();
	if($s[0] eq "annotation_db")
		{
		$db[$ndb]=$s[1];
		$ndb++;
		}
        if($s[0] eq "hier")
		{
		@hierarchy=split();
		$n=@hierarchy-1;
	for($i=1;$i<=$n;$i++)
        		{
        		$hier{$hierarchy[$i]}=$i-1;
      	  		}
		$ncat=@n;
		while(<IN>)
        		{

		        # sample string: category_number        string_preference       category_name
		        # 1     p       category name
		        # 2     m       category name
		        # 3     b       category name

		        # the classes mRNA and genome are a special classes and are treated differently in other subroutines

		        chomp();
		        @s=split(/\t/);
		        $catnames{$s[0]}=$s[2];
		        $strand{$s[0]}=$s[1];
		        $catid{$s[2]}=$s[0];
		        }
		}
	}
close(IN);
}

#--------------------#
sub nmut
{
$mm=2; 	# by default 2 mismatches
if($_[0]==0){$mm=0;} # if perfect match
else    
	{
        if($_[0]<16)
        	{
	        $mm=1;
		}
	}
return $mm;
}
