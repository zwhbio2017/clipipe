#!/home/CLIPipe_user/bin/miniconda3/envs/clipipe/bin/perl

#
# PARCLIP suite, v1.0
# Pavel Morozov, RNA Biology laboratory, HHMI/Rockefeller University, 2016
# 


#
# Version of annotation by groups for PARCLIP data, i.e. with special attention to the T>C transition
# Annotation for each length is separate.
# Uses list of annotation files for parallel processing.
# Uses list of fasta files to get statistics before annotation  
#

########################################################
# Read the configuration file and initialization

# import file with global parameters
require "parameters";

$maxlen=150;    # maximal length of read

$db=""; 	# anotation database, bowtie formatted
$hier="";	# hierarchy file
$gnm="";	# genome, bowtie formatted
$trans="";	# transcriptome, optional

read_configuration($ARGV[2]);
# init_hierarchy($hier);

# read initial reads statistics
open(INF,$ARGV[1]) || die();
while(<INF>)
        {
        chomp();
        open(IN,$_) || die("Can't open $_\n");
        while(<IN>)
                {
                chomp();
		if(/^\>/)
			{
			($nam,$cnt)=split(/[\|\-]/);
			}
		else
			{
			$l=length($_);
			if(exists $rl_count{$l})
				{$rl_count{$l}+=$cnt;}
			else 
				{$rl_count{$l}=$cnt;}
			}
		}
	close(IN);
	}
close(INF);

# Types of mismatches ##################################

# 0 for none
$mism{1}="AT"; # A>T
$mism{2}="AC"; # A>C
$mism{3}="AG"; # A>G ###
$mism{4}="TA"; # T>A
$mism{5}="TC"; # T>C ##########
$mism{6}="TG"; # T>G
$mism{7}="CA"; # C>A
$mism{8}="CT"; # C>T
$mism{9}="CG"; # C>G
$mism{10}="GA"; # G>A
$mism{11}="GT"; # G>T
$mism{12}="GC"; # G>C

$maxlen=$ARGV[5]+1;
$minlen=$ARGV[4];

# 13 for Ins / for bowtie2
# 14 for Del / for bowtie2

for($i=0;$i<$maxlen+1;$i++)
	{
	for($j=0;$j<15;$j++)
		{
		for($k=0;$k<15;$k++)
			{
			for($h=3;$h<100;$h++)
				{
				$mut_list[$h][$i][$j][$k]=0;
				}
			}
		}
	}

open(OUT,">".$ARGV[3]) || die();
open(OUT_S,">".$ARGV[3].".ctab") || die();
open(INF,$ARGV[0]) || die();

while(<INF>)
	{
	chomp();
	open(IN,$_) || die("Can't open $_\n");
	# Variant with traditional counting, just with special attention to T>C
	#
	while(<IN>)
		{
		chomp();
		($namm,$rsize,$data)=split(/\t/);
		($namme,$cnt)=split(/[\|\-]/,$namm);

		# first what is the smallest number of mutations
		$mut=3;
		@s=split(/\,/,$data);
		$n=@s;

		# First get the minimal number of mustations
		for($i=0;$i<$n;$i+=2)
			{
			# determine number of mutations
			$mn=nmut($s[$i+1]);
			if($mn<$mut){$mut=$mn;}
			}

	        %check=();
        	# now go through all categories and check if it is unique mapper, i.e. no more than one hit per category
		$flag=0;
        	for($i=0;$i<$n;$i+=2)
                	{
			if(nmut($s[$i+1])>$mut)
				{next;}
			$check{$s[$i]}++;
			if($check{$s[$i]}>1)
				{
				$flag=1;
				last;
				}
			}
		
		# now the lowerst class in hierarchy with this number of mutation
		$ind=100000;
		for($i=0;$i<$n;$i+=2)
                	{
	                if(nmut($s[$i+1])==$mut)
				{
				if($hier{$s[$i]}<$ind)
					{
					$ind=$hier{$s[$i]};
					$h_ind=$s[$i];
					$h_byte=$s[$i+1];			
					}
				}
	               }

		$last=int($h_byte/16);
	        $first=$h_byte-$last*16;

                print OUT_S "$namm\t$h_ind\t$first\t$last\n";

		$mut_list[$rsize][$h_ind][$first][$last]+=$cnt;
		}
	close(IN);
	}

# print the results, simplified

###################### forming title string ##################################
$grand_total=0;

print OUT "mapping summary\nread length\ttotal before annotation\tnon nannotated\ttotal annotated\td0\td1\td2\t";
for($i=0;$i<$ncat+100;$i++)
	{
	if(!exists $catnames{$i}){next;}
	print OUT $catnames{$i}."\t";
	}
print OUT "\n";

for($h=$minlen;$h<$maxlen;$h++) # for each length
        {
        $total=0;
        $totald0=0;
        $totald1=0;
        $totald2=0;

        for($i=1;$i<$ncat+100;$i++) # for each category
                {
                if(!exists $catnames{$i}){next;}

                # calculating totals
                $all1=0;
                $all2=0;
                for($j=1;$j<13;$j++)
                        {
                        $all1+=$mut_list[$h][$i][$j][0];
                        for($k=1;$k<13;$k++)
                                {
                                $all2+=$mut_list[$h][$i][$j][$k];
                                }
                        }
                $total+=($mut_list[$h][$i][0][0]+$all1+$all2);
                $totald0+=$mut_list[$h][$i][0][0];
                $totald1+=$all1;
                $totald2+=$all2;
                }
	$grand_total+=$total;
	$un=0;
	$total_before=0;
	
	if(exists $rl_count{$h})
		{
		$total_before=$rl_count{$h};
		$un=$total_before-$total;
		}

        print OUT "$h\t$total_before\t$un\t$total\t$totald0\t$totald1\t$totald2\t";
	
	
	# now by categories
        for($i=1;$i<$ncat+100;$i++) # for each category
                {
		$total=$mut_list[$h][$i][0][0];
                if(!exists $catnames{$i}){next;}

                # calculating totals
                for($j=1;$j<13;$j++)
                        {
                        $total+=$mut_list[$h][$i][$j][0];
                        for($k=1;$k<13;$k++)
                                {
                                $total+=$mut_list[$h][$i][$j][$k];
                                }
                        }
		print OUT "$total\t";
                }
	print OUT "\n";
	}

# output of perfect matches
print OUT "\n\nPerfect match\nRead_length\t";
for($i=0;$i<$ncat+100;$i++)
        {
        if(!exists $catnames{$i}){next;}
        print OUT $catnames{$i}."\t";
        }
print OUT "\n";
for($h=$minlen;$h<$maxlen;$h++) # for each length
        {
	print OUT "$h\t";
        for($i=1;$i<$ncat+100;$i++) # for each category
                {
                if(!exists $catnames{$i}){next;}
                print OUT "$mut_list[$h][$i][0][0]\t";
                }
	print OUT "\n";
	}

# output of one mismatch
print OUT "\n\nOne mismatch\nRead_length\t";
for($i=0;$i<$ncat+100;$i++)
        {
        if(!exists $catnames{$i}){next;}
        print OUT $catnames{$i}."\t";
        }
print OUT "\n";
        
for($h=$minlen;$h<$maxlen;$h++) # for each length
        {
        print OUT "$h\t";
        for($i=1;$i<$ncat+100;$i++) # for each category
                {
		$total=0;
                if(!exists $catnames{$i}){next;}
                for($j=1;$j<13;$j++)
                        {
                        $total+=$mut_list[$h][$i][$j][0];
                        }
                print OUT "$total\t";
                }
        print OUT "\n";
        }

# output of T>C count for one mismatch reads
print OUT "\n\nT>C transition count for one mismatch\nRead_length\t";
for($i=0;$i<$ncat+100;$i++)
        {
        if(!exists $catnames{$i}){next;}
        print OUT $catnames{$i}."\t";
        }
print OUT "\n";
        
for($h=$minlen;$h<$maxlen;$h++) # for each length
        {
        print OUT "$h\t";
        for($i=1;$i<$ncat+100;$i++) # for each category
                {
                $total=0;
                if(!exists $catnames{$i}){next;}
                print OUT "$mut_list[$h][$i][5][0]\t";
                }
        print OUT "\n";
        }


# output of additional tables for categories with over 15% annotation
for($i=0;$i<$ncat+100;$i++)
        {
	$total=0;
        if(!exists $catnames{$i}){next;}
        for($h=$minlen;$h<$maxlen;$h++) # for each length
        	{
		for($j=1;$j<13;$j++)
                        {
                        $total+=$mut_list[$h][$i][$j][0];
                        for($k=1;$k<13;$k++)
                                {
                                $total+=$mut_list[$h][$i][$j][$k];
                                }
                        }
                $total+=$mut_list[$h][$i][0][0];
		}
	$ratio=100*$total/$grand_total;
	if($ratio>3)
		{
		print OUT "\n\n$catnames{$i} ($ratio % of total)\nRead length\td0\td1 T>C\td1 other\td2\n";
	        for($h=3;$h<$maxlen;$h++) # for each length
        	        {
			print OUT "$h\t$mut_list[$h][$i][0][0]\t$mut_list[$h][$i][5][0]\t";
			$total1=0;
			$total2=0;
                	for($j=1;$j<13;$j++)
                        	{
	                        $total1+=$mut_list[$h][$i][$j][0];
        	                for($k=1;$k<13;$k++)
                	                {
                        	        $total2+=$mut_list[$h][$i][$j][$k];
                                	}
                        	}
			$total1-=$mut_list[$h][$i][5][0];
			print OUT "$total1\t$total2\n";
			}
                }
        }



#
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
                $ncat=$n;
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
