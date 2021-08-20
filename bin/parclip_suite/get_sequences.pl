#!/usr/bin/perl

#
# PARCLIP suite, v1.0
# Pavel Morozov, RNA Biology laboratory, HHMI/Rockefeller University, 2016
# 

# Types of mismatches ##################################

# 0 for none
$mism{"AT"}=1; # A>T
$mism{"AC"}=2; # A>C
$mism{"AG"}=3; # A>G ###
$mism{"TA"}=4; # T>A
$mism{"TC"}=5; # T>C ##########
$mism{"TG"}=6; # T>G
$mism{"CA"}=7; # C>A
$mism{"CT"}=8; # C>T
$mism{"CG"}=9; # C>G
$mism{"GA"}=10; # G>A
$mism{"GT"}=11; # G>T
$mism{"GC"}=12; # G>C
$mism{"A>T"}=1; # A>T
$mism{"A>C"}=2; # A>C
$mism{"A>G"}=3; # A>G ###
$mism{"T>A"}=4; # T>A
$mism{"T>C"}=5; # T>C ##########
$mism{"T>G"}=6; # T>G
$mism{"C>A"}=7; # C>A
$mism{"C>T"}=8; # C>T
$mism{"C>G"}=9; # C>G
$mism{"G>A"}=10; # G>A
$mism{"G>T"}=11; # G>T
$mism{"G>C"}=12; # G>C
$mism{"A->T"}=1; # A>T
$mism{"A->C"}=2; # A>C
$mism{"A->G"}=3; # A>G ###
$mism{"T->A"}=4; # T>A
$mism{"T->C"}=5; # T>C ##########
$mism{"T->G"}=6; # T>G
$mism{"C->A"}=7; # C>A
$mism{"C->T"}=8; # C>T
$mism{"C->G"}=9; # C>G
$mism{"G->A"}=10; # G>A
$mism{"G->T"}=11; # G>T
$mism{"G->C"}=12; # G>C
$mism{"A-T"}=1; # A>T
$mism{"A-C"}=2; # A>C
$mism{"A-G"}=3; # A>G ###
$mism{"T-A"}=4; # T>A
$mism{"T-C"}=5; # T>C ##########
$mism{"T-G"}=6; # T>G
$mism{"C-A"}=7; # C>A
$mism{"C-T"}=8; # C>T
$mism{"C-G"}=9; # C>G
$mism{"G-A"}=10; # G>A
$mism{"G-T"}=11; # G>T
$mism{"G-C"}=12; # G>C

$mism{"NONE"}=100; # no mutations 
$mism{"ANY"}=0; # any mutation
$mism{"NotAssigned"}=100;

# 13 for Ins / for bowtie2
# 14 for Del / for bowtie2

# Dealing with arguments ####
$n=@ARGV;

if($n<2)
	{
	print "\nUsage:\n./get_sequences -fasta FASTA_FILE -ctab CTAB_FILE -category CATEGORY_NAME -config CONFIG_FILE -minlen NNN -maxlen NNN -mism FIRST_MUT -mism SECOND_MUT\n\n";
        print "Parameters:\n\t-fasta\toriginal collapsed fasta file with all sequences;\n\t-ctab\tctab file from parclip annotation pipeline for this sample;\n";
	print "\t-category\tcategory to select, might be number of category or name\n\t\trefer to the configuration file used for this run\n\t-config\tconfiguration file with path\n\n";
	print "Optional parameters:\n\t-minlen\tminimal length of the read, natural number\n\t-maxlen\tmaximal length of the read, natural number\n";
	print "\t-mism\tsubstitutions, up to two can be specifyed third one and after will be ignored;\n\t\tperfect matches only if no mismatches specifyed; \n\t\tzero substitution if NONE, anything if ANY;\n\t\tmismatches can be written as AT or A->T or A>T.\n\n";
	exit(0);
	}

$mut1=0;$mut2=0;
$minlen=0;
$maxlen=1000;

$mut1="NotAssigned";
$mut2="NotAssigned";

for($i=0;$i<$n;$i++)
	{
	if(index($ARGV[$i],"-")==0) ## we do have an option name
		{
		if($ARGV[$i] eq "-fasta")
			{
			$fasta=$ARGV[$i+1];
			$i++;
			next;
			}
                if($ARGV[$i] eq "-ctab")
                        {
                        $ctab=$ARGV[$i+1];
                        $i++;
                        next;
                        }
                if($ARGV[$i] eq "-category")
                        {
                        $cat=$ARGV[$i+1];
                        $i++;
                        next;
                        }
                if($ARGV[$i] eq "-minlen")
                        {
                        $minlen=$ARGV[$i+1];
                        $i++;
                        next;
                        }
                if($ARGV[$i] eq "-maxlen")
                        {
                        $maxlen=$ARGV[$i+1];
                        $i++;
                        next;
                        }
                if($ARGV[$i] eq "-mism")
                        {
			if($mut1 eq "NotAssigned")
				{
                        	$mut1=$ARGV[$i+1];
                        	}
			else
				{
	                       if($mut2 eq "NotAssigned")
        	                        {
                	                $mut2=$ARGV[$i+1];
                        	        }       
				}
			$i++;
                        next;
                        }

                if($ARGV[$i] eq "-mut2")
                        {
                        $mut2=$ARGV[$i+1];
                        $i++;
                        next;
                        }
               if($ARGV[$i] eq "-config")
                        {
                        $config=$ARGV[$i+1];
                        $i++;
                        next;
                        }
		print "Can not recognize option \"$ARGV[$i]\"\n";
		exit(0);
		}
	else
		{
        	print "Can not recognize option \"$ARGV[$i]\"\n";
        	exit(0);
		}
	}

# clarifying mismatch
if($mut1=~/\D/)
	{
	if(exists $mism{$mut1})
		{
		$mut1=$mism{$mut1};
		}
	else
		{
		print "Can not identify mismatch $mut1\n";
		exit(0);
		}
	}

if($mut2=~/\D/)
        {
        if(exists $mism{$mut2})
                {$mut2=$mism{$mut2};}
        else
                {print "Can not identify mismatch $mut2\n";exit(0);}
        }

# read the configuration file and initialization

$db=""; 	# anotation database, bowtie formatted
$hier="";	# hierarchy file
$gnm="";	# genome, bowtie formatted
$trans="";	# transcriptome, optional

read_configuration($config);

# it could be a list, divided by commas
@ctgr=split(/\,/,$cat);
$nctgr=@ctgr;
for($f=0;$f<$nctgr;$f++)
	{
	if($ctgr[$f] eq "ANY")
		{
		for($j=0;$j<300;$j++)
			{$cats{$j}=1;}
		last;
		}
	else
		{
		if($ctgr[$f]=~/\D/) # is there nondogotal symbols?
		        {
			# yes, there are, looking for what it is
		        for($i=0;$i<200;$i++)
        		        {
                		if(exists $catnames{$i})
	        	                {
        	        	        if($catnames{$i} eq $ctgr[$f])
                                		{
		                                $cats{$i}=1;
                		                last;
                                		}
		                        }
                		}
        		}
		else
			{
			$cats{$ctgr[$f]}=1;
			}
		}
	}

@array=keys(%cats);
$ncats=@array;

if($ncats==0)
        {
	print "Can not recognize category \"$cat\"\n";
	exit(0);
	}

# Go through the ctab file and memorize right reads
open(IN,$ctab) || die("Can't open CTAB file: $ctab\n");
while(<IN>)
	{
	chomp();
	@s=split(/[\|\t\-\s]/);
	if(exists $cats{$s[2]})
		{
		if($mut1 == 100) # no mutations
			{
			if($mut2==100)
				{
				if($s[3] == 0 && $s[4]==0)
					{
					$chk{$s[0]}=0;
					}
				}
			else
				{
				if($mut2==0)
					{
					if($s[3]>0 && $s[4]==0) # any single
						{
						$chk{$s[0]}=0;
						}
					}
				else
					{
					if($s[3]==$mut2 && $s[4]==0) 
                                                        {
                                                        $chk{$s[0]}=0;
                                                        }
					}
				}
			}
		else
			{
			if($mut1 == 0) # any mutation
				{
				# two any mutations
				if($mut2==0)
					{ 
					if($s[3]>0 && $s[4]>0)
						{
						$chk{$s[0]}=0;
						}
					}
				else
					{
					# one any mutation
					if($mut2==100)
						{
						if($s[3]>0 && $s[4]==0)
							{
							$chk{$s[0]}=0;
							}
						}
					else
						{
						if(($s[3]==$mut2 && $s[4]>0) || ($s[4]==$mut2 && $s[3]>0))
                                                        {
                                                        $chk{$s[0]}=0;
                                                        }
						} 
					}
				}
			else
				{
				if($mut2==100) # first yes, but no second mutation
					{
					# no second mutation and first is equal to requested
					if($s[3] == $mut1 && $s[4]==0)
						{ 
						$chk{$s[0]}=0;
						}						
					}
				else
					{
					if($mut2==0) # any second mutation
						{
						if(($s[3]==$mut1 && $s[4]>0) || ($s[4]==$mut1 && $s[3]>0))
							{
							$chk{$s[0]}=0;
							} 
						}
					else 
						{ # both mutations
						if(($s[3]==$mut1 && $s[4]==$mut2) || ($s[4]==$mut1 && $s[3]==$mut2))
							{
                                                        $chk{$s[0]}=0;
                                                        }
						}
					}
				}
			}
		}
	}
close(IN);

# now walk through the fasta file
open(IN,$fasta) || die("Can't open FASTA file: $fasta\n");
while(<IN>)
	{
	if(/^\>/)
		{
		($nam,$cnt)=split(/[\|\-]/);
		$nam=~s/\>//;
		if(exists $chk{$nam})
			{
			$name=$_;
			}
		else
			{$name="x";}
		}
	else
		{
		if($name ne "x")
			{
			chomp();
			$l=length($_);
			if($l<=$maxlen && $l>=$minlen)
				{
				print "$name$_\n";
				}
			}
		}
	}


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

