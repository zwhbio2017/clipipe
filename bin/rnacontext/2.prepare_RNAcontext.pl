use strict;

my $file_name = $ARGV[0];

my $training_test = $ARGV[1];

my $genome_fasta = $ARGV[2];

open INPUT1,"${file_name}.${training_test}_peak.bed" or die;

open TEMP,">rename.txt" or die;

my $i;

while(<INPUT1>){
        chomp;
        $i++;
        my @line = split /\t/,$_;
        if($line[2]-$line[1]<60){
                $line[1] = $line[1] - 20;
                if($line[1] < 0){
                        $line[1] = 0;
                }
                $line[2] = $line[2] + 20;
        }
        print TEMP "$line[0]\t$line[1]\t$line[2]\t$i\t$line[4]\t$line[5]\n";
}

`bedtools getfasta -s -name -fi $genome_fasta -bed rename.txt -fo temp.fa`;

open INPUT3,"temp.fa" or die;

open OUTPUT,">$file_name.${training_test}_peak.fa" or die;

my $j;

my $seq_name;

while(<INPUT3>){
        chomp;
        $j++;
        if(($j-1)%2==0){
                $_ =~ />(.*)/;
                $seq_name = $1;
                print OUTPUT ">$seq_name\n";
        }
        elsif(($j-2)%2==0){
                $_ =~ tr/ACGTacgt/ACGUACGU/;
                print OUTPUT "$_\n";
        }
}

unlink("rename.txt");
unlink("temp.fa");
