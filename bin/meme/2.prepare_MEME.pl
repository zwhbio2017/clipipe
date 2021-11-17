use strict;

my $file_name = $ARGV[0];

my $training_test = $ARGV[1];

my $genome_fasta = $ARGV[2];

`bedtools getfasta -s -name -fi $genome_fasta -bed $file_name.${training_test}_peak.bed -fo temp.fa`;

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
                $_ =~ tr/ACGTacgt/ACGTACGT/;
                print OUTPUT "$_\n";
        }
}

unlink("temp.fa");
