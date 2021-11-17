use strict;

my $total_line = `wc -l $ARGV[0]`;

chomp($total_line);

$ARGV[0] =~ /(.*)\.all_peak\.bed/;

my $file_name = $1;

if($total_line>=1000){
        `sort -k 11rn $ARGV[0] | head -n 500 > $file_name.training_peak.bed`;
        `sort -k 11rn $ARGV[0] | head -n 1000 | tail -n 500 > $file_name.test_peak.bed`;
}
else{
        my $cutoff = int($total_line/2);
        `sort -k 11rn $ARGV[0] | head -n $cutoff > $file_name.training_peak.bed`;
        `sort -k 11rn $ARGV[0] | tail -n $cutoff > $file_name.test_peak.bed`;
}