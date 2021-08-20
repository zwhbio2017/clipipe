#!/home/CLIPipe_user/bin/miniconda3/envs/clipipe/bin/perl

#
# PARCLIP suite, v1.0
# Pavel Morozov, RNA Biology laboratory, HHMI/Rockefeller University, 2016
# 

open(IN,$ARGV[0]) || die();
while(<IN>)
        {
        chomp();
        if(/^\@/)
                {
                $title=$_;
                $seq=<IN>;
                $info=<IN>;
                $qual=<IN>;
                print "\>$title\n$seq";
                }
        }
close(IN);
