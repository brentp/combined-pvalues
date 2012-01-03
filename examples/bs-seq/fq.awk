BEGIN{FS="\t"}
(NR % 4 == 1){ name=$1; seq=""; quals=""; }
(NR % 4 == 2){ seq=$1 }
(NR % 4 == 0){
    quals=$1;
    if(length(seq) == length(quals) && length(seq) < 102){ print name"\n"seq"\n+\n"quals;
    }
}
