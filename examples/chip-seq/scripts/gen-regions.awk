BEGIN { 
    OFS="\t";

    if (!step) step = 100
}
($2 ~ /^[0-9]+$/) {

    for(i=0; i<= $2;i+=step){
        print $1,i,i+step
    }
}
