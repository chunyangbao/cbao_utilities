#!/bin/awk -f
# From: https://unix.stackexchange.com/questions/25138/how-to-print-certain-columns-by-name
# Usage: awk -f t.awk -v cols=name,age,id,name,id input

BEGIN {
    FS="\t"; OFS="\t"
    n=split(cols,out,",")
}
NR==1 {
    for (i=1; i++<NF;)
        ix[$i] = i
}
NR>1 {
    printf $ix[out[1]]
    for (i=1; i++<n;)
        printf "%s%s", "\t", $ix[out[i]]
    print ""
}
