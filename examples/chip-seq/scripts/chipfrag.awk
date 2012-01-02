#
# Reads a BED file and computes the number of overlapping sequence 5' ends
# across strands, user may specify a shift for the positions on the + strand
# 
# For a chip-seq experiment this number will be maximized when the shift is 
# equal with the sequenced DNA fragment size. Note that this is just an estimate
# but it is a good start.
#
# usage:
#
# $ awk -f chipfrag.awk data.bed 
#
# or with paramters:
#
# $ awk -v shift=100 -v max=100000 -f chipfrag.awk data.bed
#
# Outputs the shift and number of matching reads
#
#
# The script should run on a BED file that contains mostly locations
# where binding is expected to occur 
#
# A possible run that would test a range of shifts of 20, 25, 30 to 200
# would look like:
#
# $ for i in $(seq 20 5 200); do awk -v shift=$i -f chipfrag.awk long.bed; done
#
# The line with the largest count indicates the most likely fragment size
#
BEGIN { 
	FS = OFS = "\t"

    # fudge window larger numbers give a smoother plot
    if(!window || window < 0) window = 2

	# get maxsize as a power of 10
	if(!pow) pow = 5 
	
	# maximum number of records to process
	# keeps everything in memory!
	if(!max) max = 10^pow;

	# the shift that is to be applied
	if (!shift) shift = 0
}

{
	if ( $6 == "+" ) { 
		# positive strand
		fwd[$2 + shift]++
	} else { 
		# negative strand
		rev[$3]++
		
		# give it more opportunities to match
		# when in close by
        for(w=-window; w <= window; w+=1) {
            rev[$3 + w]++
        }
	} 

	if (NR > max) {
		exit
	}
}

END {
	# sum up the matching positions
	sum = 0
	for (key in fwd) {
		sum = sum + rev[key]
	}
	# prints the shift and the number of matches
	print shift, sum
}
