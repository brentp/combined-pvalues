BEGIN { 
	FS = OFS = "\t"

	if(!half) {
        print "must specify a value to shift. e.g. -v half=120"
        exit;
    }
}

{
    $4 = "."
	if ( $6 == "+" ) { 
		# positive strand
        $2 += half
        $3 += half
		fwd[$2 + shift]++
	} else { 
		# negative strand
        $2 -= half
        $3 -= half
        if($2 < 0){ next; }
	} 
    print $0
}
