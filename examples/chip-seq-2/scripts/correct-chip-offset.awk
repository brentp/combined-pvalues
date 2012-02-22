BEGIN { 
	FS = OFS = "\t"

	if(!half) {
        print "must specify a value to shift. e.g. -v half=120"
        exit;
    }
}

{
	if ( $4 == "+" ) { 
		# positive strand
        $2 += half
        $3 += half
	} else { 
		# negative strand
        $2 -= half
        if($2 < 0){ next; }
        $3 -= half
	} 
    print $0
}
