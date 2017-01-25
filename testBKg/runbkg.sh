for chan in 4e #2e2mu 4mu 4e #4mu
	do
	for cate in 4 #2 0 1 
		do
		for highmass in 0 #0 1 2
			do
				./job_bkg.lsf "$chan" $cate $highmass
			done
			done
			done

