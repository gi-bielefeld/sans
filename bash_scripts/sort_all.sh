for folder in thread_results; do # original and fingerprint results are already sorted
	cd $folder	
	for file in ./*; do
		 sort -n -k 1 -r $file -o $file
	done
	cd ..
done
