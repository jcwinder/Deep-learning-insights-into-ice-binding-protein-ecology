cat /path/to/folder/prok_pfam_filepaths.txt | while read line #collected pfam filepaths
do
	acc=$(echo "$line" | awk -F'/' '{print $9}')
	awk -v str="$acc" '{print $0 "\t" str}' "$line" >> "/path/to/folder/prok_pfam_comb.pfam"
done