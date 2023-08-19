# Function to count aligned reads in a SAM file
count_aligned_reads() {
    local file="$1"
    local count=$(bunzip2 -c "$file" | grep -vc "^@")
    echo -e "${file%.sam.bz2}\t$count" >> metaphlan_readcounts.tsv
}

rm -rf  metaphlan_readcounts.tsv

# Main script
echo -e "File\tRead Count"  >> metaphlan_readcounts.tsv
echo "----------------"

# Loop through all .sam.bz2 files in the current directory
for file in S*.sam.bz2; do
	echo $file
    if [ -e "$file" ]; then
        count_aligned_reads "$file"
    fi
done