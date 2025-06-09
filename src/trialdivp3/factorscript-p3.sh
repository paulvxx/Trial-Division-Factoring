#!/bin/bash

# Define the path to your filelist
filelist="primes_list.txt"

# Check if the file exists
if [ ! -f "$filelist" ]; then
  echo "Error: File '$filelist' not found." >&2
  exit 1
fi

while IFS= read -r prime; do
  echo "Processing $prime"
  if ./trial-const -q 3 $prime -lo 10000000 2147483648 -t 8; then
     ./trial-const -q 3 $prime 2147483648 -t 8 >> output3.txt; 
  fi
done < "$filelist"

echo "Range complete"
exit 0