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
  if ./trial-const -q 29 $prime -lo 2097152 1073741824 -t 4; then
     ./trial-const -q 29 $prime 1073741824 -t 2 >> output29.txt; 
  fi
done < "$filelist"

echo "Range complete"
exit 0
