#!/bin/sh

# Initialize variables
output=""
input_files=""

# Process command-line arguments
while [ $# -gt 0 ]; do
	case "$1" in
	-o)
		output="$2"
		shift 2
		;;
	*)
		input_files="$input_files $1"
		shift
		;;
	esac
done

# Check if output file is specified
if [ "$output" = "" ]; then
	echo "Error: Output filename (-o) is required." >&2
	exit 1
fi

# Check if there are input files
if [ "$input_files" = "" ]; then
	echo "Error: At least one input file is required." >&2
	exit 1
fi

# Run montage to stitch the images
montage "$input_files" -tile x1 -geometry +0+0 "$output"
