#!/bin/bash
set -e

FILE="../ParallelPlate.C"
DIR="build"

cd "$DIR"

for val in $(seq 1 100); do
	sed -i "" "$FILE"

	rc

	./ParallelPlate

done


