#!/bin/bash

# ==========================================
# Configuration
# ==========================================

BUILD_DIR="build"
OUTPUT_DIR="gas_table"

START=38000
END=39000
STEP=100

# ==========================================
# Create directories
# ==========================================

mkdir -p "$BUILD_DIR"
mkdir -p "$OUTPUT_DIR"

# ==========================================
# Build project
# ==========================================

echo "=================================="
echo "Building project..."
echo "=================================="

cd "$BUILD_DIR" || exit 1

cmake ..

make -j$(nproc)

# Check build success
if [ $? -ne 0 ]; then
    echo "Build failed!"
    exit 1
fi

cd ..

# ==========================================
# Run field scan
# ==========================================

CURRENT=$START

while (( CURRENT < END ))
do
    NEXT=$((CURRENT + STEP))

    echo "=================================="
    echo "Generating:"
    echo "$CURRENT -> $NEXT V/cm"
    echo "=================================="

    ./build/generate_gas $CURRENT $NEXT

    # Move generated gas file
    mv TIFRH_${CURRENT}*_*.gas "$OUTPUT_DIR/" 2>/dev/null

    CURRENT=$NEXT
done

echo "=================================="
echo "All jobs completed."
echo "Gas tables saved in:"
echo "$OUTPUT_DIR/"
echo "=================================="
