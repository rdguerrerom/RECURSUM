#!/bin/bash
# Initialize vectorclass dependency

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXTERNAL_DIR="${SCRIPT_DIR}/external"
VCL_DIR="${EXTERNAL_DIR}/vectorclass"
INCLUDE_DIR="${SCRIPT_DIR}/include/recursum"

echo "Initializing Vector Class Library..."

# Create external directory if it doesn't exist
mkdir -p "${EXTERNAL_DIR}"

# Clone vectorclass if it doesn't exist
if [ ! -d "${VCL_DIR}" ]; then
    echo "Cloning vectorclass from GitHub..."
    git clone --depth 1 git@github.com:vectorclass/version2.git "${VCL_DIR}"
    echo "✓ Cloned vectorclass"
else
    echo "✓ vectorclass already exists"
fi

# Create include directory
mkdir -p "${INCLUDE_DIR}"

# Create symlinks for vectorclass headers
echo "Creating header symlinks..."
ln -sf "${VCL_DIR}/vectorclass.h" "${INCLUDE_DIR}/vectorclass.h"
ln -sf "${VCL_DIR}/instrset.h" "${INCLUDE_DIR}/instrset.h"
ln -sf "${VCL_DIR}/instrset_detect.cpp" "${INCLUDE_DIR}/instrset_detect.cpp"
ln -sf "${VCL_DIR}/vectormath_exp.h" "${INCLUDE_DIR}/vectormath_exp.h"
ln -sf "${VCL_DIR}/vectormath_trig.h" "${INCLUDE_DIR}/vectormath_trig.h"
ln -sf "${VCL_DIR}/vectormath_hyp.h" "${INCLUDE_DIR}/vectormath_hyp.h"
ln -sf "${VCL_DIR}/vectormath_common.h" "${INCLUDE_DIR}/vectormath_common.h"

# Core vector classes
for file in "${VCL_DIR}"/vector*.h; do
    if [ -f "$file" ]; then
        basename=$(basename "$file")
        if [ ! -L "${INCLUDE_DIR}/${basename}" ]; then
            ln -sf "$file" "${INCLUDE_DIR}/${basename}"
        fi
    fi
done

echo "✓ Vector Class Library initialized successfully!"
echo ""
echo "Headers available in: ${INCLUDE_DIR}/"
