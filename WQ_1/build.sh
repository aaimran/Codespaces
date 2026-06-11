#!/usr/bin/env bash

rm -rf build

# If script was started with /bin/sh, re-exec with bash for bash-specific features
if [ -z "$BASH_VERSION" ]; then
  exec bash "$0" "$@"
fi

mkdir -p build
cd build || exit 1

echo "Detecting Fortran compiler"
if command -v gfortran >/dev/null 2>&1; then
  FC=gfortran
elif command -v ifort >/dev/null 2>&1; then
  FC=ifort
else
  FC=${FC:-f95}
fi

# Configure Fortran flags to avoid line-truncation errors and related strict warnings
FORTRAN_FLAGS="-ffree-line-length-none -Wno-error=line-truncation"
export FFLAGS="$FORTRAN_FLAGS"

echo "Running cmake, mode Release (Fortran: $FC)"
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER="$FC" -DCMAKE_Fortran_FLAGS="$FORTRAN_FLAGS" ../src

# Detect number of CPU cores (cross-platform)
if command -v nproc >/dev/null 2>&1; then
  CORES=$(nproc)
elif command -v sysctl >/dev/null 2>&1; then
  CORES=$(sysctl -n hw.ncpu 2>/dev/null || true)
else
  CORES=1
fi
CORES=${CORES:-1}

echo "Running make with $CORES cores... Logging to build.log"
# Show progress and any error/fatal/undefined messages live, but preserve the real make exit code.
make -j"$CORES" 2>&1 | tee build.log | grep --line-buffered -E '^\[\s*[0-9]{1,3}%\]|[Ee]rror|[Ff]atal|cannot|undefined reference|cannot rename|module file'
# In bash the pipeline exit statuses are in the PIPESTATUS array.
# Use the exit status of the first pipeline element (make). Default to 1 if unavailable.
MAKE_EXIT=${PIPESTATUS[0]:-1}

if [ "$MAKE_EXIT" -eq 0 ]; then
  echo "Build completed successfully!"
  cd ..
  mkdir -p simulation
  mkdir -p bin
  cp ./build/waveqlab3d bin/.
  cp ./build/pre_wql3d bin/.
else
  echo "Build failed. Re-running make to show full output:"
  make -j"$CORES"
  exit 1
fi

