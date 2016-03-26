#!/bin/sh

SRC=input.cpp
LLVM_BUILD=/h2/roshan/SourceCode/llvm/build
GGC_ROOT=/h2/roshan/SourceCode/ggc
GALOIS_ROOT=/h2/roshan/SourceCode/GaloisCpp

CXX_DEFINES="-DGALOIS_COPYRIGHT_YEAR=2015 -DGALOIS_USE_EXP -DGALOIS_VERSION=2.3.0 -DGALOIS_VERSION_MAJOR=2 -DGALOIS_VERSION_MINOR=3 -DGALOIS_VERSION_PATCH=0 -D__STDC_LIMIT_MACROS"
CXX_FLAGS="-g -Wall -fcolor-diagnostics -Wall -g -I$GALOIS_ROOT/exp/include -I/opt/apps/ossw/libraries/mpich2/mpich2-1.5/c7/clang-system/include -I/net/faraday/workspace/local/modules/tbb-4.2/include -I/opt/apps/ossw/libraries/boost/boost-1.58.0/c7/clang-system/include -I$GALOIS_ROOT/lonestar/include -I$GALOIS_ROOT/libruntime/include -I$GALOIS_ROOT/libnet/include -I$GALOIS_ROOT/libsubstrate/include -I$GALOIS_ROOT/libllvm/include -I$GALOIS_ROOT/build/debug/libllvm/include -I$GALOIS_ROOT/libgraphs/include -std=gnu++11"

CXX=$LLVM_BUILD/bin/clang++
IRGL_CXX="$CXX -Xclang -load -Xclang $LLVM_BUILD/lib/GaloisFunctions.so -Xclang -plugin -Xclang irgl"
GANALYSIS_CXX="$CXX -Xclang -load -Xclang $LLVM_BUILD/lib/GaloisFunctionsAnalysis.so -Xclang -plugin -Xclang galois-analysis"
GFUNCS_CXX="$CXX -Xclang -load -Xclang $LLVM_BUILD/lib/GaloisFunctions.so -Xclang -plugin -Xclang galois-fns"

echo "Cleaning generated files"
rm -f gen.cpp gen_cuda.py gen_cuda.cu gen_cuda.cuh gen_cuda.h
cp $SRC gen.cpp

echo "Generating IrGL code"
$IRGL_CXX $CXX_DEFINES $CXX_FLAGS -o .temp.o -c gen.cpp &>/dev/null
echo "Generating CUDA code from IrGL"
$GGC_ROOT/src/ggc -o gen_cuda.cu gen_cuda.py &>/dev/null

echo "Generating analysis information"
$GANALYSIS_CXX $CXX_DEFINES $CXX_FLAGS -o .temp.o -c gen.cpp &>/dev/null
echo "Generating communication code"
$GFUNCS_CXX $CXX_DEFINES $CXX_FLAGS -o .temp.o -c gen.cpp &>/dev/null

rm -f Entry-*.dot

