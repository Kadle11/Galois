#!/bin/bash

set -x
set -e

# This script runs the Graph Algorithm for different Partitioning Strategies

GALOIS_ROOT="$(pwd)/../../"

GRAPH_PATH="$1"
ALGO="$2"
GRAPH_NAME=$(basename "$GRAPH_PATH")
GRAPH_NAME="$(echo $GRAPH_NAME | sed 's/\..*//g')"

SKYWALKER_CMD="$GALOIS_ROOT/build/lonestar/analytics/cpu/skywalker/skywalker-cpu"
BUILD_DIR="$GALOIS_ROOT/build"
LOG_DIR="$GALOIS_ROOT/build/lonestar/analytics/cpu/skywalker/logs/$GRAPH_NAME"

mkdir -p $LOG_DIR

rm -f $BUILD_DIR/CMakeCache.txt

cd "$BUILD_DIR" || { echo "Error: Unable to change directory to $BUILD_DIR"; exit; }
cmake .. -DPARTITIONS=1 -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_BUILD_TYPE=Release || { echo "Error: CMake failed"; cd "$GALOIS_ROOT"; exit; }
make -j $(nproc) skywalker-cpu || { echo "Error: Make failed"; cd "$GALOIS_ROOT"; exit; }

cd $LOG_DIR

LOG_FILE="skywalker_$ALGO"_Baseline.log
$SKYWALKER_CMD "$GRAPH_PATH" -t=$(nproc) > $LOG_FILE || { echo "Error: skywalker failed"; cd "$GALOIS_ROOT"; exit; }


for i in 2 4 8 16 32 64 128 256
do
    rm -f $BUILD_DIR/CMakeCache.txt

    cd "$BUILD_DIR" || { echo "Error: Unable to change directory to $BUILD_DIR"; exit; }
    cmake .. -DPARTITIONS=$i  -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_BUILD_TYPE=Release || { echo "Error: CMake failed"; cd "$GALOIS_ROOT"; exit; }
    make -j $(nproc) skywalker-cpu || { echo "Error: Make failed"; cd "$GALOIS_ROOT"; exit; }

    cd $LOG_DIR

    LOG_FILE="skywalker_$ALGO"_NaivePartitiong_"$i".log
    $SKYWALKER_CMD "$GRAPH_PATH" -t=$(nproc) > $LOG_FILE || { echo "Error: skywalker failed"; cd "$GALOIS_ROOT"; exit; }

done 

for i in 2 4 8 16 32 64 128 256
do
    rm -f $BUILD_DIR/CMakeCache.txt

    cd "$BUILD_DIR" || { echo "Error: Unable to change directory to $BUILD_DIR"; exit; }
    cmake .. -DPARTITIONS=$i -DENABLE_METIS=ON -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_BUILD_TYPE=Release || { echo "Error: CMake failed"; cd "$GALOIS_ROOT"; exit; }
    make -j $(nproc) skywalker-cpu || { echo "Error: Make failed"; cd "$GALOIS_ROOT"; exit; }

    cd $LOG_DIR

    LOG_FILE="skywalker_$ALGO"_MetisPartitiong_"$i".log
    $SKYWALKER_CMD "$GRAPH_PATH" -t=$(nproc) > $LOG_FILE || { echo "Error: skywalker failed"; cd "$GALOIS_ROOT"; exit; }

done 


