#!/bin/bash

set -x
set -e

# This script runs the Graph Algorithm for different Partitioning Strategies

GALOIS_ROOT="$(pwd)/../.."
GRAPH_DIR="/proj/prismgt-PG0/vrao79/galois-graphs"
GRAPH_NAMES=("twitter7" "uk-2005")

LOG_DIR="$GALOIS_ROOT/lonestar/analytics/cpu/skywalker/logs"

CC_CMD="$GALOIS_ROOT/build/lonestar/analytics/cpu/connected-components/connected-components-cpu --noverify -t=$(nproc) --symmetricGraph"
SSSP_CMD="$GALOIS_ROOT/build/lonestar/analytics/cpu/sssp/sssp-cpu --noverify -t=$(nproc)"
PR_CMD="$GALOIS_ROOT/build/lonestar/analytics/cpu/pagerank/pagerank-push-cpu --noverify -t=$(nproc)"
# BC_CMD="$GALOIS_ROOT/build/lonestar/analytics/cpu/betweennesscentrality/betweennesscentrality-async-cpu --noverify -t=$(nproc)"
BFS_CMD="$GALOIS_ROOT/build/lonestar/analytics/cpu/bfs/bfs-cpu --noverify -t=$(nproc)"

for GRAPH_NAME in "${GRAPH_NAMES[@]}"
do
    GRAPH_PATH="$GRAPH_DIR/$GRAPH_NAME.gr"

    mkdir -p $LOG_DIR/$GRAPH_NAME
    cd $LOG_DIR/$GRAPH_NAME

    for i in {1..10}
    do
        LOG_FILE="skywalker_CC.log"
        /usr/bin/time -v $CC_CMD "$GRAPH_PATH" &>> $LOG_FILE || { echo "Error: skywalker failed"; cd "$GALOIS_ROOT"; exit; }

        LOG_FILE="skywalker_SSSP.log"
        /usr/bin/time -v $SSSP_CMD "$GRAPH_PATH" &>> $LOG_FILE || { echo "Error: skywalker failed"; cd "$GALOIS_ROOT"; exit; }

        LOG_FILE="skywalker_PR.log"
        /usr/bin/time -v $PR_CMD "$GRAPH_PATH" &>> $LOG_FILE || { echo "Error: skywalker failed"; cd "$GALOIS_ROOT"; exit; }

        # LOG_FILE="skywalker_BC.log"
        # /usr/bin/time -v $BC_CMD "$GRAPH_PATH" 2&>> $LOG_FILE || { echo "Error: skywalker failed"; cd "$GALOIS_ROOT"; exit; }

        LOG_FILE="skywalker_BFS.log"
        /usr/bin/time -v $BFS_CMD "$GRAPH_PATH" &>> $LOG_FILE || { echo "Error: skywalker failed"; cd "$GALOIS_ROOT"; exit; }
    done

    cd $GALOIS_ROOT

done
