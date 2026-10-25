#!/bin/bash
# Set up the annotation course working directory.
# Run with:  source setup_workdir.sh
# (use "source", not "bash", so the cd and variables stay in your current shell)
 
export COURSEDIR=/data/courses/genome_annotation_resources/2026
export WORKDIR=/data/users/ksales/annotation_course
 
mkdir -p "$WORKDIR/scripts" "$WORKDIR/logs"
cd "$WORKDIR" || return 1
 
echo "COURSEDIR = $COURSEDIR"
echo "WORKDIR   = $WORKDIR"
echo "Now in: $(pwd)"
ls