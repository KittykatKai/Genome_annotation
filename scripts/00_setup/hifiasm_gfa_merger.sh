#!/bin/bash
# Convert the hifiasm primary contig graph (*.bp.p_ctg.gfa) to FASTA for EDTA.
# Run with:  bash hifiasm_gfa_merger.sh

set -euo pipefail
 
HIFIASM_DIR=/data/users/ksales/assembly_annotation_course/output_dir/hifiasm
 
# Find the primary contig GFA (ignores hap1/hap2 and the unitig graphs)
GFA=$(ls "$HIFIASM_DIR"/*.bp.p_ctg.gfa 2>/dev/null | head -n 1 || true)
 
if [[ -z "$GFA" ]]; then
    echo "ERROR: no *.bp.p_ctg.gfa found in $HIFIASM_DIR. Files there:" >&2
    ls "$HIFIASM_DIR" >&2
    exit 1
fi
 
FASTA="${GFA%.gfa}.fa"
 
# S lines in a GFA hold the contigs: column 2 = name, column 3 = sequence
awk 'BEGIN{OFS="\n"} /^S/{print ">"$2, $3}' "$GFA" > "$FASTA"
 
# Quick sanity check
awk '/^>/{n++; next}{len+=length($0)} END{printf "Contigs: %d\nTotal length: %.2f Mb\n", n, len/1e6}' "$FASTA"
 
echo "Input:  $GFA"
echo "Output: $FASTA"