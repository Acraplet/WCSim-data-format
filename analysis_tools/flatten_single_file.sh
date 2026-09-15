#!/usr/bin/env bash
# flatten_one.sh - run flatten_wcsim.C on a single WCSim ROOT file.
#
# Usage:
#   ./flatten_one.sh <input.root> [output.root] [isMDT]
#
#   input.root  - WCSim (or MDT-processed) file to flatten
#   output.root - optional, defaults to "<input>_flat.root" next to the input.
#                 Refuses to run if this ends up equal to the input path:
#                 flatten_wcsim.C RECREATEs the output file while still
#                 reading from the input TFile, so overwriting the input in
#                 place is not safe - pass a genuinely different path.
#   isMDT       - optional passthrough to flatten_wcsim.C (-1 auto-detect
#                 [default], 0 force plain-WCSim truth-matching, 1 force MDT)
#
# Requires WCSIM_BUILD_DIR to be set to the WCSim install directory that has
# lib/libWCSimRoot.so (matching the ROOT version of the `root` on PATH),
# e.g.:  export WCSIM_BUILD_DIR=/eos/user/a/acraplet/WCSim_changes/build/install

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MACRO="${SCRIPT_DIR}/flatten_wcsim.C"

if [ $# -lt 1 ] || [ $# -gt 3 ]; then
    echo "Usage: $0 <input.root> [output.root] [isMDT]" >&2
    exit 1
fi

: "${WCSIM_BUILD_DIR:?Set WCSIM_BUILD_DIR to the WCSim install dir with lib/libWCSimRoot.so before running this script}"

INPUT=$1
OUTPUT=${2:-"${INPUT%.root}_flat.root"}
ISMDT=${3:--1}

if [ ! -f "$INPUT" ]; then
    echo "ERROR: input file not found: $INPUT" >&2
    exit 1
fi

if [ "$INPUT" -ef "$OUTPUT" ] 2>/dev/null || [ "$INPUT" = "$OUTPUT" ]; then
    echo "ERROR: input and output resolve to the same file ($INPUT)." >&2
    echo "flatten_wcsim.C recreates the output while reading the input, so overwriting in place is not safe." >&2
    echo "Pass a different output path as the 2nd argument." >&2
    exit 1
fi

echo "Flattening:"
echo "  input : $INPUT"
echo "  output: $OUTPUT"
echo "  isMDT : $ISMDT"

root -l -b -q "${MACRO}(\"${INPUT}\",\"${OUTPUT}\",${ISMDT})"
