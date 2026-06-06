#!/usr/bin/env bash
# Count lines of code in a given directory, grouped by file extension.
# Usage: ./tools/count_lines.sh [directory]

set -euo pipefail

ROOT="${1:-"$(cd "$(dirname "$0")/.." && pwd)"}"

declare -A ext_lines ext_files
total_lines=0
total_files=0

cd "$ROOT"

tempfile=$(mktemp)
trap 'rm -f "$tempfile"' EXIT

find . \( -name build -o -name .git -o -name .tree -o -name __pycache__ -o -name .claude -o -name vendor \) -prune -o \
    -type f \( \
        -name '*.cppm' -o -name '*.cpp' -o -name '*.cxx' -o -name '*.cc' -o -name '*.c' \
        -o -name '*.hpp' -o -name '*.hxx' -o -name '*.hh' -o -name '*.h' \
        -o -name '*.py' -o -name '*.sh' \
        -o -name '*.cmake' -o -name '*.txt' \
        -o -name '*.md' \
        -o -name '*.json' -o -name '*.yaml' -o -name '*.yml' -o -name '*.toml' \
        -o -name '*.xml' \
    \) -print0 > "$tempfile"

[[ ! -s "$tempfile" ]] && { echo "No source files found."; exit 0; }

while IFS= read -r -d '' file; do
    file="${file#./}"
    ext="${file##*.}"
    [[ "$ext" == "$file" ]] && ext="(none)"

    # handle multi-dot extensions
    case "$file" in
        *.cppm) ext="cppm" ;;
    esac

    n=$(wc -l < "$file")
    ((n == 0)) && continue

    printf "%7d  %s\n" "$n" "$file"

    ext_lines[$ext]=$(( ${ext_lines[$ext]:-0} + n ))
    ext_files[$ext]=$(( ${ext_files[$ext]:-0} + 1 ))
    total_lines=$(( total_lines + n ))
    total_files=$(( total_files + 1 ))
done < "$tempfile"

echo
printf "%-10s %5s %8s\n" "Extension" "Files" "Lines"
printf -- "-------------------------\n"
for ext in "${!ext_lines[@]}"; do
    printf "%-10s %5d %8d\n" ".$ext" "${ext_files[$ext]}" "${ext_lines[$ext]}"
done | sort -k3 -rn
printf -- "-------------------------\n"
printf "%-10s %5d %8d\n" "TOTAL" "$total_files" "$total_lines"
