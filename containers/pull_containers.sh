#!/usr/bin/env bash
# pull_from_nextflow.sh
# Run this inside your ./containers directory

set -euo pipefail

# --- check you are inside containers dir ---
if [[ ! -d . || "$(basename "$PWD")" != "containers" ]]; then
    echo "⚠️  Please cd into the 'containers/' directory before running this script."
    exit 1
fi

# --- detect project root (one level up) ---
PROJECT_DIR="$(cd .. && pwd)"

# --- detect apptainer or singularity ---
if command -v apptainer &>/dev/null; then
    PULLER="apptainer"
elif command -v singularity &>/dev/null; then
    PULLER="singularity"
else
    echo "❌ Neither apptainer nor singularity found in PATH."
    exit 1
fi

echo "🔍 Using $PULLER to pull containers..."
echo "📂 Project directory: $PROJECT_DIR"
echo

# --- extract container URIs from Nextflow config ---
CONTAINERS=$(nextflow config "$PROJECT_DIR" | grep "container =" | grep -oE "['\"]quay.io/[^'\"]+" | tr -d "'" | sort -u)

if [[ -z "$CONTAINERS" ]]; then
    echo "⚠️  No container URIs found in Nextflow config output."
    exit 0
fi

# --- pull each container ---
for IMG in $CONTAINERS; do
    NAME=$(echo "$IMG" | sed 's|/|-|g; s|:|-|g').img

    if [[ -f "$NAME" ]]; then
        echo "✅ Already exists: $NAME"
        continue
    fi

    echo "⬇️  Pulling: docker://$IMG"
    $PULLER pull "$NAME" "docker://$IMG" || echo "❌ Failed: $IMG"
done

echo
echo "✨ All done! Containers are saved in: $PWD"

