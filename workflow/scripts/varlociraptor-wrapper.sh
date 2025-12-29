#!/bin/bash
set -e

# ============================================
# Varlociraptor Wrapper for Feature Branch
# ============================================

VARLOCIRAPTOR_BIN="$1"
VARLOCIRAPTOR_DIR="$2"
SETUP_SCRIPT="$3"
shift 3

if [ -z "$VARLOCIRAPTOR_BIN" ] || [ -z "$VARLOCIRAPTOR_DIR" ] || [ -z "$SETUP_SCRIPT" ]; then
    echo "ERROR: Missing required arguments" >&2
    echo "Usage: $0 <varlociraptor_binary> <varlociraptor_dir> <setup_script> [varlociraptor_args...]" >&2
    exit 1
fi

if [ ! -f "$VARLOCIRAPTOR_BIN" ]; then
    echo "ERROR: Varlociraptor binary not found: $VARLOCIRAPTOR_BIN" >&2
    exit 1
fi

if [ ! -d "$VARLOCIRAPTOR_DIR" ]; then
    echo "ERROR: Varlociraptor directory not found: $VARLOCIRAPTOR_DIR" >&2
    exit 1
fi

if [ ! -f "$VARLOCIRAPTOR_DIR/pixi.toml" ]; then
    echo "ERROR: pixi.toml not found in: $VARLOCIRAPTOR_DIR" >&2
    exit 1
fi

if [ ! -f "$SETUP_SCRIPT" ]; then
    echo "ERROR: Setup script not found: $SETUP_SCRIPT" >&2
    exit 1
fi

echo "Activating pixi environment..." >&2
cd "$VARLOCIRAPTOR_DIR"

pixi run bash -c "
    source $SETUP_SCRIPT
    $VARLOCIRAPTOR_BIN $*
"
