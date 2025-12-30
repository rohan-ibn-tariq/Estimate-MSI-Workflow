#!/bin/bash
set -e

# ============================================
# Setup Pixi Environment for Varlociraptor 
# (local builds)
# Copy this script and source it to set up
# the environment variables for pixi in 
# your local varlociraptor build.
# ============================================

# Set PIXI_ENV to the path of your pixi environment
export PIXI_ENV=""
export CC="$PIXI_ENV/bin/cc"
export CXX="$PIXI_ENV/bin/g++"
export AR="$PIXI_ENV/bin/ar"
export PKG_CONFIG="$PIXI_ENV/bin/pkg-config"
export PKG_CONFIG_PATH="$PIXI_ENV/lib/pkgconfig"
export CFLAGS="-I$PIXI_ENV/include"
export CXXFLAGS="-I$PIXI_ENV/include"
export LDFLAGS="-L$PIXI_ENV/lib"
export BZIP2_DYNAMIC=1
export LD_LIBRARY_PATH="$PIXI_ENV/lib:$LD_LIBRARY_PATH"
