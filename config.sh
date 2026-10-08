#!/usr/bin/env bash
export DATA=/path/to/data            # root holding reads, Hi-C, and the eval reference
export POLYSPLIT=/path/to/PolySplit  # this repository
export SUBPHASER=/path/to/SubPhaser  # only for the scaffold-first baseline
export FLYE_BIN=/path/to/flye/bin    # dir containing flye + flye-minimap2 + flye-samtools
export KMC_BIN=/path/to/kmc/bin      # dir containing kmc + kmc_tools + kmc_dump (read-direct control)
