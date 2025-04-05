#!/bin/bash

# make EnsembleDesign
make

# Go to ./tools/LinearPartition and run make
cd ./tools/LinearPartition
make
cd ../..

# Go to ./tools/LinearDesign and run make
cd ./tools/LinearDesign
make
cd ../..