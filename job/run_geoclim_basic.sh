#!/bin/bash

# Run geoclim model
./geoclim.exe 0 1 3 0 0

# Recursive resubmission
test $? -eq 0 && ./submit_chain.sh

