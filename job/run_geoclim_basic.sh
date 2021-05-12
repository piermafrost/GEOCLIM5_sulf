#!/bin/bash

# Run geoclim model
./geoclim.exe

# Recursive resubmission
test $? -eq 0 && ./submit_chain.sh

