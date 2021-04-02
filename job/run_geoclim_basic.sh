#!/bin/bash

# Run geoclim model
./geoclim.exe 0 1 3 0 0

# Recursive resubmission
./submit_chain.sh CONTINUE_RUN

