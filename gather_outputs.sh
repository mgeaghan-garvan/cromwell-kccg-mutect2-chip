#!/bin/bash
cd workflow_out
find . -type f -exec mv {} . \;
rm -r Mutect2
