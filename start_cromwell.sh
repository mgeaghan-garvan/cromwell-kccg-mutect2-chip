#!/bin/bash

taskset -c 40-60 java -Xmx16G -Xms16G -Dconfig.file=./workflow/mutect2.conf -jar CROMWELL_JAR_TO_SED server
