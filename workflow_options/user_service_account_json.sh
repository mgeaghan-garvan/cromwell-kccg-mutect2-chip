#!/bin/bash

echo "The service account JSON provided is as follows:"
cat ${1}
echo ""

echo "To supply a user service account when running a workflow through Cromwell, supply the following string to the 'user_service_account_json' key in the workflow options JSON file:"
echo ""

printf '"'
cat ${1} | awk '{gsub("\"", "\\\"", $0); printf "%s", $0 "\\n"}'
printf '"'
echo ""