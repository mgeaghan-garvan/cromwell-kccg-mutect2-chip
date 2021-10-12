#!/bin/bash
curl -X POST "http://localhost:CROMWELL_PORT_TO_SED/api/workflows/v1/${1}/abort"