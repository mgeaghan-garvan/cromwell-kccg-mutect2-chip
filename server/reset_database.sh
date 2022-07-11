#!/bin/bash
echo 'DROP DATABASE DBNAME_TO_SED;
CREATE DATABASE DBNAME_TO_SED;
GRANT ALL PRIVILEGES ON DBNAME_TO_SED.* TO '\''cromwell_user'\''@'\''localhost'\'' WITH GRANT OPTION;
GRANT ALL PRIVILEGES ON DBNAME_TO_SED.* TO '\''cromwell_user'\''@'\''%'\'' WITH GRANT OPTION;' \
| MYSQL_TO_SED/bin/mysql -u root --password=password --socket=MYSQL_RUNDIR_TO_SED/socket