# cromwell-kccg-mutect2

Run the KCCG GATK4 Mutect2 somatic variant calling pipeline using the Cromwell workflow engine.

## Usage:

### Configuring MySQL server
If there is no MySQL server currently running on the cluster, run the following commands:
```bash
# Choose an empty port number.
MYSQL=/home/glsai/mysql/mysql-5.7.27-linux-glibc2.12-x86_64/
MYSQL_RUNDIR=/home/glsai/mysql/mysql-5.7.27-linux-glibc2.12-x86_64/
sqlport=40008
 
screen -S MySQLPort${sqlport}
 
# Unfortunately have to run these lines again inside screen session
MYSQL=/home/glsai/mysql/mysql-5.7.27-linux-glibc2.12-x86_64/
MYSQL_RUNDIR=/home/glsai/mysql/mysql-5.7.27-linux-glibc2.12-x86_64/
sqlport=40008
 
${MYSQL}/bin/mysqld --no-defaults --user=glsai --datadir=$MYSQL_RUNDIR/data/ --basedir=$MYSQL_RUNDIR/ --log-error=$MYSQL_RUNDIR/log/mysql.err --pid-file=$MYSQL_RUNDIR/mysql.pid --socket=$MYSQL_RUNDIR/socket --port=${sqlport}
 
# Then exit screen session once our SQL is running. When its running it says or does nothing.
# Ctrl + A and D
 
${MYSQL}/bin/mysql -u root -p --socket=$MYSQL_RUNDIR/socket
# password: password
 
# if cromwell_user not already created else skip these two lines
# you can check if it already exists with the following query: "SELECT user FROM mysql.user;"
CREATE USER 'cromwell_user'@'localhost' IDENTIFIED BY '123456781!aA';
CREATE USER 'cromwell_user'@'%' IDENTIFIED BY '12345678!1aA';
```

### Configuring Cromwell server

```bash
# Clone the repository into your working directory
git pull https://git.gimr.garvan.org.au/micgea/cromwell-kccg-mutect2.git
cd cromwell-kccg-mutect2

# Configure environment variables
MYSQL_HOST="0.0.0.0"  # Need to use '0.0.0.0' rather than 'localhost'
MYSQL_PORT=40008
${PROJECT_NAME}="project_name"
CROMWELL_PORT=8000

# Run configuration script
./configure_cromwell.sh ${MYSQL_HOST} ${MYSQL_PORT} ${PROJECT_NAME} ${CROMWELL_PORT}
```

### Running Cromwell server

```bash
# Create a screen session
screen -S "CromwellPort${CROMWELL_PORT}"

# Ensure that Java 11 is the active Java version
java -version
# If it isn't check the currently active modules
module list
# If an older version of Java is loaded, unload it. This should return the default Java version to 8/v1.8
module unload centos6.10/ccg/java/1.7.0_25
# Now load Java 11
module load centos6.10/shacar/java/jdk-11.0.2

# Start the Cromwell server
./start_cromwell.sh
```

Once running, you can detatch the session with Ctrl + A, then D.

### Configuring workflow

Edit ./workflow/inputs.json in your favourite text editor. The input file is a series of key: value pairs. The template provided has values of "REQUIRED" and "OPTIONAL" for required and optional fields, respectively. If not using an optional parameter, simply delete the entire line. Some defaults are also provided and can be left as-is or changed if desired.

### Running the workflow

Once everything is set up and configured, run the workflow as follows:

```bash
# Create a screen session for the run
screen -S ${PROJECT_NAME}

# Submit the workflow to Cromwell
./run.sh
```

Once running, you can detatch the session with Ctrl + A, then D.

You can check on the state of the run by re-attaching the Cromwell screen session:

```bash
screen -r "CromwellPort${CROMWELL_PORT}"
```

Once finished, you can run gather_outputs.sh to move all the relevant output files into workflow_out

```bash
./gather_outputs.sh
```
