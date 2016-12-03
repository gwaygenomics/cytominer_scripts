#!/bin/bash

function CHECK_PATH {

    STATE=$1

    PATHVAR=$2

    if [[ ($STATE != "EXISTS") && ($STATE != "NOT_EXISTS") ]]; then	
	(>&2 echo "Unknown state: $STATE")

	exit 1
    fi

    tstamp=`date`

    if [[ $STATE == "EXISTS" && ! -a "$PATHVAR" ]]; then
	MESSAGE="not created / does not exist."

	(>&2 echo "[$tstamp] ${PATHVAR}" $MESSAGE "Exiting.")
	    
	exit 1
    fi
    
    if [[ $STATE == "NOT_EXISTS" && -a "$PATHVAR" ]]; then
	    MESSAGE="exists."

	    (>&2 echo "[$tstamp] ${PATHVAR} $MESSAGE")

	    while true; do
		read -p "Overwrite? (Y/N)" yn
		case $yn in
		    [Yy]* ) return 1 ;;
		    [Nn]* ) return 2 ;;
		    * ) echo "Please answer yes or no.";;
		esac
	    done

    fi

    return 0
}

function INFO {
    tstamp=`date`

    (>&2 echo "[$tstamp] $1")
}

function COMPARE_MD5() {

    A=$1

    B=$2

    MD5_A=`md5sum ${A} | cut -f1 -d' '`

    MD5_B=`md5sum ${B} | cut -f1 -d' '`
    
    if [[ $MD5_A != $MD5_B ]]; then
	
	tstamp=`date`
	
	(>&2 echo "[$tstamp] md5 check failed: ${A} and ${B} are different")

	exit 1

    fi
}

function CREATE_AND_CHECK_DIR () {

    DNAME=$1

    mkdir -p ${DNAME}

    CHECK_PATH EXISTS ${DNAME}
    
    echo `readlink -e $DNAME`
}

function CREATE_BACKEND_FILE() {

    INFO "Creating ${BACKEND_FILE}"

    CHECK_PATH NOT_EXISTS $BACKEND_FILE

    CHECK_RESULT=$?

    if [[ $CHECK_RESULT = 0 || $CHECK_RESULT = 1 ]]; then
	rm -rf ${BACKEND_FILE}

	time ingest $PLATE_DIR -o sqlite:///${BACKEND_FILE} -c ingest_config.ini --no-munge
    fi

    CHECK_PATH EXISTS $BACKEND_FILE

    INFO "Indexing ${BACKEND_FILE}"

    time sqlite3 ${BACKEND_FILE} < indices.sql

}

function CREATE_AGGREGATED_FILE () {

    INFO "Aggregating ${BACKEND_FILE}"

    CHECK_PATH NOT_EXISTS $AGGREGATED_FILE

    CHECK_RESULT=$?

    if [[ $CHECK_RESULT == 0 || $CHECK_RESULT == 1 ]]; then
	rm -rf $AGGREGATED_FILE

	time Rscript -e "extends <- methods::extends; source('create_profiles.R')" ${BACKEND_FILE} ${AGGREGATED_FILE}

    fi

    CHECK_PATH EXISTS $AGGREGATED_FILE
}

function ARCHIVE_BACKEND_FILE () {

    INFO "Moving ${BACKEND_FILE} to ${BACKEND_ARCHIVE_FILE}"

    rsync -a ${BACKEND_FILE} ${BACKEND_ARCHIVE_FILE}

    COMPARE_MD5 ${BACKEND_FILE} ${BACKEND_ARCHIVE_FILE}

    rm ${BACKEND_FILE}

}

function ARCHIVE_AGGREGATED_FILE () {

    INFO "Moving ${AGGREGATED_FILE} to ${AGGREGATED_ARCHIVE_FILE}"

    rsync -a ${AGGREGATED_FILE} ${AGGREGATED_ARCHIVE_FILE}

    COMPARE_MD5 ${AGGREGATED_FILE} ${AGGREGATED_ARCHIVE_FILE}

    rm ${AGGREGATED_FILE}
}

programname=$0

ingest --help >/dev/null 2>&1 || { echo >&2 "ingest required but not installed. Exiting."; exit 1; }

while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -b|--batchid)
    BATCH_ID="$2"
    shift
    ;;
    -p|--plate)
    PLATE_ID="$2"
    shift
    ;;
    -t|--tmpdir)
    TMP_DIR="$2"
    shift
    ;;
    *)
    echo Unknown option
    ;;
esac
shift
done

TMP_DIR="${TMP_DIR:-/tmp}"

PIPELINE="${PIPELINE:-analysis}"

for var in BATCH_ID PIPELINE PLATE_ID TMP_DIR;
do 
    if [[  -z "${!var}"  ]];
    then
        echo ${var} not defined.
        exit 1
    fi
done

PLATE_DIR=`readlink -e ../../analysis/${BATCH_ID}/${PLATE_ID}/${PIPELINE}/`

BACKEND_DIR=${TMP_DIR}/${BATCH_ID}/${PLATE_ID}/

BACKEND_FILE=${BACKEND_DIR}/${PLATE_ID}.sqlite

AGGREGATED_FILE=${BACKEND_DIR}/${PLATE_ID}.csv

BACKEND_ARCHIVE_DIR="../../backend/${BATCH_ID}/${PLATE_ID}/"

BACKEND_ARCHIVE_FILE=${BACKEND_ARCHIVE_DIR}/${PLATE_ID}.sqlite

AGGREGATED_ARCHIVE_FILE=${BACKEND_ARCHIVE_DIR}/${PLATE_ID}.csv

AGGREGATED_WITH_METADATA_ARCHIVE_FILE=${BACKEND_ARCHIVE_DIR}/${PLATE_ID}_augmented.csv

if [[ (! -e $BACKEND_ARCHIVE_FILE) || (! -e $AGGREGATED_ARCHIVE_FILE) ]]; then

    CHECK_PATH EXISTS ${PLATE_DIR}
    
    BACKEND_DIR=$(CREATE_AND_CHECK_DIR $BACKEND_DIR)

    CREATE_BACKEND_FILE

    CREATE_AGGREGATED_FILE

    BACKEND_ARCHIVE_DIR=$(CREATE_AND_CHECK_DIR $BACKEND_ARCHIVE_DIR)

    ARCHIVE_BACKEND_FILE

    ARCHIVE_AGGREGATED_FILE

else
    
    INFO "$BACKEND_ARCHIVE_FILE and $AGGREGATED_ARCHIVE_FILE already exist. Will not recreate them."

fi

# join with metadata

Rscript -e "extends <- methods::extends; source('join_metadata.R')" $BATCH_ID $PLATE_ID

CHECK_PATH EXISTS $AGGREGATED_WITH_METADATA_ARCHIVE_FILE

