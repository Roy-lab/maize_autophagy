#!/bin/bash

echo "Running here: $HOSTNAME"
echo "arguments $1"
tar -xvzf sharedlib.tgz
export LD_LIBRARY_PATH=sharedlib

TAR="run.${1}.tar.gz"
OUTPUT_DIR="run.${1}"


REGULATORS=$2
CLUSTER_ASSIGNMENTS=$3


#first, if the output is .tar.gz file exists we extract it 
#so we can resume the incomplete run
if [ -f $TAR ]
then
	tar xvzf $TAR
else
	#if it doesn't exist, we create an empty one, 
	#in case that the job get evicted before creating the .tar.gz
	#(if we get evicted and .tar.gz file doesn't exist, we go on hold)
	mkdir $OUTPUT_DIR
	tar cvzf $TAR $OUTPUT_DIR
fi

# The following is suggested by Emile at CHTC so jobs
# don't start over at 0 iter after eviction:
clean_up() {
  echo "Job is exiting; tar result"
  tar -cvzf $TAR $OUTPUT_DIR
  rm -rf $OUTPUT_DIR
  rm -rf shredlib*
  exit $1
}

trap clean_up SIGTERM

#here we check if we already have a prediction
#note that the name of this file check based on -k 
#(I was using 300 by default, here I used 100)
if [ -f $OUTPUT_DIR/fold0/prediction_k100.txt ]
then
	#if we have a prediction, it is incomplete
	#so we resume using "-a"

	echo "./merlin -e 1 -d dataset${1}.txt -l ${REGULATORS} -o run.${1}/ -c ${CLUSTER_ASSIGNMENTS} -v 1 -h0.6 -k100 -p -5 -r 4 -a 1"
	./merlin -e 1 -d dataset${1}.txt -l ${REGULATORS} -o run.${1}/ -c ${CLUSTER_ASSIGNMENTS} -v 1 -h 0.6 -k100 -p -5 -r 4 -a 1 &
        PID=$!
        wait $PID
        EXIT_STATUS=$?

else
#if the file doesn't exist, we start from scratch
	echo "./merlin -e 1 -d dataset${1}.txt -l ${REGULATORS} -o run.${1}/ -c ${CLUSTER_ASSIGNMENTS} -v 1 -h0.6 -k100 -p -5 -r 4"
	./merlin -e 1 -d dataset${1}.txt -l ${REGULATORS} -o run.${1}/ -c ${CLUSTER_ASSIGNMENTS} -v 1 -h 0.6 -k100 -p -5 -r 4 &
        PID=$!
        wait $PID
        EXIT_STATUS=$?

fi

# cleanup
# check the exit code of our merlin execute
if [ $EXIT_STATUS -ne 0]
then
  echo "there was an error: $EXIT_STATUS"
fi


clean_up $EXIT_STATUS

