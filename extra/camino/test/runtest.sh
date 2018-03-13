#!/bin/bash

UNAME=`uname`
CYGWIN=`expr ${UNAME} : "CYGWIN"`

SCRIPTDIR=${0%/*}

if [ $CYGWIN == '0' ]; then
    export CLASSPATH=.:${SCRIPTDIR}/test:${SCRIPTDIR}/junit.jar:${SCRIPTDIR}/..
else
    CAMINO=`cygpath -w ${SCRIPTDIR}/..`
    TEST=`cygpath -w ${SCRIPTDIR}/test`
    JUNIT=`cygpath -w ${SCRIPTDIR}/junit.jar`
    export CLASSPATH=".;${CAMINO};${TEST};${JUNIT}"
fi


java -Djava.util.logging.config.file=./logging.properties -Xmx250M -classpath $CLASSPATH AllTests $*

