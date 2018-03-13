#!/bin/csh

/bin/sh `sed -n ${SGE_TASK_ID}p < $1`

