#!/bin/bash

FIRST=$1
LAST=$2

qdel `seq -f "%.0f" $FIRST $LAST`
