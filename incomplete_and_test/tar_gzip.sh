#!/bin/bash

DIR=$1
tar -cvzf "${DIR}.tar.gz" $DIR && rm -R $DIR
