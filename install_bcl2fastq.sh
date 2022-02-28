#!bin/bash

#install bcl2fastq v2.18.0.12
# sudo this

VERSION=bcl2fastq2-v2.18.0.12
export TMP=/slipstream/home/jonathan/bin/tmp
export SOURCE=${TMP}/bcl2fastq2-v2.18.0.12
export BUILD=${TMP}/bcl2fastq2-v2.18.0.12-build
export INSTALL_DIR=/slipstream/home/jonathan/bin/bcl2fastq2-v2.18.0.12

cd ${TMP}
tar -xvzf /slipstream/home/jonathan/bin/bcl2fastq2-v2.18.0.12.tar.gz #set $VERSION later
mv bcl2fastq bcl2fastq2-v2.18.0.12
chmod -R a+x bcl2fastq*

mkdir ${BUILD}
cd ${BUILD}
${SOURCE}/src/configure --prefix=${INSTALL_DIR}

make

make install

#run ~/bin/bcl2fastq2-v2.18.0.12/bin/bcl2fastq

#bootstrap


#install older version

cd /slipstream/home/jonathan/bin/
wget ftp://webdata:webdata@ussd-ftp.illumina.com/Downloads/Software/bcl2fastq/bcl2fastq-1.8.4.tar.bz2
tar -xjvf bcl2fastq-1.8.4.tar.bz2

export TMP=/slipstream/home/jonathan/bin/tmp
export SOURCE=${TMP}/bcl2fastq-v1.8.4
export BUILD=${TMP}/bcl2fastq-v1.8.4-build
export INSTALL=/usr/local/bcl2fastq-v1.8.4
export BOOST_ROOT=/usr/include/boost
mv bcl2fastq ${TMP}/bcl2fastq-v1.8.4
cd ${TMP}

chmod -R a+x bcl2fastq*

mkdir ${BUILD}
cd ${BUILD}
${SOURCE}/src/configure --prefix=${INSTALL_DIR}





