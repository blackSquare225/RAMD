#!/bin/bash

wget -c https://github.com/Oshlack/JAFFA/releases/download/version-2.3/JAFFA-version-2.3.tar.gz -O resources/JAFFA-version-2.3.tar.gz

cd resources/

tar -zxvf JAFFA-version-2.3.tar.gz

cd JAFFA-version-2.3

wget -v https://figshare.com/ndownloader/files/25410494 -O JAFFA_REFERENCE_FILES_HG38_GENCODE22.V2.tar.gz

tar -zxvf JAFFA_REFERENCE_FILES_HG38_GENCODE22.V2.tar.gz

./install_linux64.sh

cd ../../

mkdir -p repos/ 

mv resources/JAFFA-version-2.3.tar.gz resources/JAFFA-version-2.3/JAFFA_REFERENCE_FILES_HG38_GENCODE22.V2.tar.gz repos/
