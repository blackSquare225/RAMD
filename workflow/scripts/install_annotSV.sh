#!/bin/bash

cd resources/	
git clone https://github.com/lgmgeo/AnnotSV.git

cd AnnotSV/
make PREFIX=. install
make PREFIX=. install-human-annotation
