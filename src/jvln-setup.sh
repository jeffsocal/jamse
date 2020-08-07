#!/bin/bash

chmod +x ~/git/jvln-psm/exe/database/*.php
chmod +x ~/git/jvln-psm/exe/server/*.php
# link the php scripts
cd /usr/local/bin/
ln -s ~/git/jvln-psm/exe/database/* ./
ln -s ~/git/jvln-psm/exe/server/* ./
for f in *.php; do mv $f `basename $f .php`; done;
# link the Rscripts
cd /var/jvlnse/bin/
ln -s ~/git/jvln-psm/src/Rscript/* ./
