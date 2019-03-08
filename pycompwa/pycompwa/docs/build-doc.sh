#! /bin/sh

sphinx-apidoc -f -o source/ ../
make html

exit 0