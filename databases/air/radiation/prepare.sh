#!/bin/bash
file1=$1"_nist"
file2=$1
sed -e 's/[][()|]//g' -e 's/ - / /g' -e 's/*/star/g'  -e 's/1\/2/0.5/g' -e 's/3\/2/1.5/g' -e 's/5\/2/2.5/g' -e 's/7\/2/3.5/g' -e 's/9\/2/4.5/g' -e 's/11\/2/5.5/g' \
    -e 's/13\/2/6.5/g' -e 's/15\/2/7.5/g' -e 's/17\/2/8.5/g' -e 's/19\/2/9.5/g' -e 's/21\/2/10.5/g' -e 's/23\/2/11.5/g'            \
    -e 's/25\/2/12.5/g' -e 's/27\/2/13.5/g' -e 's/29\/2/14.5/g' -e 's/31\/2/15.5/g' -e 's/33\/2/16.5/g' -e 's/35\/2/17.5/g'        \
    -e 's/37\/2/18.5/g' -e 's/39\/2/19.5/g'  $file1 | awk 'NF' >$file2
