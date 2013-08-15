#!/bin/sh
#need a wrapper skript because BOSS and ComPWA use different versions of libstdc++
#you should source this scrpt in order to avoid opening a new shell
LD_LIBRARY_PATH=/usr/lib64:/home/weidenka/opt/lib/root:/home/weidenka/rootAnalysis/ComPWA/lib:/home/weidenka/work/external/boost_1_48_0/build/lib:/home/weidenka/work/external/qft++/lib:/home/weidenka/work/external/Minuit2-5.28.00/build/lib:$LD_LIBRARY_PATH /home/weidenka/rootAnalysis/ComPWA/bin/ampSvc $@ |grep "rEsUlT"|awk '{print $2}END{}'
