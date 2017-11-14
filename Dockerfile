FROM ubuntu:xenial
ENV ROOTSYS /rootbuild/root
RUN echo "Before installation" && \
	rm /bin/sh && ln -s /bin/bash /bin/sh &&\
    apt-get update -qq && \
    apt-get install -y wget libboost-all-dev libgsl0-dev git cmake g++ pkg-config &&\
    apt-get install -y libxmp4 &&\
    wget https://root.cern.ch/download/root_v6.08.04.Linux-ubuntu16-x86_64-gcc5.4.tar.gz &&\
	mkdir rootbuild && cd rootbuild &&\
    tar xpvfz ../root_*.tar.gz &&\
	rm ../root_*.tar.gz &&\
	cd / &&\
	git clone --depth=10 https://github.com/ComPWA/ComPWA.git &&\
	echo $ROOTSYS &&\
	mkdir ComPWA_build &&\
	cd ComPWA_build &&\
	cmake ../ComPWA &&\
	#cmake -DROOTSYS=/rootbuild/root ../ComPWA &&\
	make -j5 &&\
	echo "ComPWA ready!"
CMD echo "ComPWA start!"
