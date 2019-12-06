FROM nvidia/cuda

MAINTAINER Dejun Teng<tengdj@gmail.com>
ENV JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64/ \
	USER=root
	
RUN apt-get update && \
	DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      openjdk-8-jdk \
      net-tools \
      curl \
      netcat \
      gnupg \
      apt-utils \
      gcc g++ cmake libboost-all-dev libgeos++-dev \
      freeglut3-dev libcgal-dev software-properties-common \
      wget git vim

# install eigen
RUN wget http://bitbucket.org/eigen/eigen/get/3.3.7.tar.bz2 && \
    tar xvf 3.3.7.tar.bz2 && \
    mv eigen-eigen-323c052e1731/Eigen /usr/local/include && \
    rm -rf eigen-eigen-323c052e1731 && \
    rm 3.3.7.tar.bz2



