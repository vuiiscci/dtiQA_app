FROM nvidia/cuda:8.0-runtime-ubuntu16.04

RUN apt-get -y update && \ 
    apt-get -y upgrade && \
    apt-get -y install wget && \
    apt-get -y install zip unzip && \
    apt-get -y install gzip && \
    apt-get -y install imagemagick && \
    apt-get -y install xvfb && \
    apt-get -y install dc && \
    apt-get -y install libxrandr2 && \
    apt-get -y install eog && \
    apt-get -y install evince && \
    apt-get -y install default-jdk && \
    apt-get -y install python3 && \
    apt-get -y install python3-pip && \
    apt-get -y install vim && \
    apt-get -y install tar && \
    pip3 install pyyaml && \
    mkdir /tmp/gs && \
    cd /tmp/gs && \
    wget https://github.com/ArtifexSoftware/ghostpdl-downloads/releases/download/gs922/ghostscript-9.22-linux-x86_64.tgz --no-check-certificate && \
    tar xvf ghostscript-9.22-linux-x86_64.tgz && \
    mv ghostscript-9.22-linux-x86_64/gs-922-linux-x86_64 /usr/bin && \
    cd .. && \
    rm -rf /tmp/gs && \
    mkdir /tmp/matlab_mcr && \
    cd /tmp/matlab_mcr/ && \
    wget http://ssd.mathworks.com/supportfiles/downloads/R2017a/deployment_files/R2017a/installers/glnxa64/MCR_R2017a_glnxa64_installer.zip && \
    unzip MCR_R2017a_glnxa64_installer.zip && \
    ./install -agreeToLicense yes -mode silent && \
    cd .. && \
    rm -rf /tmp/matlab_mcr && \
    mkdir /INPUTS && \
    mkdir /OUTPUTS && \
    mkdir /extra

# Copy binaries and other files
ADD extra /extra

# Set home directory in OUTPUTS folder
ENV HOME /OUTPUTS/.local/home

# Set tmp directory in OUTPUTS folder
ENV TMPDIR /OUTPUTS/.local/tmp
ENV TEMP /OUTPUTS/.local/tmp
ENV TMP /OUTPUTS/.local/tmp

# Set up LD_LIBRARY_PATH
ENV LD_LIBRARY_PATH "/lib/x86_64-linux-gnu:/usr/lib/x86_64-linux-gnu:/usr/local/cuda/lib64"

# Set environment for MATLAB MCR
ENV LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:/usr/local/MATLAB/MATLAB_Runtime/v92/runtime/glnxa64:/usr/local/MATLAB/MATLAB_Runtime/v92/bin/glnxa64:/usr/local/MATLAB/MATLAB_Runtime/v92/sys/os/glnxa64:/usr/local/MATLAB/MATLAB_Runtime/v92/sys/opengl/lib/glnxa64"

# Set up FSL
ENV PATH "${PATH}:/extra/fsl_5_0_10_eddy_5_0_11/bin"
ENV FSLDIR /extra/fsl_5_0_10_eddy_5_0_11

# Set CMD
RUN ln -sf bash /bin/sh
CMD /extra/pipeline.sh
