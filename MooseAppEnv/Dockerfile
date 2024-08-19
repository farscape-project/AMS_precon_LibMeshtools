#Build image from bas Ubuntu image
FROM ubuntu:24.04

#Allow MPI and OpenMP bindings in Docker
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

#Compile on multicore environment
ARG compile_cores=4

#Set the working directory of the solvers
ARG WORKDIR="opt"

#Install software dependencies
RUN apt-get update && \
    apt-get install -y \
    build-essential \
    libtirpc-dev \
    bison \
    flex \
    clang-format \
    clang-tidy \
    curl \
    doxygen \
    git \
    ninja-build \
    pre-commit \
    python3 \
    python3-dev \
    python3-distutils-extra \
    python-is-python3 \
    python3-pip \
    libboost-all-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libglew-dev \
    libglm-dev \
    libhdf5-mpich-dev \
    libnetcdf-dev \
    libomp-14-dev \
    libpng-dev \
    libsdl2-dev \
    libssl-dev \
    rsync \
    xxd

#Install python packages that should be from Pip
RUN pip install packaging --break-system-packages && \
	pip install pyyaml --break-system-packages && \
	pip install setuptools --break-system-packages && \
	pip install jinja2 --break-system-packages

#Install CMake
RUN cd /$WORKDIR && \
    git clone https://github.com/Kitware/CMake.git && \
    cd CMake/ && \
    git checkout 78ae759 && \
    ./bootstrap --parallel="$compile_cores" && \
    make -j"$compile_cores" && \
    make install && \
    alias cmake='/$WORKDIR/CMake/bin/cmake'


#Configure PETSc
RUN     cd /$WORKDIR && \
	mkdir PETSc && cd PETSc && \
	git clone https://github.com/petsc/petsc.git && \
	cd petsc && \
        git checkout main && \
	./configure \
	--prefix=/$WORKDIR/PETSc \
	CXXOPTFLAGS='-O3 -funroll-loops' COPTFLAGS='-O3 -funroll-loops' FOPTFLAGS='-O3' \
	--with-debugging=0 \
	--with-ssl=0 \
	--with-pic=1 \
	--with-openmp=1 \
	--with-mpi=1 \
	--with-shared-libraries=1 \
	--with-fortran-bindings=0 \
	--with-sowing=0 \
	--download-hypre=1 \
	--download-fblaslapack=1 \
	--download-metis=1 \
	--download-ptscotch=1 \
	--download-parmetis=1 \
	--download-superlu_dist=1\
	--download-scalapack=1 \
	--download-mumps=1 \
	--download-slepc=1 \
	--with-64-bit-indices=1 \
	PETSC_DIR=`pwd` PETSC_ARCH=arch-linux-c-opt && \
	make

#Build PETSc
RUN cd /$WORKDIR/PETSc/petsc && \
    make PETSC_DIR=/$WORKDIR/PETSc/petsc PETSC_ARCH=arch-linux-c-opt install

#Test Build
RUN cd /$WORKDIR/PETSc/petsc && \
    make PETSC_DIR=/$WORKDIR/PETSc/petsc PETSC_ARCH=arch-linux-c-opt test &&\
    cd .. && \
    cd .. && \
    export PETSC_DIR=/$WORKDIR/PETSc


#Install LibMesh, WASP and get/configure my own MOOSE
RUN export MOOSE_JOBS=$compile_cores && \
    cd /$WORKDIR && \
    git clone https://github.com/farscape-project/mooseNewPrecon.git && \
    mv mooseNewPrecon moose && cd moose && \
    git checkout next && \
    export PETSC_DIR=/$WORKDIR/PETSc && \
    export PETSC_ARCH=arch-linux-c-opt && \
    export CC=mpicc && \
    export CXX=mpicxx && \
    export F90=mpif90 && \
    export F77=mpif77 && \
    export FC=mpif90  && \
    METHODS=opt ./scripts/update_and_rebuild_libmesh.sh --with-mpi && \
    ./configure --with-derivative-size=200

RUN cd /$WORKDIR/moose/scripts && \
    METHODS=opt ./update_and_rebuild_wasp.sh

# Build MOOSE
RUN cd /$WORKDIR/moose/framework && \
    METHOD=opt make -j$compile_cores

RUN cd /$WORKDIR/moose/test && \
    METHOD=opt make -j$compile_cores

# This is needed or it mpiexec complains because docker runs as root
# Discussion on this issue https://github.com/open-mpi/ompi/issues/4451
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
