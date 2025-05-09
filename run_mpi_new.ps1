g++ $args[1] -I $env:MSMPI_INC\ -L $env:MSMPI_LIB64\ -lmsmpi -o Build/mpiRun
mpiexec -n $args[0] Build/mpiRun