INCLUDE_PATH=  -fext-numeric-literals  -I. -I/home/ubuntu/openfpm/install/openfpm_numerics/include -I/home/ubuntu/openfpm/install/openfpm_pdata/include/config -I/home/ubuntu/openfpm/install/openfpm_pdata/include -I/home/ubuntu/openfpm/install/openfpm_data/include -I/home/ubuntu/openfpm/install/openfpm_vcluster/include -I/home/ubuntu/openfpm/install/openfpm_io/include -I/home/ubuntu/openfpm/install/openfpm_devices/include -I/home/ubuntu/openfpm/install/METIS/include -I/home/ubuntu/openfpm/install/PARMETIS/include -I/home/ubuntu/openfpm/install/BOOST/include -I/home/ubuntu/openfpm/install/HDF5/include -I/home/ubuntu/openfpm/install/LIBHILBERT/include   -I/home/ubuntu/openfpm/install/PETSC/include -I/home/ubuntu/openfpm/install/HYPRE/include -I/home/ubuntu/openfpm/install/MUMPS/include -I/home/ubuntu/openfpm/install/OPENBLAS/include -I/home/ubuntu/openfpm/install/SCALAPACK/include -I/home/ubuntu/openfpm/install/SUPERLU_DIST/include -I/home/ubuntu/openfpm/install/TRILINOS/include -I/home/ubuntu/openfpm/install/SUITESPARSE/include -I/home/ubuntu/openfpm/install/EIGEN
LIBS_PATH=  -fext-numeric-literals -L/home/ubuntu/openfpm/install/openfpm_devices/lib -L/home/ubuntu/openfpm/install/openfpm_pdata/lib  -L/home/ubuntu/openfpm/install/openfpm_vcluster/lib -L/home/ubuntu/openfpm/install/METIS/lib -L/home/ubuntu/openfpm/install/PARMETIS/lib  -L/home/ubuntu/openfpm/install/BOOST/lib -L/home/ubuntu/openfpm/install/HDF5/lib -L/home/ubuntu/openfpm/install/LIBHILBERT/lib   -L/home/ubuntu/openfpm/install/PETSC/lib -L/home/ubuntu/openfpm/install/HYPRE/lib -L/home/ubuntu/openfpm/install/MUMPS/lib -L/home/ubuntu/openfpm/install/OPENBLAS/lib -L/home/ubuntu/openfpm/install/SCALAPACK/lib -L/home/ubuntu/openfpm/install/SUPERLU_DIST/lib -L/home/ubuntu/openfpm/install/TRILINOS/lib -L/home/ubuntu/openfpm/install/SUITESPARSE/lib
LIBS=-lvcluster -lofpm_pdata -lofpmmemory -lparmetis -lmetis -lboost_iostreams -lhdf5 -llibhilbert   -lrt -lpetsc -lopenblas -lumfpack -lamd -lbtf -lcamd -lccolamd -lcholmod -lcolamd -lcxsparse -lklu -ldl -lrbio -lspqr -lsuitesparseconfig -lgomp
LIBS_SE2=-lvcluster -lofpmmemory_se2 -lparmetis -lmetis -lboost_iostreams -lhdf5 -llibhilbert  -lrt -lpetsc -lopenblas -lumfpack -lamd -lbtf -lcamd -lccolamd -lcholmod -lcolamd -lcxsparse -lklu -ldl -lrbio -lspqr -lsuitesparseconfig -lgomp
