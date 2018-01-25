mkdir -p src
ln -s $LBM/src/*.f90 src/
ln -s $LBM/src/Makefile src/
cd src
make
mv lbm.x ../
cd ../
