mkdir -p src
ln -s $LBM/src/*.f90 src/
ln -s $LBM/src/Makefile src/
ln -s $LBM/python/analysis.py .
cd src
make
mv lbm.x ../
cd ../
