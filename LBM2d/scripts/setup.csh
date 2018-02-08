mkdir -p src
ln -s /home/noemi/Documents/LBM/LBM2d/src/*.f90 src/
ln -s /home/noemi/Documents/LBM/LBM2d/src/Makefile src/
cd src
make
mv lbm.x ../
cd ../
