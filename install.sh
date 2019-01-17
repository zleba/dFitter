wget -qO-  https://www.nikhef.nl/~h24/qcdnum-files/download/qcdnum170114.tar.gz | tar zxv
mkdir -f qcdnum
pwd=$PWD
cd qcdnum-17-01-14
./configure --prefix=$pwd/qcdnum
make -j`nproc`
make install
