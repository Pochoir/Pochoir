#!/bin/bash

set -x
cp -rf /home/yuantang/Git/Pochoir/Release/Examples/*.cpp Examples/.
cp -rf /home/yuantang/Git/Pochoir/Release/Examples/*.hpp Examples/.
cp -rf /home/yuantang/Git/Pochoir/Release/Examples/*.h Examples/.
cp -rf /home/yuantang/Git/Pochoir/Release/Examples/Makefile Examples/.
cp -rf /home/yuantang/Git/Pochoir/Release/Examples/*.sh Examples/.
cp -rf /home/yuantang/Git/Pochoir/Release/Examples/LBM_pochoir2 ./Examples/LBM_pochoir2
cp /home/yuantang/Git/Pochoir/Release/run.sh .
cp /home/yuantang/Git/Pochoir/Release/pochoir .
cp /home/yuantang/Git/Pochoir/Release/*.hpp .
cp /home/yuantang/Git/Pochoir/Release/*.cpp .
cp /home/yuantang/Git/Pochoir/Release/*.hs .
cp /home/yuantang/Git/Pochoir/Release/Makefile .
cp /home/yuantang/Git/Pochoir/Release/LICENSE .
set +x
