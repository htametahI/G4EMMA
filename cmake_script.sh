
rm Makefile
rm CMakeCache.txt
rm -r CMakeFiles
rm cmake_install.cmake
cmake \
  -DGeant4_DIR=/Users/mikeqiu/opt/geant4-11.4.0 \
  -DROOT_DIR=/opt/homebrew/Cellar/root/6.36.06_1 \
  -DCMAKE_POSITION_INDEPENDENT_CODE=OFF \
  -DCMAKE_EXE_LINKER_FLAGS="-no-pie" \
  /Users/mikeqiu/GEANT4/G4EMMA
