@echo off
TITLE Vector Slicer installation

SET mypath=%~dp0
set pwd=%mypath:~0,-1%
IF EXIST "C:\src\vcpkg\" (
    cd C:\src\vcpkg\
    git pull
) else (
     mkdir C:\src
     cd C:\src
     git clone https://github.com/Microsoft/vcpkg.git
     cd vcpkg && bootstrap-vcpkg.bat
)

vcpkg install boost:x64-windows eigen3:x64-windows vtk:x64-windows
vcpkg integrate install

set cmake_path="C:/src/vcpkg/scripts/buildsystems/vcpkg.cmake"

echo Installing Morphoshell-OMP
cd %pwd%

IF EXIST ".git" (
    git pull
) else (
    echo Warning: Git directory does not exist, which will prevent you from easily updating the project in the future.
    echo Please consider cloning the GitHub repository, by going into the desired parent directory and running:
    echo git clone https://github.com/Zmmyslony/MorphoShell-OMP.git
)

cmake -S ./ -B ./build -DCMAKE_TOOLCHAIN_FILE=%cmake_path%
cmake --build ./build --config Release -j4

setx MORPHOSHELL "%pwd%\build\Release\Morphoshell.exe"

cmd /k