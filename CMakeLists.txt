cmake_minimum_required(VERSION 3.20)
project(ShelloMorph)

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(CMAKE_CXX_STANDARD 17)

find_path(LIBCONFIG_DIR libconfig.h)
find_library(LIBCONFIGPP_LIB libconfig++)

set(SOURCE_FILES
        source/main.cpp source/advanceDynamics.cpp source/advanceDynamics.hpp source/calculations/calc_nonVertexPatchNodes_and_MatForPatchDerivs.cpp
        source/calculations/calc_nonVertexPatchNodes_and_MatForPatchDerivs.hpp source/calculations/calcCurvatures.cpp source/calculations/calcCurvatures.hpp source/calculations/calcDeformationForces.cpp
        source/calculations/calcDeformationForces.hpp source/calculations/calcEnergiesAndStresses.cpp source/calculations/calcEnergiesAndStresses.hpp source/calculations/calcNodeNeighbours.cpp
        source/calculations/calcNodeNeighbours.hpp source/calculations/calcNonDeformationForces_and_ImposeBCS.cpp source/calculations/calcNonDeformationForces_and_ImposeBCS.hpp
        source/calculations/calcSecFFsAndRelatedQuantities.cpp source/calculations/calcSecFFsAndRelatedQuantities.hpp source/calculations/calcTriangleAdjacencies_And_Edges.cpp
        source/calculations/calcTriangleAdjacencies_And_Edges.hpp source/calculations/calcTriangleGeometries_and_DialledProgTensors.cpp source/calculations/calcTriangleGeometries_and_DialledProgTensors.hpp
        source/calculations/calcTrianglesIncidentOnNodes.cpp source/calculations/calcTrianglesIncidentOnNodes.hpp source/CustomOutStreamClass.cpp source/CustomOutStreamClass.hpp source/Edge.cpp
        source/Edge.hpp source/equilCheck.cpp source/equilCheck.hpp source/functions/extract_Just_Filename.cpp source/functions/extract_Just_Filename.hpp source/functions/getRealTime.cpp
        source/functions/getRealTime.hpp source/initialDirAndFileHandling.cpp source/initialDirAndFileHandling.hpp source/functions/kahanSum.cpp source/functions/kahanSum.hpp source/Node.cpp
        source/Node.hpp source/functions/perturbInitialPositionsWithRandomNoise.cpp source/functions/perturbInitialPositionsWithRandomNoise.hpp source/readSettingsFile.cpp
        source/readSettingsFile.hpp source/readVTKData.cpp source/readVTKData.hpp source/setRemainingInitCond_and_NodeMasses.cpp source/setRemainingInitCond_and_NodeMasses.hpp
        source/SettingsStruct.hpp source/StatusEnum.hpp source/Triangle.cpp source/Triangle.hpp source/writeVTKDataOutput.cpp source/writeVTKDataOutput.hpp source/functions/zeroForces.cpp
        source/functions/zeroForces.hpp source/SettingsStruct.cpp)

add_executable(ShelloMorph ${SOURCE_FILES})

target_include_directories(ShelloMorph PRIVATE ${LIBCONFIG_DIR})
target_link_libraries(ShelloMorph PRIVATE ${LIBCONFIGPP_LIB})



