cmake_minimum_required(VERSION 3.16)
project(Morphoshell)

find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED COMPONENTS system filesystem)
include_directories(${Boost_INCLUDE_DIRS})

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(VTK COMPONENTS)

set(CMAKE_CXX_STANDARD 17)

find_path(LIBCONFIG_DIR libconfig.h REQUIRED)
find_library(LIBCONFIGPP_LIB libconfig++)
if (NOT LIBCONFIGPP_LIB)
    find_library(LIBCONFIGPP_LIB config++)
endif ()

set(SOURCE_FILES
        source/main.cpp source/advanceDynamics.cpp source/advanceDynamics.hpp source/calculations/calc_nonVertexPatchNodes_and_MatForPatchDerivs.cpp
        source/calculations/calc_nonVertexPatchNodes_and_MatForPatchDerivs.hpp source/calculations/calcCurvatures.cpp source/calculations/calcCurvatures.hpp source/calculations/calcDeformationForces.cpp
        source/calculations/calcDeformationForces.hpp source/calculations/calcEnergiesAndStresses.cpp source/calculations/calcEnergiesAndStresses.hpp source/calculations/calcNodeNeighbours.cpp
        source/calculations/calcNodeNeighbours.hpp
        source/calculations/calcTriangleAdjacencies_And_Edges.cpp
        source/calculations/calcTriangleAdjacencies_And_Edges.hpp source/calculations/calcTriangleGeometries_and_DialledProgTensors.cpp source/calculations/calcTriangleGeometries_and_DialledProgTensors.hpp
        source/calculations/calcTrianglesIncidentOnNodes.cpp source/calculations/calcTrianglesIncidentOnNodes.hpp source/Edge.cpp
        source/Edge.hpp source/EquilibriumCheck.cpp source/EquilibriumCheck.hpp source/functions/extract_Just_Filename.cpp source/functions/extract_Just_Filename.hpp source/functions/getRealTime.cpp
        source/functions/getRealTime.hpp source/initialDirAndFileHandling.cpp source/initialDirAndFileHandling.hpp source/functions/kahanSum.cpp source/functions/kahanSum.hpp source/Node.cpp
        source/Node.hpp source/functions/perturbInitialPositionsWithRandomNoise.cpp source/functions/perturbInitialPositionsWithRandomNoise.hpp
        source/readVTKData.cpp source/readVTKData.hpp source/setRemainingInitCond_and_NodeMasses.cpp source/setRemainingInitCond_and_NodeMasses.hpp
        source/SimulationStatus.hpp source/Triangle.cpp source/Triangle.hpp source/writeVTKDataOutput.cpp source/writeVTKDataOutput.hpp source/functions/zeroForces.cpp
        source/functions/zeroForces.hpp
        source/simulation.cpp
        source/simulation.h
        source/configuration/gravity_config.cpp
        source/configuration/gravity_config.h
        source/configuration/config_base.cpp
        source/configuration/config_base.h
        source/configuration/core_config.cpp
        source/configuration/core_config.h
        source/settings_new.cpp
        source/settings_new.h
        source/configuration/slide_configuration.cpp
        source/configuration/slide_configuration.h
        source/physics/slide.cpp
        source/physics/slide.h
        source/configuration/derived_configuration.cpp
        source/configuration/derived_configuration.h
        source/physics/cone.cpp
        source/physics/cone.h
        source/physics/rigid_body.cpp
        source/physics/rigid_body.h)

add_executable(Morphoshell ${SOURCE_FILES})

target_link_libraries(Morphoshell PRIVATE ${Boost_LIBRARIES})
target_include_directories(Morphoshell PRIVATE ${LIBCONFIG_DIR})
target_link_libraries(Morphoshell PRIVATE ${LIBCONFIGPP_LIB})
target_link_libraries(Morphoshell PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(Morphoshell PRIVATE Eigen3::Eigen)
target_link_libraries(Morphoshell PRIVATE ${VTK_LIBRARIES})
#
#vtk_module_autoinit(
#        TARGETS Morphoshell
#        MODULES ${VTK_LIBRARIES}
#)

