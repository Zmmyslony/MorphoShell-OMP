cmake_minimum_required(VERSION 3.28)
project(Morphoshell)

set(CMAKE_CXX_STANDARD 17)

# Try to use ccache if available to speed up repeated builds
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    message(STATUS "ccache found: ${CCACHE_PROGRAM} - enabling compiler launcher")
    set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
endif()

add_compile_options("$<$<CONFIG:RELEASE>:-fno-math-errno;-march=native;-mfma>")

find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED CONFIG COMPONENTS system filesystem iostreams)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(VTK REQUIRED)

set(SOURCE_FILES
    source/main.cpp
    source/calculations/calc_nonVertexPatchNodes_and_MatForPatchDerivs.cpp
    source/calculations/calcCurvatures.cpp
    source/calculations/calcDeformationForces.cpp
    source/calculations/calcEnergiesAndStresses.cpp
    source/calculations/calcNodeNeighbours.cpp
    source/calculations/calcTriangleAdjacencies_And_Edges.cpp
    source/calculations/calcTrianglesIncidentOnNodes.cpp
    source/Edge.cpp
    source/functions/extract_Just_Filename.cpp
    source/functions/getRealTime.cpp
    source/initialDirAndFileHandling.cpp
    source/functions/kahanSum.cpp
    source/Node.cpp
    source/functions/perturbInitialPositionsWithRandomNoise.cpp
    source/readVtk.cpp
    source/setRemainingInitCond_and_NodeMasses.cpp
    source/SimulationStatus.hpp
    source/Triangle.cpp
    source/exportVtk.cpp
    source/simulation.cpp
    source/configuration/gravity_config.cpp
    source/configuration/config_base.cpp
    source/configuration/core_config.cpp
    source/settings.cpp
    source/configuration/slide_configuration.cpp
    source/physics/slide.cpp
    source/physics/cone.cpp
    source/physics/rigid_body.cpp
    source/physics/magnetic_field.cpp
)

add_executable(Morphoshell ${SOURCE_FILES})

target_include_directories(Morphoshell PRIVATE ${Boost_INCLUDE_DIRS})

# Use precompiled headers to reduce per-translation-unit overhead for heavy headers
# Only enable PCH for compilers with reliable support on this platform
if(DEFINED CMAKE_CXX_COMPILER_ID)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        target_precompile_headers(Morphoshell PRIVATE "source/pch.hpp")
    else()
        message(STATUS "Skipping PCH: unsupported or problematic compiler (${CMAKE_CXX_COMPILER_ID}) on this platform")
    endif()
endif()

if(MSVC)
    target_compile_options(Morphoshell PRIVATE "-openmp:llvm")
endif()

target_link_libraries(Morphoshell PRIVATE
    Boost::filesystem
    Boost::system
    Boost::iostreams
    OpenMP::OpenMP_CXX
    Eigen3::Eigen
    ${VTK_LIBRARIES}
)

if(WIN32 AND CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    # MinGW usually uses standard math libs or specific MinGW-MKL ports
else()
    find_package(MKL)
    if(MKL_FOUND)
        target_link_libraries(Morphoshell PRIVATE MKL::MKL)
    endif()
endif()

vtk_module_autoinit(
    TARGETS Morphoshell
    MODULES ${VTK_LIBRARIES}
)