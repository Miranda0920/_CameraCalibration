cmake_minimum_required(VERSION 3.12)
project(CameraCalibration)

set(CMAKE_CXX_STANDARD 11)
find_package(Ceres REQUIRED)
find_package(OpenCV 3.0 QUIET)

include_directories(${OpenCV_INCLUDE_DIRS} ${CERES_INCLUDE_DIRS}
        ${PROJECT_SOURCE_DIR}/include)
add_executable(CameraCalibration src/main.cpp
        include/calibration.h
        src/calibration.cpp
        include/matrixSolve.h
        src/matrixSolve.cpp
        src/CeresOptimize.cpp
        include/CeresOptimize.h)
target_link_libraries(CameraCalibration ${OpenCV_LIBS} ${CERES_LIBRARIES})