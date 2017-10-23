# This is the root ITK CMakeLists file.
cmake_minimum_required(VERSION 2.8.9)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif()


# This project is designed to be built outside the Insight source tree.
project(MasoumehCode)

# Find ITK.
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable( Multistage2D Multistage2D.cpp )

target_link_libraries(Multistage2D ${ITK_LIBRARIES})

add_executable( testMulti2D testMulti2D.cpp )

target_link_libraries(testMulti2D ${ITK_LIBRARIES})
