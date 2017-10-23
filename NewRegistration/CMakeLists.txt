# This is the root ITK CMakeLists file.
cmake_minimum_required(VERSION 2.8.9)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif()


# This project is designed to be built outside the Insight source tree.
project(NewRegistration)

# Find ITK.
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable( Stage4Registration Stage4Registration.cpp )

target_link_libraries(Stage4Registration ${ITK_LIBRARIES})

add_executable( testStage4 testStage4.cpp )

target_link_libraries(testStage4 ${ITK_LIBRARIES})

add_executable( Stage4Bspline Stage4Bspline.cpp )

target_link_libraries(Stage4Bspline ${ITK_LIBRARIES})

add_executable( testStage4Bspline testStage4Bspline.cpp )

target_link_libraries(testStage4Bspline ${ITK_LIBRARIES})
