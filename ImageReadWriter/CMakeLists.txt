# This is the root ITK CMakeLists file.
cmake_minimum_required(VERSION 2.8.9)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif()


# This project is designed to be built outside the Insight source tree.
project(ImageReaderWriter)

# Find ITK.
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable( ImageRW ImageRW.cxx )

target_link_libraries(ImageRW ${ITK_LIBRARIES})

add_executable( DICOM2Volume DICOM2Volume.cxx )

target_link_libraries(DICOM2Volume ${ITK_LIBRARIES})

add_executable( JPEG2Volume JPEG2Volume.cxx )

target_link_libraries(JPEG2Volume ${ITK_LIBRARIES})

add_executable( ReorientAtlas ReorientAtlas.cxx)

target_link_libraries(ReorientAtlas ${ITK_LIBRARIES})

add_executable( CalibrationReader CalibrationReader.cxx)

target_link_libraries( CalibrationReader ${ITK_LIBRARIES})

add_executable( CombineLabelMap CombineLabelMap.cxx)

target_link_libraries( CombineLabelMap ${ITK_LIBRARIES})

