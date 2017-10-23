# This is the root ITK CMakeLists file.
cmake_minimum_required(VERSION 2.8.9)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif()

# This project is designed to be built outside the Insight source tree.
project(Registration)

# Find ITK.
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable(RegistrationRigid3D RegistrationRigid3D.cxx )

target_link_libraries(RegistrationRigid3D ${ITK_LIBRARIES})

#add_executable(RegistrationAffine3D RegistrationAffine3D.cxx )
#add_executable(RegistrationSimilarity3D RegistrationSimilarity3D.cxx )
#add_executable(RegistrationSimilarity3DQN RegistrationSimilarity3DQN.cxx )
#add_executable(RegistrationSimilarity3DGDLS RegistrationSimilarity3DGDLS.cxx )
#add_executable(RegistrationSimilarity3DLBFGSB RegistrationSimilarity3DLBFGSB.cxx )
#add_executable(PreRegistration PreRegistration.cxx)
add_executable(DeformableRegistrationBSpline3D DeformableRegistrationBSpline3D.cxx)
#add_executable(MultiResAffine3D MultiResAffine3D.cxx)
#add_executable(MultiResAffine3D2 MultiResAffine3D2.cxx)
add_executable(MultiStage MultiStage.cxx)
add_executable(MultiStage2 MultiStage2.cxx)
add_executable(MultiStage3 MultiStage3.cxx)

#target_link_libraries(RegistrationSimilarity3D ${ITK_LIBRARIES})
#target_link_libraries(RegistrationSimilarity3DQN ${ITK_LIBRARIES})
#target_link_libraries(RegistrationSimilarity3DGDLS ${ITK_LIBRARIES})
#target_link_libraries(RegistrationSimilarity3DLBFGSB ${ITK_LIBRARIES})
#target_link_libraries(RegistrationSimilarity3DLBFGSB ${ITK_LIBRARIES})
#target_link_libraries(RegistrationAffine3D ${ITK_LIBRARIES})
#target_link_libraries(PreRegistration ${ITK_LIBRARIES})
target_link_libraries(DeformableRegistrationBSpline3D ${ITK_LIBRARIES})
#target_link_libraries(MultiResAffine3D ${ITK_LIBRARIES})
#target_link_libraries(MultiResAffine3D2 ${ITK_LIBRARIES})
target_link_libraries(MultiStage ${ITK_LIBRARIES})
target_link_libraries(MultiStage2 ${ITK_LIBRARIES})
target_link_libraries(MultiStage3 ${ITK_LIBRARIES})

add_executable(testRigid testRigid.cxx)
#add_executable(testAffine testAffine.cxx)
#add_executable(testSimilarity testSimilarity.cxx)
add_executable(testBSpline testBSpline.cxx)
#add_executable(testMultiResAffine testMultiResAffine.cxx)
#add_executable(testMultiResAffine2 testMultiResAffine2.cxx)
add_executable(testMultiStage testMultiStage.cxx)
add_executable(testMulti2 testMulti2.cxx)
add_executable(testMulti3 testMulti3.cxx)