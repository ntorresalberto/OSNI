cmake_minimum_required(VERSION 3.5)

project(OSNI LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# external utilities to configure the package
include(GNUInstallDirs)
set(INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_INCLUDEDIR})
set(LIB_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR})


find_package(Eigen3 3.4 NO_MODULE REQUIRED)


include_directories(include)


add_library(${PROJECT_NAME}_Chebyshev include/chebyshev_differentiation.hpp
                                      src/chebyshev_differentiation.cpp)
target_link_libraries(${PROJECT_NAME}_Chebyshev Eigen3::Eigen)


add_library(${PROJECT_NAME}     include/ode_solver_interface.h
                                #src/ode_solver_interface.cpp
                                include/ode_a.h
                                src/ode_a.cpp
                                #include/ode_ab.h
                                #src/ode_ab.cpp
                                )
target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_Chebyshev Eigen3::Eigen)


#add_library(ForwardKinematics   include/quaternion_integrator.h
#                                src/quaternion_integrator.cpp
#                                #include/twist_integrator.h
#                                #src/twist_integrator.cpp
#                                )
#target_link_libraries(ForwardKinematics ${PROJECT_NAME} Eigen3::Eigen)


#add_executable(OSNI_test main.cpp)
#target_link_libraries(OSNI_test ${PROJECT_NAME}  ForwardKinematics)

# install headers
install(DIRECTORY include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})


install (TARGETS ${PROJECT_NAME}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)

# construct .cmake files
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}ConfigVersion.cmake
  VERSION 0.0.1
  COMPATIBILITY SameMajorVersion)

# generates OSNIConfig.cmake file to use our package in other projects
configure_package_config_file(
  Config.cmake.in # input template
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake # output config file
  INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake # where to put the config file during install
  PATH_VARS INCLUDE_INSTALL_DIR LIB_INSTALL_DIR # paths to be used
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

# when running make install, config files should be copied as well
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake
)
