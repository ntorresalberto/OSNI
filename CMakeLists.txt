cmake_minimum_required(VERSION 3.5)

project(OSNI LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


find_package(Eigen3 3.4 NO_MODULE REQUIRED)


include_directories(include)


add_library(${PROJECT_NAME}_Chebyshev include/chebyshev_differentiation.hpp
                                      src/chebyshev_differentiation.cpp)
target_link_libraries(${PROJECT_NAME}_Chebyshev Eigen3::Eigen)


add_library(${PROJECT_NAME}     include/spectral_integration_base.h
                                src/spectral_integration_base.cpp
                                include/ode_a.h
                                src/ode_a.cpp
                                include/ode_ab.h
                                src/ode_ab.cpp
                                )
target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_Chebyshev Eigen3::Eigen)


add_library(ForwardKinematics   include/quaternion_integrator.h
                                src/quaternion_integrator.cpp
                                #include/twist_integrator.h
                                #src/twist_integrator.cpp
                                )
target_link_libraries(ForwardKinematics ${PROJECT_NAME} Eigen3::Eigen)


add_executable(OSNI_test main.cpp)
target_link_libraries(OSNI_test ${PROJECT_NAME}  ForwardKinematics)
