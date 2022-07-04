cmake_minimum_required(VERSION 3.5)

project(osni_test LANGUAGES CXX)

find_package(OSNI REQUIRED)
# find_package(Boost REQUIRED COMPONENTS math)
find_package(Boost REQUIRED)
find_package(Eigen3 3.4 NO_MODULE REQUIRED)


message("- CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
message("- CMAKE_INSTALL_LIBDIR: ${CMAKE_INSTALL_LIBDIR}")
message("- CMAKE_INSTALL_INCLUDEDIR: ${CMAKE_INSTALL_INCLUDEDIR}")
message("- OSNI_INCLUDE_DIRS: ${OSNI_INCLUDE_DIRS}")
message("- OSNI_LIBRARIES: ${OSNI_LIBRARIES}")
message("- Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
message("- Boost_LIBRARIES: ${Boost_LIBRARIES}")


include_directories(
  ${Boost_INCLUDE_DIRS}
  )

add_executable(osni_test main.cpp)
# target_link_libraries(osni_test OSNI Eigen3::Eigen ${Boost_LIBRARIES})
target_link_libraries(osni_test OSNI Eigen3::Eigen Boost::boost)
