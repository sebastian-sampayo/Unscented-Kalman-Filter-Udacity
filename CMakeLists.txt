project(UnscentedKF)

cmake_minimum_required (VERSION 3.5)

add_definitions(-std=c++0x)

set(sources
   src/ukf.cpp
   src/main.cpp
   src/tools.cpp)

add_executable(UnscentedKF ${sources})

# Add an option for debug version. Allows gdb debugging and defines macro DEBUG.
# If the option is not set, define a macro NDEBUG, so the asserts won't take effect.
option(DEFINE_DEBUG
  "Build the project using debugging code"
  OFF)
IF(DEFINE_DEBUG)
  message("Adding Debug flag...")
  add_definitions(-DDEBUG)
  set(CMAKE_BUILD_TYPE Debug)
  message("Build type is " ${CMAKE_BUILD_TYPE})
ELSE(DEFINE_DEBUG)
  message("Adding NotDebug flag...")
  add_definitions(-D NDEBUG)
ENDIF(DEFINE_DEBUG)

