cmake_minimum_required(VERSION 3.16)
project(plasma0d LANGUAGES CXX)

# Plain C++17, standard library only, so a classmate can build it with any
# compiler without installing anything.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Off by default. Turn on if a built program crashes at start-up on Windows with two
# MinGW installs (the loader can pick the other install's libstdc++ DLL).
option(PLASMA_STATIC_RUNTIME "Link the C++ runtime statically (fixes a DLL mix-up on some Windows setups)" OFF)

add_executable(plasma0d src/main.cpp)
if(MSVC)
  target_compile_options(plasma0d PRIVATE /W4)
else()
  target_compile_options(plasma0d PRIVATE -Wall -Wextra)
  if(PLASMA_STATIC_RUNTIME)
    target_link_options(plasma0d PRIVATE -static)
  endif()
endif()

# ---- v0.2.0: 2D fluid model and its verification checks -----------------------
# The checks are plain programs that print PASS or FAIL; `ctest` runs them all
# (from the project folder, since they read config/).
function(plasma_target name src)
  add_executable(${name} ${src})
  if(MSVC)
    target_compile_options(${name} PRIVATE /W4)
  else()
    target_compile_options(${name} PRIVATE -Wall -Wextra)
    if(PLASMA_STATIC_RUNTIME)
      target_link_options(${name} PRIVATE -static)
    endif()
  endif()
endfunction()

plasma_target(plasma2d src/main2d.cpp)

enable_testing()
foreach(chk check1_reduction check2_diffusion check3_drift check4_sgsteady check5_poisson check6_mms)
  plasma_target(${chk} tests/${chk}.cpp)
  add_test(NAME ${chk} COMMAND ${chk} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endforeach()
