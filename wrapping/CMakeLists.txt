cmake_minimum_required(VERSION 3.7...3.22)

option(BUILD_PYTHON_MODULE "Build documentation" ON)
if(BUILD_PYTHON_MODULE)
  add_subdirectory(Python)
endif()
