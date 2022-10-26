cmake_minimum_required(VERSION 3.7...3.22)

################### PROJECT SETUP #################
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project("PSSolver" LANGUAGES CXX VERSION 0.0.1)

### Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
  message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file).
    Please create a separate build subdirectory.")
endif()

# Create a compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

################### PROJECT CONFIG ##############
# Set c++17 as the required standard, and do not allow
# compiler-specific extensions
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# set default build type
SET(DEFAULT_BUILD_TYPE "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}'")
  set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
      "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set(PROJECT_ROOT_NAME ${PROJECT_NAME})
set(PROJECT_ROOT_VERSION ${PROJECT_VERSION})
set(PROJECT_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
set(PROJECT_ROOT_BINARY_DIR ${CMAKE_BINARY_DIR})

# tell VS to export all symbols to its dll files
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  SET(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE CACHE BOOL "Export all symbols")
endif()

# When using clang, always use libc++
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libc++ -lc++abi")
endif()


################### DEPENDENCIES ################
# Download conan cmake config and use it to automatically fetch all required packages
if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
   message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
   file(DOWNLOAD "https://raw.githubusercontent.com/conan-io/cmake-conan/master/conan.cmake"
                  "${CMAKE_BINARY_DIR}/conan.cmake")
endif()

# run conan install and create cmake targets
include(${CMAKE_BINARY_DIR}/conan.cmake)
conan_cmake_run(CONANFILE conanfile.py
                BUILD missing
                BASIC_SETUP CMAKE_TARGETS)

# include cmake_paths generator file to get package-supplied PackageConfig.cmake
include(${CMAKE_BINARY_DIR}/conan_paths.cmake)


################### BUILDING ####################
include(GNUInstallDirs)

# get all sources and headers
file(GLOB CPP_FILES "${CMAKE_CURRENT_SOURCE_DIR}/lib/*.cpp")
file(GLOB_RECURSE HPP_FILES "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")

# build shared library
add_library(${PROJECT_NAME} SHARED ${CPP_FILES})
# Adding the install interface generator expression makes sure that cmake 
# chooses the correct include folder for current builds or already installed targets
target_include_directories(${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# set library postfix if it is a debug build
set_target_properties(${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "d")

# set target properties
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${HPP_FILES}")
set_target_properties(${PROJECT_NAME} PROPERTIES VERSION ${PROJECT_VERSION})
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)

# enable strict standard, extra warnings
target_compile_options(${PROJECT_NAME} PRIVATE
  $<$<CXX_COMPILER_ID:MSVC>:/W4 /WX>
  $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic -Werror>
)


################### INSTALLING ##################
# # set the correct paths for installation
set(LOCAL_CONFIG_DIR "lib/cmake/${PROJECT_NAME}")

# Adding the install interface generator expression makes sure that the include
# files are installed to the proper location (provided by GNUInstallDirs)
include(CMakePackageConfigHelpers)
write_basic_package_version_file("${PROJECT_NAME}ConfigVersion.cmake"
                                 VERSION ${PROJECT_VERSION}
                               COMPATIBILITY AnyNewerVersion)

# automatically generate config files
configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION
  ${LOCAL_CONFIG_DIR})

# install library target to lib location and headers to include location
install(TARGETS ${PROJECT_NAME}
        EXPORT ${PROJECT_NAME}_Targets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")

# install exported targets
install(EXPORT ${PROJECT_NAME}_Targets
        FILE ${PROJECT_NAME}Targets.cmake
        DESTINATION ${LOCAL_CONFIG_DIR})

# Install the generated config files
install(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
              "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION ${LOCAL_CONFIG_DIR})


################### WRAPPING ####################
# wrap library for use in different languages
option(BUILD_PYTHON "Build python module" OFF)
if(BUILD_PYTHON)
  add_subdirectory(wrapping)
endif()

################## DOCUMENTATION ################
option(BUILD_DOCS "Build documentation" OFF)
if(BUILD_DOCS)
  add_subdirectory(docs)
endif()


find_package(Eigen3 REQUIRED)
if(Eigen3_FOUND)
  target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)
endif()


