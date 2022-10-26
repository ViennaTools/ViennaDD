cmake_minimum_required(VERSION 3.7...3.22)

project("${PROJECT_ROOT_NAME}Python"
        VERSION ${PROJECT_ROOT_VERSION})

# include library
set(PYTHON_MODULE_SOURCE "PythonWrapping.cpp")
set(PYTHON_MODULE_LIBRARIES ${PROJECT_ROOT_NAME})
set(PYTHON_MODULE_NAME "${PROJECT_ROOT_NAME}")
set(PYTHON_PYBIND "${PROJECT_NAME}")

message(STATUS "Building python module: ${PYTHON_MODULE_NAME}")

# find dependencies, first one is our pre-built shared library
find_package(${PROJECT_ROOT_NAME} REQUIRED PATHS ${PROJECT_ROOT_BINARY_DIR} NO_DEFAULT_PATH)
find_package(pybind11 REQUIRED)

# build python module
pybind11_add_module(${PYTHON_PYBIND} ${PYTHON_MODULE_SOURCE})
target_link_libraries(${PYTHON_PYBIND} PRIVATE ${PYTHON_MODULE_LIBRARIES})
# provide module name to the compiler
target_compile_definitions(${PYTHON_PYBIND} PRIVATE -DMODULE_NAME=${PYTHON_MODULE_NAME})
# set the output name of the python module
set_target_properties(${PYTHON_PYBIND} PROPERTIES OUTPUT_NAME ${PYTHON_MODULE_NAME})

# set rpath so dynamic loading works as expected
set_target_properties(${PYTHON_PYBIND} PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)


############################# PACKAGING #############################
# Build python wheel package from pybind library
set(PYTHON_MODULE_INDIR "${CMAKE_CURRENT_LIST_DIR}/PythonLib")
set(PYTHON_MODULE_OUTDIR "${CMAKE_BINARY_DIR}/wheel/${PROJECT_ROOT_NAME}")
set(PYTHON_SETUP_INFILE "${CMAKE_CURRENT_LIST_DIR}/setup.py")
set(PYTHON_SETUP_OUTFILE "${CMAKE_BINARY_DIR}/wheel/setup.py")
set(PYTHON_INIT_INFILE "${PYTHON_MODULE_INDIR}/__init__.py")
set(PYTHON_INIT_OUTFILE "${PYTHON_MODULE_OUTDIR}/__init__.py")

# read README to put as docstring for the python module
file(READ "${CMAKE_SOURCE_DIR}/README.md" README_TEXT)

# Setup the python configuration for our module
configure_file(${PYTHON_SETUP_INFILE} ${PYTHON_SETUP_OUTFILE} @ONLY)
# Setup up __init__ file for python package
configure_file(${PYTHON_INIT_INFILE} ${PYTHON_INIT_OUTFILE} @ONLY)
configure_file(${CMAKE_CURRENT_LIST_DIR}/pyproject.toml  ${CMAKE_BINARY_DIR}/wheel/pyproject.toml COPYONLY)

# copy all shared libraries to module path
add_custom_command(TARGET ${PYTHON_PYBIND} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/* ${PYTHON_MODULE_OUTDIR}
)
# # copy all python files to module path
# configure_file(${CMAKE_CURRENT_LIST_DIR}/PythonLib/IntegrableFunctions.py  ${PYTHON_MODULE_OUTDIR} COPYONLY)

# Add command to build the module
set(PYTHON_WHEEL "${PROJECT_NAME}Wheel")
add_custom_target(${PYTHON_WHEEL}
        COMMAND "${PYTHON_EXECUTABLE}" -m pip wheel .
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/wheel"
        COMMENT "Build python wheel package"
        )
add_dependencies(${PYTHON_WHEEL} ${PYTHON_PYBIND})


# # install python module in library destination
# install(TARGETS ${PYTHON_PYBIND}
#   LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

