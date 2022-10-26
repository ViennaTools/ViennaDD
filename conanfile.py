from conans import ConanFile, CMake, tools
import os
import fnmatch
import shutil

class ConanSetup(ConanFile):
    name = "PSSolver"
    version = "0.1.0"
    license = "GNU Lesser General Public License version 3 (LGPL-3.0)"
    author = "Xaver klemenschits@iue.tuwien.ac.at, Markus kampl@iue.tuwien.ac.at"
    url = ""
    description = "IntegrationLib is a example for demonstrating interfaces."
    topics = ("Example", "Template", "Project")
    settings = "os", "compiler", "build_type", "arch"
    # options = {"shared": [True, False], "fPIC": [True, False]}
    # default_options = {"shared": False, "fPIC": True}
    generators = "cmake", "cmake_find_package", "cmake_paths"

    exports = "LICENSE"
    exports_sources = ("include/*", "CMakeLists.txt" , "cmake/*", "lib/*")
    _build_folder = "build/"

    requires = "doxygen/1.9.2", "pybind11/2.9.1", "eigen/3.4.0"

    def _configure_cmake(self):
        cmake = CMake(self)
        cmake.definitions["BUILD_TESTING"] = "OFF"
        cmake.configure(build_folder=self._build_folder)
        return cmake

    def _prep_build_folder(self):
        """
        Our CMakeLists.txt enforces out of source builds for safety reasons so we need to first create the 
        build folder and then copy the Find*.cmake files into it.
        """
        buildRoot = self.build_folder
        build = os.path.join(buildRoot, self._build_folder)
        if os.path.exists(build):
            shutil.rmtree(build)
        os.mkdir(build)
        shutil.copy2("conanbuildinfo.cmake", build)
        shutil.copy2("conan_paths.cmake", build)
        with os.scandir(buildRoot) as entries:
            for entry in entries:
                if entry.is_file() and fnmatch.fnmatch(entry.name, "Find*.cmake"):
                    shutil.copy2(entry.path, build)

    def build(self):
        # copy conan generated files into build folder
        self._prep_build_folder()

        cmake = self._configure_cmake()
        cmake.build()

    def package(self):
        self.copy(pattern="LICENSE", dst="licenses")
        cmake = self._configure_cmake()
        cmake.install()

    def package_info(self):
        lib_suffix = "d" if self.settings.build_type == "Debug" else ""

        self.cpp_info.names["cmake"] = self.name
        self.cpp_info.names["cmake_paths"] = self.name
        self.cpp_info.names["cmake_find_package"] = self.name
        self.cpp_info.names["cmake_find_package_multi"] = self.name

        self.cpp_info.includedirs = ["include"]
        self.cpp_info.libs = [self.name + lib_suffix]
