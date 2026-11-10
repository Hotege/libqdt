cmake_minimum_required(VERSION 3.16)

project(libqdt LANGUAGES C CXX)

include(GNUInstallDirs)

# Follow CMake's conventional BUILD_SHARED_LIBS switch. Static is the default.
option(BUILD_SHARED_LIBS "Build libqdt as a shared library" OFF)
option(QDT_BUILD_BENCHMARK "Build the qdt_bench benchmark executable" OFF)

find_package(Threads REQUIRED)

add_library(qdt
    libqdt/qdt.cpp
    libqdt/GuibasStolfiDwyer_Wrapper.cpp
    libqdt/GuibasStolfiDwyer_Threads.cpp
)
add_library(QDT::qdt ALIAS qdt)

set_target_properties(qdt PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    POSITION_INDEPENDENT_CODE ON
    EXPORT_NAME qdt
)

target_include_directories(qdt
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/libqdt>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(qdt PRIVATE Threads::Threads)

# qdt.h uses QDT_DLL to select __declspec(dllexport) on Windows.
# For static builds, the same define must also reach consumers so the header
# does not declare the API as __declspec(dllimport).
if(WIN32)
    if(BUILD_SHARED_LIBS)
        target_compile_definitions(qdt PRIVATE QDT_DLL)
    else()
        target_compile_definitions(qdt PUBLIC QDT_DLL)
    endif()
endif()

# The robust geometric predicates depend on strict IEEE-754 floating-point
# behavior. Do not compile this library with -ffast-math or /fp:fast.
if(MSVC)
    target_compile_options(qdt PRIVATE /fp:precise)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(qdt PRIVATE -fno-fast-math)
endif()

# Install the library and its public C API header.
install(TARGETS qdt
    EXPORT qdtTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(FILES libqdt/qdt.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Export a usable package: find_package(qdt CONFIG REQUIRED), then link QDT::qdt.
install(EXPORT qdtTargets
    FILE qdtTargets.cmake
    NAMESPACE QDT::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/qdt
)
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/qdtConfig.cmake" [=[
include(CMakeFindDependencyMacro)
find_dependency(Threads)
include("${CMAKE_CURRENT_LIST_DIR}/qdtTargets.cmake")
]=])
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/qdtConfig.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/qdt
)

if(QDT_BUILD_BENCHMARK)
    add_executable(qdt_bench test/qdt_bench.cpp)
    target_compile_features(qdt_bench PRIVATE cxx_std_14)
    target_link_libraries(qdt_bench PRIVATE QDT::qdt Threads::Threads)

    # MinGW does not honor the MSVC #pragma comment(lib, "psapi.lib").
    if(MINGW)
        target_link_libraries(qdt_bench PRIVATE psapi)
    endif()
endif()
