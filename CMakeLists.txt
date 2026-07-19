cmake_minimum_required(VERSION 3.16)
project(LVN LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Fetch ImGui
include(FetchContent)
FetchContent_Declare(imgui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG v1.91.5
    GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(imgui)

set(IMGUI_SOURCES
    ${imgui_SOURCE_DIR}/imgui.cpp
    ${imgui_SOURCE_DIR}/imgui_draw.cpp
    ${imgui_SOURCE_DIR}/imgui_tables.cpp
    ${imgui_SOURCE_DIR}/imgui_widgets.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_sdl2.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp)

# Fetch ImPlot
FetchContent_Declare(implot
    GIT_REPOSITORY https://github.com/epezent/implot.git
    GIT_TAG v0.16
    GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(implot)

set(IMPLOT_SOURCES
    ${implot_SOURCE_DIR}/implot.cpp
    ${implot_SOURCE_DIR}/implot_items.cpp)

# Fetch Eigen
FetchContent_Declare(eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
    GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(eigen)

# Find SDL2 and OpenGL
find_package(SDL2 REQUIRED)
if(APPLE)
    find_library(OPENGL_FRAMEWORK OpenGL REQUIRED)
    set(OPENGL_LIB ${OPENGL_FRAMEWORK})
else()
    find_package(OpenGL REQUIRED)
    set(OPENGL_LIB OpenGL::GL)
endif()

file(GLOB_RECURSE SOURCES src/*.cpp)

add_executable(${PROJECT_NAME} ${SOURCES} ${IMGUI_SOURCES} ${IMPLOT_SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE
    src
    ${imgui_SOURCE_DIR}
    ${imgui_SOURCE_DIR}/backends
    ${implot_SOURCE_DIR}
    ${eigen_SOURCE_DIR}
)

target_link_libraries(${PROJECT_NAME} PRIVATE
    $<TARGET_NAME_IF_EXISTS:SDL2::SDL2main>
    SDL2::SDL2
    ${OPENGL_LIB}
)

target_compile_definitions(${PROJECT_NAME} PRIVATE GL_SILENCE_DEPRECATION)

if(MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE /W3 /GR- /EHs-c-)
else()
    target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -O3 -fno-rtti -fno-exceptions)
endif()
