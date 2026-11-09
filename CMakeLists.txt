cmake_minimum_required(VERSION 3.16)
project(TokyoGoppaExperiments C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS ON)

function(add_experiment target source)
    add_executable(${target} ${source})
    target_compile_options(${target} PRIVATE -O2 -Wall -Wextra)
endfunction()

add_experiment(tokyo ymo.c)
add_experiment(bm_gf2m_lyon bm_gf2m/lyon.c)
add_experiment(daikoku_xeuclide daikoku/Xeuclide.c)
add_experiment(gf2m_peter_ptar gf2m-peter/p-tar.c)
add_experiment(marimo_vmo marimo/vmo.c)

if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(daikoku_xeuclide PRIVATE -mavx2)
    target_compile_options(marimo_vmo PRIVATE -mavx2)
endif()

find_package(OpenMP COMPONENTS C)
if(OpenMP_C_FOUND)
    target_link_libraries(daikoku_xeuclide PRIVATE OpenMP::OpenMP_C)
endif()