cmake_minimum_required(VERSION 3.20)
project(sgk_quadric CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Locate the vendor SGGK 1.7 SDK via this project's own find-module
# (order: -DSGGK_SDK_DIR= -> env SGGK_ROOT -> <project>/../SGGK -> third_party/SGGK -> ./SGGK;
#  see README). Outputs: SGGK_SDK_DIR / SGGK_INCLUDE_DIR / SGGK_LIB_DIR / SGGK_BIN_DIR / SGGK_LICENSE_FILE
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(SGGK MODULE REQUIRED)

# Stage the SDK runtime once: copy the DLLs and sggk.lic into the build output root.
# Doing this per-target in POST_BUILD made parallel Ninja builds race on the same directory
# (copy errors / files in use), so it runs exactly once and every target depends on it.
# Applies to single-config generators (build.bat pins Ninja); for multi-config, stage per config.
# sggk.lic is staged only when it sits at the SDK root (<SDK>/sggk.lic), the layout the SDK ships;
# a licence kept anywhere else must be copied next to the executables by hand.
set(SGGK_STAGE_STAMP "${CMAKE_BINARY_DIR}/.sggk_runtime_staged")
set(_stage_cmds
    COMMAND ${CMAKE_COMMAND} -E copy_directory "${SGGK_BIN_DIR}" "${CMAKE_BINARY_DIR}"
)
if(EXISTS "${SGGK_SDK_DIR}/sggk.lic")
    list(APPEND _stage_cmds
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${SGGK_SDK_DIR}/sggk.lic" "${CMAKE_BINARY_DIR}/sggk.lic"
    )
endif()
list(APPEND _stage_cmds COMMAND ${CMAKE_COMMAND} -E touch "${SGGK_STAGE_STAMP}")

add_custom_command(
    OUTPUT  "${SGGK_STAGE_STAMP}"
    ${_stage_cmds}
    COMMENT "Stage SGGK runtime DLLs and sggk.lic (once)"
    VERBATIM
)
add_custom_target(stage_runtime DEPENDS "${SGGK_STAGE_STAMP}")

# All 23 import libraries, listed explicitly (no globbing).
set(SGGK_LIBS
    SGGK_BSplineAlgo
    SGGK_Boolean
    SGGK_DataBase
    SGGK_DataExchangeBase
    SGGK_DxfReader
    SGGK_Foundation
    SGGK_GeomAlgo
    SGGK_GeomBase
    SGGK_GeomInt
    SGGK_GeomProject
    SGGK_Geometry
    SGGK_Heal
    SGGK_IgesExchange
    SGGK_InteropUtil
    SGGK_Math
    SGGK_ModelAnalysis
    SGGK_ModelingBase
    SGGK_ModelingPrim
    SGGK_Offset
    SGGK_Projection
    SGGK_ShrinkWrap
    SGGK_StepExchange
    SGGK_Topology
)

# These compile switches are hard requirements, not preferences: /utf-8 (the SDK headers are UTF-8
# with BOM while a zh-CN locale is CP936, which otherwise yields C4819 plus syntax errors), /MD (the
# SDK exports embed MSVC STL types, so a static CRT would give this program and the SDK separate
# heaps), c++17 and /EHsc.
set(SGKQ_COMPILE_OPTIONS $<$<CXX_COMPILER_ID:MSVC>:/utf-8 /EHsc /MD>)

# --- sgkqcore: numeric core + the ingest/detect/fit/pipeline modules.
# The math core uses SDK geometry types only; the rest uses the SDK surface/topology API.
add_library(sgkqcore STATIC
    src/math/math.cpp
    src/io/io.cpp
    src/detect/detect.cpp
    src/fit/fit.cpp
    src/pipeline/pipeline.cpp
)

target_include_directories(sgkqcore PUBLIC
    "${SGGK_INCLUDE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)
set_target_properties(sgkqcore PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)
target_compile_options(sgkqcore PRIVATE ${SGKQ_COMPILE_OPTIONS})
target_link_directories(sgkqcore PUBLIC "${SGGK_LIB_DIR}")
target_link_libraries(sgkqcore PUBLIC ${SGGK_LIBS})

# --- sgkqtestkit: evaluation-harness support (corpus manifest parser + executable-relative data dir).
# Kept outside the library so the shipped ingestion API (io.hpp) covers only real CAD data.
# Only the corpus/benchmark tools link it.
add_library(sgkqtestkit STATIC src/testkit/testkit.cpp)

target_include_directories(sgkqtestkit PUBLIC
    "${SGGK_INCLUDE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)
set_target_properties(sgkqtestkit PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)
target_compile_options(sgkqtestkit PRIVATE ${SGKQ_COMPILE_OPTIONS})
target_link_directories(sgkqtestkit PUBLIC "${SGGK_LIB_DIR}")
target_link_libraries(sgkqtestkit PUBLIC sgkqcore ${SGGK_LIBS})

# --- Target helpers. One per source layout, so a target's directory is derived rather than repeated:
#   apps/probes/<target>/   kernel-only programs (no project headers, no sgkqcore)
#   apps/checks/<target>/   known-answer self-checks
#   apps/harness/<target>/  corpus / benchmark programs (adds sgkqtestkit)
#   tools/<target>/         corpus generator (adds sgkqtestkit)
# The hard switches, include dirs, link dirs and the stage_runtime dependency are applied here, once.
function(sgkq_configure_app tgt src with_harness)
    add_executable(${tgt} "${src}")
    target_include_directories(${tgt} PRIVATE
        "${SGGK_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
    )
    set_target_properties(${tgt} PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
    )
    target_compile_options(${tgt} PRIVATE ${SGKQ_COMPILE_OPTIONS})
    target_link_directories(${tgt} PRIVATE "${SGGK_LIB_DIR}")
    if(with_harness)
        target_link_libraries(${tgt} PRIVATE sgkqcore sgkqtestkit ${SGGK_LIBS})
    else()
        target_link_libraries(${tgt} PRIVATE sgkqcore ${SGGK_LIBS})
    endif()
    add_dependencies(${tgt} stage_runtime)
endfunction()

function(sgkq_add_probe tgt)
    add_executable(${tgt} "apps/probes/${tgt}/main.cpp")
    target_include_directories(${tgt} PRIVATE "${SGGK_INCLUDE_DIR}")
    set_target_properties(${tgt} PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
    )
    target_compile_options(${tgt} PRIVATE ${SGKQ_COMPILE_OPTIONS})
    target_link_directories(${tgt} PRIVATE "${SGGK_LIB_DIR}")
    target_link_libraries(${tgt} PRIVATE ${SGGK_LIBS})
    add_dependencies(${tgt} stage_runtime)
endfunction()

function(sgkq_add_check tgt)
    sgkq_configure_app(${tgt} "apps/checks/${tgt}/main.cpp" FALSE)
endfunction()

function(sgkq_add_harness tgt)
    sgkq_configure_app(${tgt} "apps/harness/${tgt}/main.cpp" TRUE)
endfunction()

function(sgkq_add_tool_harness tgt)
    sgkq_configure_app(${tgt} "tools/${tgt}/main.cpp" TRUE)
endfunction()

# --- Targets: 3 probes + 7 checks + 3 harness tools + 1 tool = 14 executables, plus 2 static libs ---
foreach(_tgt IN ITEMS hello probe_writeback probe_speed)
    sgkq_add_probe(${_tgt})
endforeach()

foreach(_tgt IN ITEMS
        test_math
        test_io
        test_detect
        test_fit
        test_pipeline
        test_exceptions
        test_g1_tangent)
    sgkq_add_check(${_tgt})
endforeach()

foreach(_tgt IN ITEMS run_corpus bench_perf batch_stability)
    sgkq_add_harness(${_tgt})
endforeach()

sgkq_add_tool_harness(gen_corpus)
