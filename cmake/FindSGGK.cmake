# =============================================================================================
# FindSGGK.cmake - locate the prebuilt binary "SGGK" SDK (China Shipbuilding HaiZhou).
#
# Why this module exists:
#   SGGK ships only as headers + MSVC import libraries + DLLs + a license file. The vendor
#   provides no CMake package config (no SGGKConfig.cmake) and the material is not
#   redistributable, so this project carries its own find-module following the CMake
#   Find<Package> convention - the standard approach for a binary-only dependency.
#
# A candidate directory qualifies as the SDK root only if it contains BOTH:
#     <root>/include/Foundation/SGGK_Version.h
#     <root>/x64-win/lib/SGGK_Foundation.lib
#
# Search order (first hit wins):
#     1. -DSGGK_SDK_DIR=<root>      explicit, highest priority (typical in CI)
#     2. env var SGGK_ROOT          machine-wide install location
#     3. <project>/../SGGK          native layout when this component sits inside the SGK repo
#     4. <project>/third_party/SGGK recommended drop-in location for a standalone checkout
#     5. <project>/SGGK             SDK embedded inside this project
#
# Outputs:
#     SGGK_FOUND / SGGK_VERSION
#     SGGK_SDK_DIR       SDK root
#     SGGK_INCLUDE_DIR   <root>/include
#     SGGK_LIB_DIR       <root>/x64-win/lib   (Release import libraries)
#     SGGK_BIN_DIR       <root>/x64-win/bin   (Release runtime DLLs)
#     SGGK_LICENSE_FILE  <root>/sggk.lic      (when present)
#
# On success the resolved root is written back into the cache variable SGGK_SDK_DIR, so a stale
# cached value is corrected automatically instead of being silently reused.
# =============================================================================================

include_guard(GLOBAL)

set(SGGK_SDK_DIR "" CACHE PATH "SGGK 1.7 SDK root (must contain include/ and x64-win/)")

# Build the candidate list: explicit override, environment, then project-relative layouts.
set(_sggk_candidates "")

if(SGGK_SDK_DIR)
    list(APPEND _sggk_candidates "${SGGK_SDK_DIR}")
endif()

if(DEFINED ENV{SGGK_ROOT} AND NOT "$ENV{SGGK_ROOT}" STREQUAL "")
    list(APPEND _sggk_candidates "$ENV{SGGK_ROOT}")
endif()

list(APPEND _sggk_candidates
    "${CMAKE_CURRENT_SOURCE_DIR}/../SGGK"
    "${CMAKE_CURRENT_SOURCE_DIR}/third_party/SGGK"
    "${CMAKE_CURRENT_SOURCE_DIR}/SGGK"
)

# Validate candidates one by one; first qualifying directory wins.
set(SGGK_FOUND FALSE)

foreach(_cand IN LISTS _sggk_candidates)
    if(_cand STREQUAL "")
        continue()
    endif()

    get_filename_component(_root "${_cand}" ABSOLUTE)

    if(EXISTS "${_root}/include/Foundation/SGGK_Version.h"
       AND EXISTS "${_root}/x64-win/lib/SGGK_Foundation.lib")

        set(SGGK_SDK_DIR     "${_root}"                 CACHE PATH "SGGK 1.7 SDK root (must contain include/ and x64-win/)" FORCE)
        set(SGGK_INCLUDE_DIR "${_root}/include")
        set(SGGK_LIB_DIR     "${_root}/x64-win/lib")
        set(SGGK_BIN_DIR     "${_root}/x64-win/bin")

        if(EXISTS "${_root}/sggk.lic")
            set(SGGK_LICENSE_FILE "${_root}/sggk.lic")
        endif()

        # Version: extract SGGK_VERSION_{MAJOR,MINOR,PATCH} from SGGK_Version.h
        foreach(_part MAJOR MINOR PATCH)
            file(STRINGS "${_root}/include/Foundation/SGGK_Version.h" _ver_line
                 REGEX "^#define[ \t]+SGGK_VERSION_${_part}[ \t]+[0-9]+")
            if(_ver_line MATCHES "SGGK_VERSION_${_part}[ \t]+([0-9]+)")
                set(_sggk_ver_${_part} "${CMAKE_MATCH_1}")
            endif()
        endforeach()
        if(DEFINED _sggk_ver_MAJOR)
            set(SGGK_VERSION "${_sggk_ver_MAJOR}.${_sggk_ver_MINOR}.${_sggk_ver_PATCH}")
        endif()

        set(SGGK_FOUND TRUE)
        break()
    endif()
endforeach()

# Report the result; on failure list every candidate that was tried.
if(NOT SGGK_FOUND)
    message(FATAL_ERROR
        "\n[FindSGGK] SGGK 1.7 SDK not found.\n"
        "\nProvide it in one of these ways (in priority order):\n"
        "  1) explicit:   cmake -S . -B build -DSGGK_SDK_DIR=<SDK root>\n"
        "  2) env var:    set SGGK_ROOT=<SDK root>        (Linux/macOS: export SGGK_ROOT=...)\n"
        "  3) native:     place this component inside the SGK repo, next to the vendor SGGK/ dir\n"
        "  4) preferred:  copy the SDK to <this project>/third_party/SGGK/\n"
        "\nThe SDK root must contain BOTH:\n"
        "    include/Foundation/SGGK_Version.h\n"
        "    x64-win/lib/SGGK_Foundation.lib\n"
        "    (optional) sggk.lic\n"
        "\nCandidate directories tried:\n"
        "    ${_sggk_candidates}\n")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SGGK
    REQUIRED_VARS SGGK_SDK_DIR SGGK_INCLUDE_DIR SGGK_LIB_DIR
    VERSION_VAR   SGGK_VERSION)
