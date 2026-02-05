# Patch Catch2 for C++20: std::uncaught_exception() was removed.
# Run from FetchContent PATCH_COMMAND with WORKING_DIRECTORY = catch2 source.
# Replaces "return std::uncaught_exception();" with "return std::uncaught_exceptions() > 0;"
set(FILE "src/catch2/internal/catch_uncaught_exceptions.cpp")
if(NOT EXISTS "${FILE}")
  message(FATAL_ERROR "patch_catch2_uncaught: ${FILE} not found (cwd=${CMAKE_CURRENT_LIST_DIR})")
endif()
file(READ "${FILE}" content)
string(REPLACE "return std::uncaught_exception();" "return std::uncaught_exceptions() > 0;" content "${content}")
file(WRITE "${FILE}" "${content}")
message(STATUS "Catch2 C++20 patch applied to ${FILE}")
