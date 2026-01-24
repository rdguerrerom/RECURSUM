# Build single combined Python extension module
set(BINDING_FILE "${CMAKE_CURRENT_SOURCE_DIR}/bindings/recursum_bindings.cpp")

if(EXISTS ${BINDING_FILE})
    message(STATUS "Building combined RECURSUM module")

    pybind11_add_module(_recursum ${BINDING_FILE})

    target_include_directories(_recursum PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/src/generated
    )

    # Set output directory to recursum package
    set_target_properties(_recursum PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/recursum"
        OUTPUT_NAME "_recursum"
    )

    # Enable compile-time optimizations
    # Note: pybind11 requires exceptions and RTTI
    target_compile_options(_recursum PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang>:-O3 -ffast-math>
        $<$<CXX_COMPILER_ID:MSVC>:/O2 /EHsc>
    )

    # Link time optimization
    set_target_properties(_recursum PROPERTIES
        INTERPROCEDURAL_OPTIMIZATION TRUE
    )
else()
    message(FATAL_ERROR "Binding file not found: ${BINDING_FILE}")
endif()

# Print configuration summary
message(STATUS "RECURSUM Configuration:")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  Compiler: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  SIMD native arch: ${RECURSUM_USE_NATIVE_ARCH}")
