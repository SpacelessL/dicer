set(SOURCES
	dicer.h
	dicer/dice.h
	dicer/symbol_set.h
	dicer/monomial.h
	dicer/polynomial.h
	dicer/fft.h
	dicer/fft.cpp
	dicer/hash.h
	dicer/utils.h
)

add_subdirectory(sus)

add_library(dicer ${SOURCES})
set_target_properties(dicer PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(dicer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(dicer PUBLIC "/Zc:preprocessor" "/openmp")

target_link_libraries(dicer PUBLIC sus)

if (NOT DEFINED USE_ONEMKL)
	set(USE_ONEMKL true)
endif()

if (NOT DEFINED USE_MIMALLOC)
	set(USE_MIMALLOC true)
endif()

if (USE_ONEMKL)
	add_definitions(-DUSE_ONEMKL)

	find_package(MKL CONFIG REQUIRED)
	target_compile_options(dicer PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
	target_include_directories(dicer PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
	target_link_libraries(dicer PUBLIC $<LINK_ONLY:MKL::MKL>)
endif()

if (USE_MIMALLOC)
	add_definitions(-DUSE_MIMALLOC)

	set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}" "C:\\Program Files (x86)\\libmimalloc")
	find_package(mimalloc 2.0 REQUIRED)
	message(STATUS "Found mimalloc installed at: ${MIMALLOC_LIBRARY_DIR} (${MIMALLOC_VERSION_DIR})")

	target_link_libraries(dicer PUBLIC mimalloc)
endif()
