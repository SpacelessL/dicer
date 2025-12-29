add_subdirectory(sus)

find_package(TBB CONFIG REQUIRED COMPONENTS tbb tbbmalloc)
find_package(unordered_dense CONFIG REQUIRED)

set(TBB_tbb_FOUND ${TBB_FOUND})
set(MKL_LINK static)
set(MKL_THREADING tbb_thread)
find_package(MKL CONFIG REQUIRED)

add_library(dicer
	dicer.h
	dicer/dice.h
	dicer/symbol_set.h
	dicer/monomial.h
	dicer/polynomial.h
	dicer/fft.h
	dicer/fft.cpp
	dicer/hash.hpp
	dicer/utils.hpp
)

set_target_properties(dicer PROPERTIES LINKER_LANGUAGE CXX)

if(MSVC)
	target_compile_options(dicer PUBLIC /Zc:preprocessor)
endif()

target_compile_definitions(dicer PUBLIC DICER_USE_ONEMKL)
target_include_directories(dicer PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}
	$<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>
)
target_compile_options(dicer PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
target_link_libraries(dicer PUBLIC sus TBB::tbb MKL::MKL unordered_dense::unordered_dense)
