cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

project (gbons2)
find_package(CGAL REQUIRED Core)
find_package(Boost REQUIRED)

add_subdirectory(vendor/libdts2 libdts2)

include_directories(
	src/
	src/algorithms/
	src/geofunctions/
	src/helpers/
	src/io/
	src/primitives/
	src/storage/
	${LIBDTS2_INCLUDE_DIR}
	)

set(MY_C_FLAGS "${CMAKE_C_FLAGS} -frounding-math -std=c99 -Wall -Wextra -Wno-unknown-pragmas -Wstrict-overflow=4")
set(MY_CXX_FLAGS "${CMAKE_CXX_FLAGS} -frounding-math -std=c++17 -Wall -Wextra -Wno-unknown-pragmas -Wstrict-overflow=4 -Wshadow")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MY_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MY_CXX_FLAGS}")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${DEBUG_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${DEBUG_FLAGS}")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} ${DEBUG_FLAGS}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${DEBUG_FLAGS}")

SET(MY_LIB_LINK_LIBS
	${LIBRT_LIBRARIES}
	${LIBDTS2_LIBRARIES}
	${LIBDTS2_LINK_LIBRARIES}
	${CGAL_LIBRARIES}
	${CGAL_3RD_PARTY_LIBRARIES}
	)

#set(GBS_ON_S2_SOURCES_CPP
#	src/gbons2.cpp
#)
#
#add_executable(${PROJECT_NAME} ${GBS_ON_S2_SOURCES_CPP})

set(GROWING_BALLS_SOURCES
	src/growing_balls.cpp
	src/helpers/timer.cpp
	src/io/io.cpp
	src/primitives/pointofinterest.cpp
	)

add_executable(growing_balls ${GROWING_BALLS_SOURCES})

#target_link_libraries(${PROJECT_NAME} ${MY_LIB_LINK_LIBS})
#add_dependencies(${PROJECT_NAME} dts2)
#target_compile_definitions(${PROJECT_NAME} PUBLIC ${LIBDTS2_COMPILE_DEFINITIONS})

target_link_libraries(growing_balls ${MY_LIB_LINK_LIBS})
add_dependencies(growing_balls dts2)
target_compile_definitions(growing_balls PUBLIC ${LIBDTS2_COMPILE_DEFINITIONS})
