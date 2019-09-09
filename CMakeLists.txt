cmake_minimum_required(VERSION 3.9)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

project (growing_balls)
find_package(CGAL REQUIRED Core)

add_subdirectory(vendor/libdts2 libdts2)

set(MY_C_FLAGS "${CMAKE_C_FLAGS} -frounding-math -std=c99 -Wall -Wextra -Wno-unknown-pragmas -Wstrict-overflow=4")
set(MY_CXX_FLAGS "${CMAKE_CXX_FLAGS} -frounding-math -std=c++14 -Wall -Wextra -Wno-unknown-pragmas -Wstrict-overflow=4 -Wshadow")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MY_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MY_CXX_FLAGS}")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${DEBUG_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${DEBUG_FLAGS}")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} ${DEBUG_FLAGS}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${DEBUG_FLAGS}")

set(GROWING_BALLS_SOURCES
	src/growing_balls.cpp
	src/helpers/timer.cpp
	src/io/io.cpp
	src/primitives/pointofinterest.cpp
)
set(MY_INCLUDE_DIRS
	src/algorithms
	src/geofunctions
	src/helpers
	src/io
	src/primitives
	src/storage
)

add_executable(${PROJECT_NAME} ${GROWING_BALLS_SOURCES})
target_link_libraries(${PROJECT_NAME} dts2)
target_include_directories(${PROJECT_NAME} PUBLIC ${MY_INCLUDE_DIRS})
