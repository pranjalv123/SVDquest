cmake_minimum_required (VERSION 2.6)
project (SVDQuest)

file(GLOB SOURCES *.cpp)

add_definitions(-std=c++14 -g -O3 -Wall)

add_executable(SVDquest ${SOURCES} SVDQuest.cpp)
target_link_libraries(SVDquest phylokit phylonaut)


# SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
# SET(BUILD_SHARED_LIBRARIES OFF)
# SET(CMAKE_EXE_LINKER_FLAGS "-static")
# add_executable(SVDQuest.dist ${SOURCES} SVDQuest.cpp)
# target_link_libraries(SVDQuest.dist phylokit phylonaut)

install(TARGETS SVDquest 
		RUNTIME DESTINATION bin)
