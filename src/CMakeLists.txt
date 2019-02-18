cmake_minimum_required (VERSION 2.6)
project (SVDQuest)

file(GLOB SOURCES *.cpp)

if(WIN32)
add_definitions(/std:c++14)
INCLUDE_DIRECTORIES("C:/Users/pranjal/src/phylokit/install/x64-Release/include")
INCLUDE_DIRECTORIES("C:/Users/pranjal/src/phylonaut/install/x64-Release/include")
INCLUDE_DIRECTORIES("boost/")


if(WIN32)
LINK_DIRECTORIES("C:/Users/pranjal/src/phylokit/install/x64-Release/lib/")

LINK_DIRECTORIES("C:/Users/pranjal/src/phylonaut/install/x64-Release/lib/")
endif()

else()
add_definitions(-std=c++14 -g -O3 -Wall)
endif()

add_executable(SVDquest ${SOURCES} SVDQuest.cpp)
target_link_libraries(SVDquest phylokit phylonaut)


# SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
# SET(BUILD_SHARED_LIBRARIES OFF)
# SET(CMAKE_EXE_LINKER_FLAGS "-static")
# add_executable(SVDQuest.dist ${SOURCES} SVDQuest.cpp)
# target_link_libraries(SVDQuest.dist phylokit phylonaut)

install(TARGETS SVDquest 
		RUNTIME DESTINATION bin)
