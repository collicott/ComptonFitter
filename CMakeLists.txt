project(ComptonFitter)
cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# look for APLCON fitter package
# don't forget to build the library  
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/" ${CMAKE_MODULE_PATH})
find_package(APLCON REQUIRED)
find_package(ROOT "5.30" REQUIRED)

include_directories(${APLCON_INCLUDE_DIR})

# since all subprojects need ROOT, set that up here
include_directories(${ROOT_INCLUDES})
link_directories(${ROOT_LIBRARY_DIR})

include_directories(inc)
set(ROOT_LIBRARIES "${ROOT_LIBRARIES} -lEG")

set(SRC_LIST src/CEFTwrapper.cpp src/DataHandling.cpp)

add_executable(phil_fit ${SRC_LIST} src/phil_fit.cpp)

    target_link_libraries(phil_fit
    	${APLCON_LIBRARIES}
    	${ROOT_LIBRARIES}
    )

add_executable(SP_fit ${SRC_LIST} src/SP_fit.cpp)

    target_link_libraries(SP_fit
    	${APLCON_LIBRARIES}
    	${ROOT_LIBRARIES}
    )

