add_executable(skywalker-cpu skywalker.cpp)
add_dependencies(apps skywalker-cpu)
target_link_libraries(skywalker-cpu PRIVATE Galois::shmem lonestar)

target_include_directories(skywalker-cpu PRIVATE "/proj/prismgt-PG0/vrao79/metis-install/include")

find_library(GKLIB GKlib "/proj/prismgt-PG0/vrao79/metis-install/lib")

if(NOT GKLIB)
  message(FATAL_ERROR "GKLib not found")
endif()

find_library(METIS metis "/proj/prismgt-PG0/vrao79/metis-install/lib")

if(NOT METIS)
  message(FATAL_ERROR "METIS not found")
endif()

target_link_libraries(skywalker-cpu PRIVATE "${GKLIB}")
target_link_libraries(skywalker-cpu PRIVATE "${METIS}")

install(TARGETS skywalker-cpu DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT apps EXCLUDE_FROM_ALL)

add_test_scale(small1 skywalker-cpu "${BASEINPUT}/reference/structured/rome99.gr" -delta 8)
add_test_scale(small2 skywalker-cpu "${BASEINPUT}/scalefree/rmat10.gr" -delta 8)
