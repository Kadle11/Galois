add_executable(skywalker-cpu skywalker.cpp)
add_dependencies(apps skywalker-cpu)
target_link_libraries(skywalker-cpu PRIVATE Galois::shmem lonestar)
install(TARGETS skywalker-cpu DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT apps EXCLUDE_FROM_ALL)

add_test_scale(small1 skywalker-cpu "${BASEINPUT}/reference/structured/rome99.gr" -delta 8)
add_test_scale(small2 skywalker-cpu "${BASEINPUT}/scalefree/rmat10.gr" -delta 8)
