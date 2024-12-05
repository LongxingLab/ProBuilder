project(gzip)
aux_source_directory(. DIR_GZIP_SRCS)
add_library(${PROJECT_NAME} STATIC ${DIR_GZIP_SRCS})
target_link_libraries(${PROJECT_NAME} z)
