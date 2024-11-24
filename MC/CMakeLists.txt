cmake_minimum_required(VERSION 3.16)
#we gave the name for the project
project(mc_sim)
#cpp standard is needed to avoid dtupid errors
set(CMAKE_CXX_STANDARD 14)

#find root
find_package(ROOT REQUIRED COMPONENTS RIO Net Hist Gpad Tree)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

#setting variables with info about our library
set(SOURCE_LIB generator.cpp)
add_library(generator SHARED ${SOURCE_LIB}) 
target_link_libraries(generator ROOT::RIO ROOT::Net ROOT::Hist ROOT::Gpad ROOT::Tree)

#main line where we say tat we need an exe file which to be created from .C
add_executable(main.exe montecarlosim_exp_b.C)

target_link_libraries(main.exe generator ROOT::RIO)#ROOT::RIO ROOT::Net ROOT::Hist ROOT::Gpad ROOT::Tree)

