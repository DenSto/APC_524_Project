#include <stdio.h>
#include <assert.h>
#include <string.h>

#if USE_MPI
#include "mpi.h"
#endif

#include "input.hpp"
#include "../globals.hpp"

Input::Input(void){
    input_info_ = new Input_Info_t; 
}

Input::~Input(void){
    delete input_info_;
}

Input_Info_t* Input::getinfo(void){
    return input_info_;
} 

