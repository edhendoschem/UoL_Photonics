/*--------------Vector definitions, #2--------------*/
#include "definitions.h"
Flat_vec::Flat_vec(unsigned long long lengthh)
: length{lengthh} {

    vec_size = length;
    width = 0;
    depth = 0;

    //Initialize the vector
    for (unsigned long long i = 0; i < length; ++i) {
        vec.push_back(0.0);
    }
}

Flat_vec::Flat_vec(unsigned long long lengthh, unsigned long long widthh)
: length{lengthh}, width{widthh} {

    vec_size = length * width;
    depth = 0;

    //Initialize the vector
    for (unsigned long long i = 0; i < length; ++i) {
        for (unsigned long long j = 0; j < width; ++j) {
            vec.push_back(0.0);
        }
    }
}

Flat_vec::Flat_vec(unsigned long long lengthh, unsigned long long widthh,unsigned long long depthh)
: length{lengthh}, width{widthh}, depth{depthh} {

    vec_size = length * width * depth;

    //Initialize the vector
    for (unsigned long long i = 0; i < length; ++i) {
        for (unsigned long long j = 0; j < width; ++j) {
            for (unsigned long long k = 0; k < depth; ++k) {
                vec.push_back(0.0);
            }
        }
    }
}

unsigned long long Flat_vec::get_ind_3D(unsigned long long i, unsigned long long j, unsigned long long k) {
    unsigned long long s = k + depth * (j + width * i);
    if (s >= vec_size) {
        cout<<"Warning invalid index, returning 0"<<'\n';
        return 0;
    }else {
        return s;
    }
}

unsigned long long Flat_vec::get_ind_2D(unsigned long long i, unsigned long long j) {
    unsigned long long s = j + width * i;
    if (s >= vec_size) {
        cout<<"Warning invalid index, returning 0"<<'\n';
        return 0;
    }else {
        return s;
    }
}
