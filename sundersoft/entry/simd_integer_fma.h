namespace simd_integer_namespace {


//inputs are not sign-extended; assumed that the vector is nonnegative
template<int size> void matrix_vector_multiply(
    array<int64, size*size> c_matrix,
    array<simd_integer*, size> out_vector, array<simd_integer*, size> mul_vector, array<simd_integer*, size> add_vector,
    int integer_start, int integer_end, int shift, int truncate=-1
) {
    assert(size>=2);

    int integer_size=mul_vector[0]->current_size();

    for (int x=0;x<size;++x) {
        assert(add_vector[x]==nullptr || add_vector[x]->current_size()==integer_size);
        assert(mul_vector[x]->current_size()==integer_size);
        assert(out_vector[x]->current_size()==integer_size);
    }

    for (int x=integer_end-1;x>=integer_start;--x) {
        vector<uint64> vs;
        for (int y=0;y<size;++y) {
            int pos=x-shift;
            uint64 v=(pos<0)? 0 : mul_vector[y]->memory.at(pos);

            if (truncate!=-1 && pos>=truncate) {
                v=0;
            }

            vs.push_back(v);
        }

        vector<uint64> outs;
        for (int y=0;y<size;++y) {
            uint64 out=(add_vector[y]==nullptr)? 0 : add_vector[y]->memory.at(x);
            for (int z=0;z<size;++z) {
                out+=int64(c_matrix[y*size+z])*int64(vs[z]);
            }
            outs.push_back(out);
        }

        for (int y=0;y<size;++y) {
            out_vector[y]->memory.at(x)=outs[y];
        }
    }
}

void integer_multiply(simd_integer& c, simd_integer& a, simd_integer& b) {
    c.clear();
    int current_num_adds=0;

    //don't care what order this is in since only a single carry pass is done
    for (int x=0;x<b.current_size();++x) {
        int num=1; //min(b.current_size()-x, 8);

        vector<uint64> bs;
        for (int y=0;y<num;++y) {
            bs.push_back(b.memory.at(x+y));
        }

        if (current_num_adds+num>max_num_adds) {
            c.calculate_carry(0, true);
            current_num_adds=0;
        }

        //fma_asm(c, c, a, bs[0], x);
        c.fma(&c, a, bs, x, true, x+num==b.current_size());
        current_num_adds+=num;
    }

    c.calculate_carry();
}


}