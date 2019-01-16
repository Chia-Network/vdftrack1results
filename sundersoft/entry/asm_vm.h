namespace simd_integer_namespace {


string vpermq_mask(array<int, 4> lanes) {
    int res=0;
    for (int x=0;x<4;++x) {
        int lane=lanes[x];
        assert(lane>=0 && lane<4);
        res|=lane << (2*x);
    }
    return to_hex(res);
}

string vpblendd_mask_4(array<int, 4> lanes) {
    int res=0;
    for (int x=0;x<4;++x) {
        int lane=lanes[x];
        assert(lane>=0 && lane<2);
        res|=((lane==1)? 3 : 0) << (2*x);
    }
    return to_hex(res);
}

string vpblendd_mask_8(array<int, 8> lanes) {
    int res=0;
    for (int x=0;x<8;++x) {
        int lane=lanes[x];
        assert(lane>=0 && lane<2);
        res|=((lane==1)? 1 : 0) << x;
    }
    return to_hex(res);
}

//layout: [immediates] [spilled registers] |0| [divide table] [integers]
//start of memory should be 4k aligned
struct asm_memory_type {
    const int initial_pos=(num_immediates+num_spilled_registers)*32;

    reg_scalar memory_base; //global register; 32-aligned
    int allocated_size=initial_pos;
    map<int, uint64> initial_data;

    vector<char> memory_buffer;
    char* memory_buffer_base=nullptr;

    void init_buffer() {
        assert(memory_buffer_base==nullptr && memory_buffer.empty());
        assert(initial_pos%4096==0);

        memory_buffer.resize(allocated_size+4096, 0);
        
        char* c_pos=&memory_buffer[0];
        while (uint64(c_pos)%4096!=0) {
            ++c_pos;
        }
        assert((c_pos-&memory_buffer[0])+allocated_size<=memory_buffer.size());

        memory_buffer_base=c_pos+initial_pos;
    }

    void assign_initial_value(int pos, int size, const char* initial_value) {
        pos+=initial_pos;

        assert(size%8==0 && size>0);
        assert(pos%8==0 && pos>=0);

        for (int x=0;x<size;x+=8) {
            assert(pos+x+8<=allocated_size);
            initial_data[pos+x]=*(uint64*)(initial_value+x);
        }
    }

    int alloc(int size, int alignment, const char* initial_value=nullptr) {
        assert(size%8==0);
        assert(alignment==1 || alignment==8 || alignment==32 || alignment==4096);

        int res=allocated_size-initial_pos;
        while (res%alignment!=0) {
            ++res;
        }

        allocated_size=initial_pos+res+size;

        if (initial_value!=nullptr) {
            assign_initial_value(res, size, initial_value);
        }

        return res;
    }

    string operator()(reg_scalar c_memory_base, int pos) {
        assert(initial_pos+pos>=0 && initial_pos+pos<allocated_size); //assumes c_memory_base==memory_base
        return str( "[#+#]", get_scalar_register_name(c_memory_base.value, 64), to_hex(pos) );
    }

    void assign_initial_data(reg_alloc regs) {
        EXPAND_MACROS_SCOPE;

        reg_scalar tmp=regs.bind_scalar(m, "tmp");

        for (auto c : initial_data) {
            if (c.second==0) {
                continue;
            }

            APPEND_M(str( "MOV `tmp, #", to_hex(c.second) ));
            APPEND_M(str( "MOV #, `tmp", (*this)(memory_base, c.first-initial_pos) ));
        }
    }
};
asm_memory_type asm_memory;

struct asm_immediate_type {
    map<array<uint64, 4>, int> values;
    int next_immediate=0;

    int lookup(array<uint64, 4> value) {
        auto i=values.find(value);
        if (i!=values.end()) {
            return i->second;
        }

        assert(next_immediate<num_immediates);

        int res=-(next_immediate+num_spilled_registers+1)*32;
        ++next_immediate;

        asm_memory.assign_initial_value(res, 32, (char*)&value[0]);

        values[value]=res;
        return res;
    }

    //this can change the flags register
    void assign(reg_scalar reg, uint64 value) {
        EXPAND_MACROS_SCOPE;
        m.bind(reg, "reg");

        if (value==0) {
            APPEND_M(str( "XOR `reg, `reg" )); //changes the flags register
        } else
        if (uint32(value)==uint64(value)) {
            APPEND_M(str( "MOV `reg_32, #", to_hex(value) ));
        } else
        if (int32(value)==int64(value)) {
            APPEND_M(str( "MOV `reg, #", to_hex(int64(value)) ));
        } else {
            int pos=lookup({value, value, value, value});
            APPEND_M(str( "MOV `reg, #", asm_memory(asm_memory.memory_base, pos) ));
        }
    }

    void assign(reg_vector reg, array<uint64, 4> value) {
        EXPAND_MACROS_SCOPE;
        m.bind(reg, "reg");

        if (value==array<uint64, 4>({0, 0, 0, 0})) {
            APPEND_M(str( "VPXOR `reg, `reg, `reg" ));
        } else
        if (value==array<uint64, 4>({~uint64(0), ~uint64(0), ~uint64(0), ~uint64(0)})) {
            APPEND_M(str( "VPCMPEQB `reg, `reg, `reg" ));
        } else {
            int pos=lookup(value);
            APPEND_M(str( "VMOVDQA `reg, #", asm_memory(asm_memory.memory_base, pos) ));
        }
    }

    void assign(reg_vector reg, uint64 value) {
        assign(reg, {value, value, value, value});
    }

    string operator()(array<uint64, 4> value) {
        return asm_memory(asm_memory.memory_base, lookup(value));
    }

    //address pointing to 4 copies of the value
    string operator()(uint64 value) {
        return (*this)({value, value, value, value});
    }
};
asm_immediate_type asm_immediate;

struct asm_function {
    string name;

    asm_function(string t_name) {
        EXPAND_MACROS_SCOPE;

        static bool outputted_header=false;
        if (!outputted_header) {
            APPEND_M(str( ".intel_syntax noprefix" ));
            outputted_header=true;
        }

        //first argument is in RDI
        assert(asm_memory.memory_base.value==reg_rdi.value);

        name=t_name;

        APPEND_M(str( ".global asm_func_#", t_name ));
        APPEND_M(str( "asm_func_#:", t_name ));

        APPEND_M(str( "PUSH RBX" ));
        APPEND_M(str( "PUSH RBP" ));
        APPEND_M(str( "PUSH R12" ));
        APPEND_M(str( "PUSH R13" ));
        APPEND_M(str( "PUSH R14" ));
        APPEND_M(str( "PUSH R15" ));
    }

    ~asm_function() {
        EXPAND_MACROS_SCOPE;

        asm_immediate.assign(reg_rax, 0);

        string end_label=m.alloc_label();
        APPEND_M(str( "#:", end_label ));
        APPEND_M(str( "POP R15" ));
        APPEND_M(str( "POP R14" ));
        APPEND_M(str( "POP R13" ));
        APPEND_M(str( "POP R12" ));
        APPEND_M(str( "POP RBP" ));
        APPEND_M(str( "POP RBX" ));
        APPEND_M(str( "RET" ));

        while (m.next_output_error_label_id<m.next_error_label_id) {
            APPEND_M(str( "label_error_#:", m.next_output_error_label_id ));

            assert(m.next_output_error_label_id!=0);
            asm_immediate.assign(reg_rax, m.next_output_error_label_id);
            APPEND_M(str( "JMP #", end_label ));

            ++m.next_output_error_label_id;
        }
    }
};


}