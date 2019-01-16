typedef gpu_form<half_d_num_bits> gform;
const int table_num_words=bits_to_words(table_num_bits);

struct gpu_generator_table {
    gform values[table_num_bits];

    void assign(const generator_table& t) {
        assert(t.values.size()==table_num_bits);
        for (int x=0;x<table_num_bits;++x) {
            values[x]=t.values[x];
            assert(form(values[x])==t.values[x]);
        }
    }
};

struct gpu_exponent {
    gpu_integer<table_num_words, false> bits;
    gform value;

    void assign_to(exponent& e) const {
        e.bits=integer(bits);
        e.value=value;
    }

    void assign(const exponent& e) {
        bits=e.bits;
        value=e.value;
    }

    GFUNC void add_bit(int index, const gpu_generator_table* c_table) {
        gform f=c_table->values[index];
        gform new_value=multiply(value, f);

        gpu_integer<table_num_words, false> delta_bits;
        int index_copy=index;

        #pragma unroll
        for (int x=0;x<table_num_words;++x) {
            if (index_copy>=0 && index_copy<32) {
                delta_bits[x]|=1u << index_copy;
            }
            index_copy-=32;
        }

        auto new_bits=bits + delta_bits;

        #ifdef TEST_GPU_CODE
        {
            integer new_bits_int=integer(bits) + (integer(1)<<index);
            assert(new_bits_int==integer(new_bits));
        }
        #endif

        value=new_value;
        bits=new_bits;
    }
};

struct gpu_thread_state {
    gpu_exponent state;
    int thread_number=-1;
    int64 sequence_number=0;

    void assign_to(thread_state& s) const {
        state.assign_to(s.state);
        s.thread_number=thread_number;
        s.sequence_number=sequence_number;
    }

    void assign(const thread_state& s) {
        state.assign(s.state);
        thread_number=s.thread_number;
        sequence_number=s.sequence_number;
    }

    GFUNC bool is_distinguished() const {
        int new_hash=state.value.hash();
        return (new_hash & ((1<<distinguished_num_bits)-1)) == 0;
    }
};

struct gpu_rho_thread {
    gpu_thread_state state;
    gpu_thread_state distinguished_states[distinguished_buffer_size];
    int num_distinguished_states=0;

    void assign(const rho_thread& t) {
        state.assign(t.state);

        num_distinguished_states=t.distinguished_states.size();
        for (int x=0;x<num_distinguished_states;++x) {
            distinguished_states[x].assign(t.distinguished_states[x]);
        }
    }

    void assign_to(rho_thread& t) {
        state.assign_to(t.state);

        //gpu indicates overflow by inrementing num_distinguished_states but not doing the write if it exceeds the buffer size
        assert(num_distinguished_states>=0 && num_distinguished_states<distinguished_buffer_size);

        t.distinguished_states.resize(num_distinguished_states, t.state);
        for (int x=0;x<num_distinguished_states;++x) {
            distinguished_states[x].assign_to(t.distinguished_states[x]);
        }
    }

    GFUNC void advance(const gpu_generator_table* c_table) {
        int c_hash=state.state.value.hash();
        int add_index=c_hash%exponent_num_bits;
        state.state.add_bit(add_index, c_table);
        ++state.sequence_number;

        if (state.is_distinguished()) {
            //could theoretically just discard the extra states if there are too many of them and everything will still work
            //anyway, if this happens then distinguished_buffer_size should be increased and the program restarted from the last
            // checkpoint
            //assert(num_distinguished_states>=0 && num_distinguished_states<distinguished_buffer_size);
            //not using asserts on the gpu because they have some kind of issues; this will get detected on the cpu if it fails
            if (num_distinguished_states>=0 && num_distinguished_states<distinguished_buffer_size) {
                distinguished_states[num_distinguished_states]=state;
            }
            ++num_distinguished_states;
        }
    }
};

struct gpu_rho_thread_array {
    gpu_rho_thread state[num_threads];

    void assign(const vector<rho_thread>& t) {
        assert(t.size()==num_threads);
        for (int x=0;x<num_threads;++x) {
            state[x].assign(t[x]);
        }
    }

    void assign_to(vector<rho_thread>& t) {
        assert(t.size()==num_threads);
        for (int x=0;x<num_threads;++x) {
            state[x].assign_to(t[x]);
        }
    }
};

int gpu_simulated_thread_number=-1;

#ifdef USE_GPU
    __global__
#endif
void gpu_advance(
    const gpu_rho_thread_array* in_data, gpu_rho_thread_array* out_data, const gpu_generator_table* c_table, int iterations
) {
    #ifdef USE_GPU
        int thread_number=blockIdx.x*blockDim.x + threadIdx.x;
    #else
        int thread_number=gpu_simulated_thread_number;
    #endif

    const gpu_rho_thread& in_state=in_data->state[thread_number];
    gpu_rho_thread& out_state=out_data->state[thread_number];

    if (&in_state!=&out_state) {
        out_state.state=in_state.state;
        out_state.num_distinguished_states=0;
    }

    #pragma unroll 1
    for (int x=0;x<iterations;++x) {
        out_state.advance(c_table);

        #ifdef USE_GPU
            //#pragma unroll 1
            //for (int x=0;x<100000000;++x) { out_state.state.sequence_number*=(out_state.state.sequence_number+1); }

            //prevent warps from becoming desynchronized and causing icache misses
            //scheduler ought to keep them synchronized for one iteration
            __syncthreads();
        #endif
    }
}

struct gpu_state {
    #ifdef USE_GPU
        cudaStream_t c_stream;
    #endif

    gpu_generator_table c_table;

    static const int num_buffers=2;
    gpu_rho_thread_array buffers[num_buffers];

    int64 read_sequence_number=0;
    int64 write_sequence_number=0;
    volatile int64 write_finished_sequence_number=0; //this gets written to from another thread

    gpu_rho_thread_array* device_copy[num_buffers];
    gpu_generator_table* device_copy_c_table;

    #ifdef USE_GPU
        static void CUDART_CB stream_finished_callback(cudaStream_t event, cudaError_t status, void *data) {
            gpu_state& this_v=*(gpu_state*)data;
            ++this_v.write_finished_sequence_number;
        }
    #endif

    bool run() {
        assert(gpu_num_blocks*gpu_threads_per_block==num_threads);

        assert(read_sequence_number<=write_sequence_number);

        int64 pending_read_buffers=write_sequence_number-read_sequence_number;
        assert(pending_read_buffers<=num_buffers);

        if (pending_read_buffers==num_buffers) {
            return false;
        }

        int in_buffer=write_sequence_number%num_buffers;
        int out_buffer=(write_sequence_number+1)%num_buffers;

        for (int x=0;x<checkpoint_interval;x+=gpu_iterations_per_call) {
            int num_iterations=checkpoint_interval-x;
            if (num_iterations>gpu_iterations_per_call) {
                num_iterations=gpu_iterations_per_call;
            }

            int c_in_buffer=(x==0)? in_buffer : out_buffer;

            #ifdef USE_GPU
                gpu_advance<<<gpu_num_blocks, gpu_threads_per_block, 0, c_stream>>>(
                    device_copy[c_in_buffer],
                    device_copy[out_buffer],
                    device_copy_c_table,
                    num_iterations
                );
                error_check(cudaGetLastError());
            #else
                for (int y=0;y<num_threads;++y) {
                    gpu_simulated_thread_number=y;
                    gpu_advance(&buffers[c_in_buffer], &buffers[out_buffer], &c_table, num_iterations);
                }
            #endif
        }

        #ifdef USE_GPU
            error_check(cudaStreamAddCallback(c_stream, stream_finished_callback, this, 0));

            if (!gpu_is_async) {
                error_check(cudaStreamSynchronize(c_stream));
                assert(write_finished_sequence_number==write_sequence_number+1);
            }
        #else
            ++write_finished_sequence_number;
        #endif

        ++write_sequence_number;

        return true;
    }

    void begin(const vector<rho_thread>& t, const generator_table& t_table) {
        c_table.assign(t_table);

        buffers[0].assign(t);

        #ifdef USE_GPU
            for (int x=0;x<num_buffers;++x) {
                error_check(cudaMalloc((void**)&(device_copy[x]), sizeof(buffers[x])));
                error_check(cudaMemcpy(device_copy[x], &(buffers[x]), sizeof(buffers[x]), cudaMemcpyHostToDevice));
            }

            error_check(cudaMalloc((void**)&(device_copy_c_table), sizeof(c_table)));
            error_check(cudaMemcpy(device_copy_c_table, &c_table, sizeof(c_table), cudaMemcpyHostToDevice));

            error_check(cudaStreamCreateWithFlags(&c_stream, cudaStreamNonBlocking));
        #endif

        read_sequence_number=0;
        write_sequence_number=0;

        while (run());
    }

    void read_next(vector<rho_thread>& t) {
        assert(read_sequence_number<write_sequence_number);

        #ifdef USE_GPU
            while (read_sequence_number>=write_finished_sequence_number) {
                usleep(10000); //good enough
            }
        #endif

        assert(read_sequence_number<write_finished_sequence_number);

        int in_buffer=(read_sequence_number+1)%num_buffers;

        #ifdef USE_GPU
            error_check(cudaMemcpy(
                &(buffers[in_buffer]), device_copy[in_buffer], sizeof(buffers[in_buffer]), cudaMemcpyDeviceToHost
            ));
        #endif
        buffers[in_buffer].assign_to(t);

        ++read_sequence_number;

        while (run());
    }
};
