#include "include.h"




//*
const bool performance_test_gpu=false;
const bool test_gpu=false;
const bool solve_linear_congruence_test_gcd=false;
const bool use_gpu_advance=true;
const bool gpu_is_async=true;
#define USE_GPU

//cuda streams have a buffer limit of 1024 pending operations after which it will block
//this isn't used because compute preemption seems to work good enough
const int gpu_iterations_per_call=1000000000; //10 = 5ms

//half_d_num_bits has the same number of bits as the order; increasing it by 2 will double the runtime
//can do about 26 million iterations per second
//cpu code is 50k iterations per second on a single sandy bridge core at 3ghz
//expected runtime is sqrt(pi*n/2) iterations where n is the order
//if increasing the exponent by 1 randomizes the quadratic form, then the order will be about the square root of the classgroup order
// because of the birthday paradox. the classgroup order is bounded by the discriminant size because a and b are half the discriminant
// size and c is uniquely determined by a and b
// n=2^half_d_num_bits ; iterations=sqrt(pi/2)*2^(half_d_num_bits/2)*0.75
//
//this has a birthday problem distribution; the above times are for a 50% probability
//if the time is doubled, the probability is 97%
//if the time is tripled, the probability is 99.97%
//the actual times are about 75% of birthday paradox times; this was observed for other pollard rho implementations that use large
// numbers of subsets. could also be because the order is an upper bound on the actual order
const int half_d_num_bits=86;

// 70 bits:
// -963475864336883842909521796527706961553559
// -785552272377590740002371572305022106240719
// -1242863421552049670288326799371773042627263

// 83 bits:
// -64657770747066387997612462548213558686444281012583
// -52717520611878693618478513523283894958126407695671
// -83407152327511170244624163966595602523359911573207

// 80 bits:
// -1010277667922912312462694727315836854475696256727
// -823711259560604587788726773801310858720727224239
// -1303236755117362035072252561978056289427500699839

// 85.5 bits:
// -2069048663906124415923598801542833877966216965923727
// -1686960659580118195791312432745084638660044959065487
// -2669028874480357447827973246931059280747517146469511

// 85 bits:
// -1034524331953062207961799400771416938983108482273447
// -843480329790059097895656216372542319330022479468847
// -1334514437240178723913986623465529640373758570232911

//leave empty to randomly generate a value
const string d_value = "-52717520611878693618478513523283894958126407695671";

//18 = 0.5MB file size increase per 1 minute; increasing by 1 will halve this
//     also, will generate a checkpoint once per 2 minutes; increasing by 1 will double this
//anything less than 18 will result in the cpu bottlenecking the gpu (since it has to verify each state)
//(could probably just verify the current state and assume the distinguished states are correct)
//
//the number of distinguished states per thread is a binomial distribution with n=num_checkpoints*checkpoint_interval and
// p=2^-distinguished_num_bits
//distance between two consecutive distinguished states is a geometric distribution
//if the distance is very long then the collision resolution might take a long time since it has to replay all of the states between
// the two distinguished states. if that happens it is probably better to just ignore the state

#ifndef THREADS_PER_BLOCK_FACTOR
    //have 1024 registers per block
    //this is 64 registers per thread; need to compile the kernel with that setting
    #define THREADS_PER_BLOCK_FACTOR 16
#endif

//22 seconds to validate a checkpoint with these settings. also need to validate each distinguished point
const int distinguished_num_bits=20; //22; //18
const int gpu_num_blocks=28*2;
const int gpu_threads_per_block=32*THREADS_PER_BLOCK_FACTOR;
const int save_checkpoints_interval=1;

const int checkpoint_interval=(performance_test_gpu)? 1<<18 : 1<<19;
//***/




/*
const bool performance_test_gpu=false;
const bool test_gpu=true;
const bool solve_linear_congruence_test_gcd=false;
const bool use_gpu_advance=false;
const bool gpu_is_async=true;

const int half_d_num_bits=40;
const int distinguished_num_bits=10;
const int num_threads=10;
const int gpu_num_blocks=1;
const int gpu_threads_per_block=10;
const int save_checkpoints_interval=0;
const int checkpoint_interval=1<<distinguished_num_bits;
const string d_value;
//***/




const int distinguished_buffer_size=32;
const int d_seed=0;
const int d_num_bits=2*half_d_num_bits;
const int exponent_num_bits=d_num_bits;
const int table_num_bits=exponent_num_bits+64; //need some extra bits
const int max_process_collision_iterations=(1<<(distinguished_num_bits)) * 20;
const int num_threads=gpu_num_blocks*gpu_threads_per_block;
//const double total_required_states_factor=0.75; //num required for birthday paradox is multiplied by this
const double total_required_states_factor=1;
const bool enable_track_cycles=false;

#ifndef __CUDA_ARCH__
    #define TEST_GPU_CODE
#endif

#ifdef USE_GPU
    #include <cuda_runtime.h>

    void error_check_impl(cudaError_t v, const char* file, int line) {
        if (v==cudaSuccess) {
            return;
        }

        print( "error_check failed. v=", v, "file=", file, "line=", line );

        const char* name=cudaGetErrorName(v);
        if (name!=nullptr) {
            print(name);
        }

        const char* description=cudaGetErrorString(v);
        if (description!=nullptr) {
            print(description);
        }

        assert(false);
    }
    #define error_check(v) error_check_impl((v), __FILE__, __LINE__)
#endif

#include "bit_manipulation.h"
#include "double_utility.h"

#include "integer.h"
#include "vdf_new.h"

#include "gpu_integer.h"
#include "gpu_integer_divide.h"
#include "gpu_integer_gcd.h"
#include "gpu_integer_reduce.h"
#include "gpu_integer_vdf.h"

#include "rho.h"
#include "gpu_integer_rho.h"

int main(int argc, char** argv) {
    set_rounding_mode();

    //todo //have to use one of the allowed discriminants for the final program; assign it here
    //integer d=generate_discriminant(d_num_bits, d_seed);

    integer d;
    if (d_value.empty()) {
        print( "Generating random value for d" );
        d=generate_discriminant(d_num_bits, d_seed);
    } else {
        d=integer(d_value);
    }

    double d_double=-mpz_get_d(d.impl);

    test_form(d);

    database c_database(d, exponent_num_bits, table_num_bits, distinguished_num_bits, checkpoint_interval);

    vector<rho_thread> threads;
    if (argc==2) {
        ifstream in(argv[1]);
        c_database.load(in);
        threads=c_database.restore_checkpoint();
    } else {
        for (int x=0;x<num_threads;++x) {
            threads.push_back(c_database.make_thread(x));

            if (performance_test_gpu && x>=1000) {
                threads.back()=threads[rand()%1000]; //faster
            } else {
                threads.back().init();
            }
        }

        if (!performance_test_gpu) {
            c_database.add_checkpoint(threads);
        }
    }

    unique_ptr<gpu_state> c_gpu_state=make_unique<gpu_state>();
    if (use_gpu_advance) {
        c_gpu_state->begin(threads, c_database.c_table);
    }

    time_t last_checkpoint_time=time(nullptr);

    int save_checkpoints_counter=0;
    while (true) {
        if (use_gpu_advance) {
            c_gpu_state->read_next(threads);
        } else {
            for (int y=0;y<checkpoint_interval;++y) {
                for (int x=0;x<num_threads;++x) {
                    threads.at(x).advance();
                }
            }
        }

        if (performance_test_gpu) {
            time_t c_checkpoint_time=time(nullptr);
            double delta_seconds=double(c_checkpoint_time)-double(last_checkpoint_time);
            double states_per_second=double(checkpoint_interval)*double(num_threads)/delta_seconds;

            print(states_per_second);
            return 0;
        }

        c_database.add_checkpoint(threads);

        //the progress of the algorithm is based on the number of distinguished states
        //if the number of distinguished states fails to increase, then distinguished_num_bits is too high and all of the threads have
        // entered into loops without any distinguished states
        //if that happens then it ought to be possible to lower distinguished_num_bits and restart the program from the last checkpoint
        // but this will involve replaying lots of states when a collision is found
        print( "Checkpoint. Num distinguished states:", c_database.distinguished_states.size() );
        {
            double num_states_processed=double(c_database.checkpoint_sequence_number)*double(num_threads);
            time_t c_checkpoint_time=time(nullptr);

            double delta_seconds=double(c_checkpoint_time)-double(last_checkpoint_time);

            const double pi=3.14159265358979323846;

            //50% likely to get a collision by this point
            double num_states_total=sqrt(pi/2)*pow(d_double, 1.0/4)*total_required_states_factor;

            double percent_done=num_states_processed/num_states_total*100;

            double states_per_second=double(checkpoint_interval)*double(num_threads)/delta_seconds;

            double hours_remaining=(num_states_total-num_states_processed)/states_per_second/60/60;

            print( "Percent done:", percent_done, "%" );
            print( "Hours remaining:", hours_remaining );
            print( "States per second:", states_per_second / 1e6, "million" );

            last_checkpoint_time=c_checkpoint_time;
        }

        if (c_database.process_collisions()) {
            break;
        }

        if (save_checkpoints_interval>=1 && save_checkpoints_counter%save_checkpoints_interval==0) {
            c_database.save();
            ++save_checkpoints_counter;
        }
    }

    //track_max.output(d_num_bits/2);
}