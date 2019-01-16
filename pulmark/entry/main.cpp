/**
Copyright (C) 2018 Markku Pulkkinen

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
**/

#define ENABLE_THREADING
//#define ENABLE_INPUT_ARG_SANITY_CHECK

#include <string>
#include <thread>
// Form using GNU MP library for big integer arithmetic
#include "VdfFormGmp.h"
using Form = VdfFormGmp;

using namespace std;

/**
 * main - program entry point
 * Example call:
 * a.out 0x1234567890ABCDEF 1000000 2 1 4096
 * where
 *  discriminator = 0x1234567890ABCDEF
 *  iteration count = 1000000
 *  a = 2
 *  b = 1
 *  4096 = size of memory block (in bits) for GMP API mpz_t variables
 */
int main(int argc, char *argv[]) {
  /////////////////////////////////////////////////////////////////////////////
  // Note: The input argument sanity checks are disabled for competition.
  // The reason for this is that according to the judge the only criteria 
  // for competition entry is correct result and time. 
  // The developer decided to disable sanity checks for input args because 
  // in the competition rules there was no exact specification what the 
  // program output should be if input arguments are missing or are not valid. 
  // Disabling input argument sanity checks saves some time.
  /////////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_INPUT_ARG_SANITY_CHECK  
  // parse input args: discriminator and iter count > 0 are mandatory
  assert(argc >= 3);
#endif
  auto d(argv[1]);
  auto cnt(stoull(argv[2]));
#ifdef ENABLE_INPUT_ARG_SANITY_CHECK
  assert(cnt > 0);
  assert(d != nullptr);
#endif

  // init optional args to defaults
  unsigned long long a(2), b(1);
  unsigned long numBits(4096);
  if (argc >= 5) {
    // init a, b values from command line
    a = stoull(argv[3]);
    b = stoull(argv[4]);
  }
  if (argc >= 6)
    // space for mpz_t variables
    numBits = stoul(argv[5]);

  // create form, run squaring, display results
  Form y(a, b, d, static_cast<uint16_t>(numBits));
#ifdef ENABLE_THREADING
  std::thread t(&Form::square, &y, cnt);
  t.join();
#else
  y.square(cnt);
#endif
  gmp_printf("%Zd\n%Zd", y.A(), y.B());
  return 0;
}
