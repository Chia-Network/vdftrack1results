
install.sh
------------------------------------------------------

sudo apt-get install -y cmake &&                 # We use cmake to build the application. This utility generates make files.
sudo apt-get install -y m4 &&                    # m4 is required to build gmp from source code.
cd ./gmp-source &&                               # Go to the folder with gmp source code
sudo ./configure --prefix=/usr --enable-cxx &&   # Configuring the library for our computer.
sudo make &&                                     # Build gmp from source code 
cp ./.libs/libgmp.a  ./../  &&                   # Copy built libraries to project folder. This gives a significant speed up. 
cp ./.libs/libgmpxx.a  ./../  &&                 # 
cd ./../  &&                                     # Go to the project folder.
mkdir build                                      # Create the build folder.
cd build/ &&                                     # Go into it
cmake -DCMAKE_BUILD_TYPE=Release ./../ &&        # Specify to cmake that the build will be done for release, this also speeds up the application.
cmake --build . &&                               # Build the application
cp ./bin/a.out ./../ &&                          # Copy the executable file to the project folder
cd ./../                                         # Go to the project folder


CMakeLists.txt
------------------------------------------------------
                                                # This file describes the rules for cmake.
cmake_minimum_required(VERSION 3.8)             # Specify the minimum version for cmake

project(FAST_VDF)                               # Specify the name of the project

set(CMAKE_CXX_STANDARD 17)                      # Use C++ 17 standard
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(TARGET main.cpp)                            # Add cpp files to the build

include(GNUInstallDirs)                         # Add optimizing settings for GNU. it gives a speed up too.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

if(CMAKE_CXX_COMPILER_ID MATCHES GNU)            # Set c++ compiler flags for GNU
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Wall -Wno-unknown-pragmas -Wno-sign-compare -Woverloaded-virtual -Wwrite-strings -Wno-unused")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
endif()

include_directories(TARGET ${PROJECT_SOURCE_DIR}/gmp-source)
add_executable(a.out ${TARGET})

link_directories (${PROJECT_SOURCE_DIR})         # Link libraries to the project
target_link_libraries(a.out
     ${PROJECT_SOURCE_DIR}/libgmpxx.a
     ${PROJECT_SOURCE_DIR}/libgmp.a
     )


main.cpp
------------------------------------------------------
# Normalization step
# Compute r*(ar+b)+c instead of ar^2 + br + c. a*r is computed only once and used to compute b+2ar and then c+r(ar+b) accordingly. 
# Compute new C value before new B.
------------------------------------------------------
#before                                        #after

    mpz_sub(r, f.a, f.b);                       mpz_sub(r, f.a, f.b); //r=a-b

    mpz_mul_ui(denom, f.a, 2);                  mpz_mul_ui(denom, f.a, 2); //denom=2*a

    // r = (a-b) / 2ar
    mpz_fdiv_q(r, r, denom);                    mpz_fdiv_q(r, r, denom); // r = r / denom 

    mpz_set(old_b, f.b);

    mpz_mul(ra, r, f.a);                        mpz_mul(ra, r, f.a); //ra=a*r
    mpz_add(f.b, f.b, ra);                      //new C
    mpz_add(f.b, f.b, ra);                      mpz_add(ra1, ra, f.b); //ra1=a*r+b
                                                mpz_mul(ra1,ra1,r); //ra1=(a*r+b)*r
    // c += ar^2                                mpz_add(f.c, f.c, ra1); // c=c+(a*r+b)*r
    mpz_mul(ra, ra, r); 
    mpz_add(f.c, f.c, ra);                      //new B
                                                mpz_add(f.b, f.b, ra); //b=b+a*r
    // c += rb                                  mpz_add(f.b, f.b, ra); //b=b+2*a*r
    mpz_set(ra, r);
    mpz_mul(ra, ra, old_b);
    mpz_add(f.c, f.c, ra);
}
------------------------------------------------------
# Reduction step
# Compute s*(cs-b)+a instead of cs^2-bs+a. c*s is computed only once and used to compute b+2cs and then a+s(cs-b) accordingly. 
# Replace multiplication 2*cs by addition cs+cs  
------------------------------------------------------
#before                                        #after
inline void reduce(form& f) {                  inline void reduce(form& f) {
    ..........                                     ..........

    // x = 2sc
    mpz_mul(x, s, f.c);                             mpz_mul(x, s, f.c); // x=s*c
    mpz_mul_ui(x, x, 2);                            mpz_add(f.b, f.b, x); // b=-b+sc
                                                    mpz_add(f.b, f.b, x); // b=-b+2sc
    // b += 2sc
    mpz_add(f.b, f.b, x);
                                            
    // c = cs^2                                     // new C 
    mpz_mul(f.c, f.c, s);                           mpz_sub(f.c, x, old_b); //c=cs-b
    mpz_mul(f.c, f.c, s);                           mpz_mul(f.c, f.c, s); //c=(cs-b)*s
                                                    mpz_add(f.c, f.c, old_a); //c=(cs-b)*s+a
    // x = bs
    mpz_mul(x, old_b, s);                           ..............
                                                }
    // c -= bs
    mpz_sub(f.c, f.c, x);

    // c += a
    mpz_add(f.c, f.c, old_a);
        
    ..............
} 
------------------------------------------------------
# Square step
# Compute New C, New B and then New A instead of computing New A, New B, New C. It helps to remove one mpz_set call.
# Replace multiplication 2*a by addition a+a 
------------------------------------------------------
#before                                         #after
inline form* square(form &f1) {                 inline void square(form &f1) {
    ..............                                  ..............
    // New a                                        // New c
    mpz_set(old_a, f1.a);                           mpz_mul(f1.c, mu, mu); //c=mu^2 
    mpz_mul(f3.a, f1.a, f1.a);                      mpz_sub(f1.c, f1.c, x1); //c=mu^2-m

    // New b                                        // New b
    mpz_mul(a, mu, old_a);                          mpz_mul(a, mu, f1.a); // a=a*mu
    mpz_mul_ui(a, a, 2);                            mpz_sub(f1.b, f1.b, a); // b=b-a*mu
    mpz_sub(f3.b, f1.b, a);                         mpz_sub(f1.b, f1.b, a); // b=b-2*a*mu

    // New c                                        // New a
    mpz_mul(f3.c, mu, mu);                          mpz_mul(f1.a, f1.a, f1.a); //a=a^2 
    mpz_sub(f3.c, f3.c, m);
    mpz_set(f3.d, f1.d);
    reduce(f3);                                     reduce(f1);
    return &f3;
}                                               }
------------------------------------------------------
# We call square function from main instead of using repeated_square.
------------------------------------------------------
int main(int argc, char* argv[])
{
    ..............

    form x = generator_for_discriminant(&discriminant);
    for (uint64_t i=0; i < iterations; i++)
    square(x);
    ..............
}
















































