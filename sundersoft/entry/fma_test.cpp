#include "include.h"

extern "C" {
    uint64 square(void* v);
}

int main() {
    vector<char> memory;
    memory.resize(1024*1024, 0);

    char* p=&memory[0];
    while ((uint64(p)%4096)!=0) {
        ++p;
    }

    print(square(p));
}

//each input is 52 bits
/*pair<uint64, uint64> mul_int(uint64 a, uint64 b) {
    assert(a<1ull<<52);
    assert(b<1ull<<52);

    uint128 res_both=uint128(a)*uint128(b);
    uint128 res_low=res_both & ((1ull<<52)-1);
    uint128 res_high=res_both >> 52;
    return make_pair(uint64(res_high), uint64(res_low));
}

//both inputs are 52-bit unsigned integers
//output is two 52-bit integers (high and low). the high integer is shifted left by 52
//can probably make this work if the number of bits is less than 52
//
//in general for any doubles a*b+c:
// calculate the high part using fma(a,b,c)
// calculate the low part using fma(a,b,c-high); for this to work, c and high have to have the same exponent
pair<double, double> fma_mul(double a, double b, int a_shift, int b_shift) {
    {
        //check range of inputs
        int_from_double(a*d_exp2(-a_shift));
        int_from_double(b*d_exp2(-b_shift));
    }

    int shift=53+52+a_shift+b_shift;

    double res_high=fma(a, b, d_exp2(shift)); //one fma

    res_high-=d_exp2(shift); //one add
    double res_low=fma(a, b, -res_high); //x86 has a multiply-subtract instruction so this is also one fma

    //note: can use an integer add to get rid of the 52 bit shift in res_high

    {
        uint64 a_i=int_from_double(a*d_exp2(-a_shift));
        uint64 b_i=int_from_double(b*d_exp2(-b_shift));
        uint64 res_high_i=int_from_double(res_high*d_exp2(-a_shift-b_shift-52));
        uint64 res_low_i=int_from_double(res_low*d_exp2(-a_shift-b_shift));

        pair<int64, int64> res_i=mul_int(a, b);
        assert(res_high_i==res_i.first);
        assert(res_low_i==res_i.second);
    }

    return make_pair(res_high, res_low);
} */

//#include "bit_manipulation.h"
//#include "double_utility.h"

/*

uint64 int_div(uint64 a1, uint64 a2, uint64 b1, uint64 b2) {
    assert(a1_i<1ull<<52);
    assert(a2_i<1ull<<52);
    assert(b1_i==1);
    assert(b2_i<1ull<<52);

    // q=floor((a1*K + a2)/(K + b2))
    // q>1/2*a1 and q<=a1

    double b=double_from_int(b2_i)+d_exp2(52); //addition is exact

    double a1=double_from_int(a1_i)*d_exp2(52);
    double a2=double_from_int(a2_i);

    //between 0 and 2^52-1
    double q=a1/b;




    uint64 q=a1
}

int main() {
    //{
    //    uint32 b[]={0, 1, 2, 3}; //0 0 1 0 2 0 3 0
    //    uint64 a[]={1, 1};       //
    //
    //    __m128i a_simd=_mm_loadu_si128((__m128i*)a);
    //    __m128i b_simd=_mm_loadu_si128((__m128i*)b);
    //    a_simd=_mm_add_epi64(a_simd, b_simd);
    //}

    set_rounding_mode();

    while (true) {
        uint64 a; cout << "a: "; cin >> a;
        uint64 a_shift; cout << "a_shift: "; cin >> a_shift;
        uint64 b; cout << "b: "; cin >> b;
        uint64 b_shift; cout << "b_shift: "; cin >> b_shift;

        double d_a=double_from_int(a)*d_exp2(a_shift);
        double d_b=double_from_int(b)*d_exp2(b_shift);

        double d_res=d_a/d_b;
        double_bits(d_res).output(cout);
        cout << "\n";

        //pair<double, double> r=fma_mul(double_from_int(a), double_from_int(b), 0, 0);
        //cout << int_from_double(r.first*d_exp2(-52)) << ", " << int_from_double(r.second) << "\n";

        //double v1=double_from_int(a);
        //int_from_double(v1);
    }

    //double v=1.5;
    //double_bits(v).output(cout);
    //cout << "\n";

    //double_bits b;
    //b.sign=true;
    //b.set_exponent(10);
    //cout << b.to_double() << "\n";
} */