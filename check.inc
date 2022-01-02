#ifdef PREPROCESSED
/*
    Set when the test suite is coalesced into a single .c file, to be fed to
    online compilers that do not support multiple files.
    The UNEXPANDED macro avoids expanding the system #includes,
    streamlining the resulting file. This, of course, may not compile
    normally.
*/
#define UNEXPANDED(x) x
UNEXPANDED(#include <stddef.h>)
UNEXPANDED(#include <stdio.h>)
UNEXPANDED(#include <stdlib.h>)
UNEXPANDED(#include <math.h>)
UNEXPANDED(#include <time.h>)
#define dbcf_index ptrdiff_t
#define dbcf_malloc(n) malloc((size_t)(n))
#define dbcf_free(p)   free(p)
#else
#include <stdio.h>
#include <math.h>
#include <time.h>
#endif

#define DBC_FFT_IMPLEMENTATION
#include "dbc_fft.h"

#if defined(__GNUC__) && (defined(__cplusplus) || (__STDC_VERSION__>=199901L)) && defined(USE_FLOAT128)
/* Custom type example: float128. */
#define DBC_FFT_INSTANTIATION
#define DBCF_Type __float128
#define DBCF_Id q
/*
    Unfortunately, produces "non-standard suffix on floating constant"
    under -Wpedantic with no obvious way to disable it. __extension__
    only helps for C, but not C++.
*/ 
#define DBCF_LITERAL(x) (__extension__ x##q)
#include "dbc_fft.h"
#undef DBCF_Type
#undef DBCF_Id
#undef DBCF_LITERAL
#undef DBC_FFT_INSTANTIATION
#endif /* defined(__GNUC__) && (defined(__cplusplus) || (__STDC_VERSION__>=199901L)) && defined(USE_FLOAT128) */

#if (__cplusplus>=201103L) && defined(USE_FIXEDPOINT)
/* Custom type example: fixed-point. */

/*
    The NO_LITERALS version avoids handling the floating-point-like
    literals for fixed-point. In this case DBCF_LITERAL can be left
    undefined. This ensures that FFT code is purely integer (in this
    case fixed<->double conversions are only needed for test code,
    not FFT itself). Custom cexpm* functions are needed to make this
    work, since the default ones do use literals.
*/
#define NO_LITERALS

#define FRACTION_BITS 16

template<int fraction_bits>
class fixed
{
public:
    static fixed from_double(double x) {fixed ret;ret.value=(int)(x*double(1<<fraction_bits));return ret;}
    double to_double() const {return double(value)/double(1<<fraction_bits);}
    static inline fixed from_bits(int bits) {fixed ret; ret.value=bits; return ret;}
    inline int get_bits() const {return value;}
    friend inline fixed operator+(const fixed &l,const fixed &r) {fixed ret;ret.value=l.value+r.value;return ret;}
    friend inline fixed operator-(const fixed &l,const fixed &r) {fixed ret;ret.value=l.value-r.value;return ret;}
    friend inline fixed operator*(const fixed &l,const fixed &r) {fixed ret;ret.value=int(((long long)l.value*(long long)r.value)>>fraction_bits);return ret;}
    friend inline fixed operator/(const fixed &l,const fixed &r) {fixed ret;ret.value=(r.value==0?0:int((long long)l.value*(1LL<<fraction_bits)/(long long)r.value));return ret;}
    friend inline fixed operator-(const fixed &src) {fixed ret;ret.value=-src.value;return ret;}
    friend inline bool operator!=(const fixed &l,const fixed &r) {return l.value!=r.value;}
private:
    int value;
};

typedef fixed<FRACTION_BITS> fix;

#ifdef NO_LITERALS
typedef int static_assert_fraction_bits[FRACTION_BITS>=1&&FRACTION_BITS<=28?1:-1];

/* This is an example, not a high quality implementation. */
/* Compute exp(2*pi*i/n)-1. */
static void cexpm1_x(dbcf_index log2n,fix *real,fix *imag)
{
    int LUT[8][2]={
        {         0,         0},
        { 536870912,         0},
        { 268435456, 268435456},
        {  78622925, 189812531},
        {  20433432, 102725802},
        {   5157912,  52369160},
        {   1292590,  26311276},
        {    323342,  13171504}
    };
    int P=28-FRACTION_BITS;
    int M=1<<FRACTION_BITS;
    if(log2n<8) {*real=-fix::from_bits(LUT[log2n][0]>>P);*imag=fix::from_bits(LUT[log2n][1]>>P);return;}
    if(P+log2n>30) {*real=fix::from_bits(0);*imag=fix::from_bits(0);return;}
    fix x=fix::from_bits(1686629713>>(P+log2n));
    fix x2=x*x;
    *real=-x2*fix::from_bits(M/2);
    *imag=x-x*x2*fix::from_bits(M/6);
}

/* This is an example, not a high quality implementation. */
/* Compute exp(2*pi*i*(p/q))-1. */
static void cexpm1_npot_x(dbcf_index p,dbcf_index q,fix *real,fix *imag)
{
    int P=28-FRACTION_BITS;
    int M=1<<FRACTION_BITS;
    fix x=fix::from_bits(int(((1686629713LL*p))/(((long long)q)<<P)));
    fix x2=x*x;
    *real=-x2*(fix::from_bits(M/2)-x2*(fix::from_bits(M/24)-x2*(fix::from_bits(M/720)-x2*fix::from_bits(M/40320))));
    *imag=x-x*x2*(fix::from_bits(M/6)-x2*(fix::from_bits(M/120)-x2*fix::from_bits(M/5040)));
}
#endif

#define DBC_FFT_INSTANTIATION
#define DBCF_Type fix
#define DBCF_Id x
#ifdef NO_LITERALS
#undef DBCF_ZERO
#undef DBCF_ONE
#define DBCF_ZERO (fix::from_bits(0))
#define DBCF_ONE  (fix::from_bits(1<<FRACTION_BITS))
#define DBCF_LITERAL(x) /* Compile-time error if actually used. */
#define DBCF_cexpm1      cexpm1_x
#define DBCF_cexpm1_npot cexpm1_npot_x
#else
#define DBCF_LITERAL(x) (fix(x))
#endif /* NO_LITERALS */
#include "dbc_fft.h"
#undef DBCF_Type
#undef DBCF_Id
#undef DBCF_LITERAL
#ifdef NO_LITERALS
#undef DBCF_cexpm1
#undef DBCF_cexpm1_npot
#undef DBCF_ZERO
#undef DBCF_ONE
#endif /* NO_LITERALS */
#undef DBC_FFT_INSTANTIATION
#endif /* (__cplusplus>=201103L) && defined(USE_FIXEDPOINT) */

static double get_cpu_time()
{
    return (double)clock()/(double)CLOCKS_PER_SEC;
}

/* Bob Jenkins's small PRNG: http://burtleburtle.net/bob/rand/smallprng.html . */
typedef struct RNG {unsigned a,b,c,d;} RNG;

static unsigned RNG_generate(RNG *x)
{
    #define rot(x,k) (((x)<<(k))|((x)>>(32-(k))))
    unsigned e = x->a - rot(x->b, 27);
    x->a = x->b ^ rot(x->c, 17);
    x->b = x->c + x->d;
    x->c = x->d + e;
    x->d = e + x->a;
    return x->d;
    #undef rot
}

static void RNG_init(RNG *x,unsigned seed)
{
    unsigned i;
    x->a = 0xf1ea5eed, x->b = x->c = x->d = seed;
    for (i=0; i<20; ++i)
        (void)RNG_generate(x);
}

static int use_mflops=1;

#ifdef PREPROCESSED
#define MAXB ((1<<24)+64)
#else
#define MAXB ((1<<26)+64)
#endif

#ifndef DBC_FFT_NO_SIMD
DBCF_ALIGNED(64)
#endif
union Data
{
    float       buf_f[MAXB/sizeof(float)];
    double      buf_d[MAXB/sizeof(double)];
    long double buf_l[MAXB/sizeof(long double)];
#if defined(__GNUC__) && (defined(__cplusplus) || (__STDC_VERSION__>=199901L)) && defined(USE_FLOAT128)
    __float128  buf_q[MAXB/sizeof(__float128)];
#endif
#if (__cplusplus>=201103L) && defined(USE_FIXEDPOINT)
    fix buf_x[MAXB/sizeof(fix)];
#endif
    struct Mixed
    {
        float       bf[MAXB/8/sizeof(float)];
        double      bd[MAXB/8/sizeof(float)];
        long double bl[MAXB/8/sizeof(float)];
    } mixed;
} data;

static dbcf_index bitreverse_table[MAXB/sizeof(float)];

static dbcf_index bitreverse_bruteforce(dbcf_index i,dbcf_index k)
{
    dbcf_index ret=0;
    dbcf_index j=0;
    for(j=0;j<k;++j)
    {
        ret=(ret<<1)^(i&1);
        i=i>>1;
    }
    return ret;
}

#if !defined(__cplusplus) && (__STDC_VERSION__<199901L)
#define cosf(x) ((float)cos((double)(x)))
#define sinf(x) ((float)sin((double)(x)))
#define atanf(x) ((float)atan((double)(x)))
#define sqrtf(x) ((float)sqrt((double)(x)))
#define cosl(x) ((long double)cos((double)(x)))
#define sinl(x) ((long double)sin((double)(x)))
#define atanl(x) ((long double)atan((double)(x)))
#define sqrtl(x) ((float)sqrt((double)(x)))
#endif
#define cosd cos
#define sind sin
#define atand atan

#define CONCAT1(l,r) l##r
#define CONCAT(l,r) CONCAT1(l,r)
#define NAME(name) CONCAT(name,Suffix)
#define NAME2(l,r) CONCAT(l,CONCAT(Suffix,r))

#ifdef __cplusplus
template<typename T,typename S>
static inline T cast(const S &v) {return (T)v;}

template<> inline double cast<double,fix>(const fix &v) {return v.to_double();}
template<> inline fix cast<fix,double>(const double &v) {return fix::from_double(v);}
template<> inline dbcf_index cast<dbcf_index,fix>(const fix &v) {return (dbcf_index)v.to_double();}
template<> inline fix cast<fix,dbcf_index>(const dbcf_index &v) {return fix::from_double((double)v);}

#define CAST(Type,value) (cast<Type>(value))
#else
#define CAST(Type,value) ((Type)(value))
#endif

#define Type float
#define Suffix f
#include "test.inc"
#undef Type
#undef Suffix

#define Type double
#define Suffix d
#include "test.inc"
#undef Type
#undef Suffix

#define Type long double
#define Suffix l
#include "test.inc"
#undef Type
#undef Suffix

#if defined(__GNUC__) && (defined(__cplusplus) || (__STDC_VERSION__>=199901L)) && defined(USE_FLOAT128)
#include <quadmath.h>
#define Type __float128
#define Suffix q
#include "test.inc"
#undef Type
#undef Suffix
#endif

#if (__cplusplus>=201103L) && defined(USE_FIXEDPOINT)
#define STRINIFY1(x) #x
#define STRINIFY(x) STRINIFY1(x)
#define FIXED_POINT_NAME "fixed<" STRINIFY(FRACTION_BITS) ">"
static fix cosx (fix x) {return fix::from_double(cos (x.to_double()));}
static fix sinx (fix x) {return fix::from_double(sin (x.to_double()));}
static fix atanx(fix x) {return fix::from_double(atan(x.to_double()));}
#define Type fix
#define Suffix x
#include "test.inc"
#undef Type
#undef Suffix
#endif

#ifndef FIXED_POINT_NAME
#define FIXED_POINT_NAME ""
#endif

static void test_accuracy(dbcf_index maxn)
{
    dbcf_index i,q,k,n;
    dbcf_index MAX=MAXB/sizeof(float)/16;
    if(maxn<MAX) MAX=maxn;
    for(q=0;q<2;++q)
    {
        dbcf_index a=5,b=8;
        for(n=(q?a:1);n<=MAX;a+=b,b+=a,n=(q?a:n*2))
        {
            long double Ef=0.0,Ed=0.0;
            double m;
            float       *bf=data.mixed.bf;
            double      *bd=data.mixed.bd;
            long double *bl=data.mixed.bl;
            if(n<2) m=1.0;
            else    m=log((double)n)/log(2.0);
            for(k=0;k<20;++k)
            {
                long double r,e;
                r=0.0L;
                e=0.0L;
                generate_f(n*37+k*17,n,bf,bf+n);
                for(i=0;i<n;++i)
                {
                    bl[i  ]=bf[i  ];
                    bl[n+i]=bf[n+i];
                }
                dbc_fft_fc(n,bf,bf+n,bf,bf+n,1.0f);
                dbc_fft_lc(n,bl,bl+n,bl,bl+n,1.0L);
                for(i=0;i<n;++i)
                {
                    r+=bl[i  ]*bl[i  ];
                    r+=bl[n+i]*bl[n+i];
                    e+=(bf[i  ]-bl[i  ])*(bf[i  ]-bl[i  ]);
                    e+=(bf[n+i]-bl[n+i])*(bf[n+i]-bl[n+i]);
                }
                e/=r;
                if(e>Ef) Ef=e;
                r=0.0L;
                e=0.0L;
                generate_d(n*37+k*17,n,bd,bd+n);
                for(i=0;i<n;++i)
                {
                    bl[i  ]=bd[i  ];
                    bl[n+i]=bd[n+i];
                }
                dbc_fft_dc(n,bd,bd+n,bd,bd+n,1.0 );
                dbc_fft_lc(n,bl,bl+n,bl,bl+n,1.0L);
                for(i=0;i<n;++i)
                {
                    r+=bl[i  ]*bl[i  ];
                    r+=bl[n+i]*bl[n+i];
                    e+=(bd[i  ]-bl[i  ])*(bd[i  ]-bl[i  ]);
                    e+=(bd[n+i]-bl[n+i])*(bd[n+i]-bl[n+i]);
                }
                e/=r;
                if(e>Ed) Ed=e;
            }
            Ef=sqrtl(Ef);
            Ed=sqrtl(Ed);
            printf("%10.0f| %10.3e | %7.3f | %10.3e | %7.3f\n",(double)n,(double)Ef,(double)Ef*pow(2.0,23.0)/m,(double)Ed,(double)Ed*pow(2.0,52.0)/m);
        }
    }
}

int main()
{
    int simd_flags;
    static const char *types[5]={"float","double","long double","__float128",FIXED_POINT_NAME};
    dbc_fft_fi(0,0,0,0.0f); /* Initialize if neccessary . */
    simd_flags=dbcf_detect_simd();
    printf("Detected SIMD: \n");
    printf("        |x2 |x4 |x8 |x16\n");
    printf("float:  ");
    printf("| %s ",(0                          ?"+":"-"));
    printf("| %s ",(simd_flags&DBCF_HAS_SIMD4F ?"+":"-"));
    printf("| %s ",(simd_flags&DBCF_HAS_SIMD8F ?"+":"-"));
    printf("| %s ",(simd_flags&DBCF_HAS_SIMD16F?"+":"-"));
    printf("\n");
    printf("double: ");
    printf("| %s ",(simd_flags&DBCF_HAS_SIMD2D ?"+":"-"));
    printf("| %s ",(simd_flags&DBCF_HAS_SIMD4D ?"+":"-"));
    printf("| %s ",(simd_flags&DBCF_HAS_SIMD8D ?"+":"-"));
    printf("| %s ",(0                          ?"+":"-"));
    printf("\n");
    printf("\n");
    if(0)
    {
        printf("Testing permutation speed.\n");
        printf("Time is in ns/element.\n");
        printf("        %s:\n",types[0]);
        test_bitreversal_permutation_f();
        printf("        %s:\n",types[1]);
        test_bitreversal_permutation_d();
        printf("        %s:\n",types[2]);
        test_bitreversal_permutation_l();
#if defined(__GNUC__) && (defined(__cplusplus) || (__STDC_VERSION__>=199901L)) && defined(USE_FLOAT128)
        printf("        %s:\n",types[3]);
        test_bitreversal_permutation_q();
#endif
#if (__cplusplus>=201103L) && defined(USE_FIXEDPOINT)
        printf("        %s:\n",types[4]);
        test_bitreversal_permutation_x();
#endif
        printf("\n");
    }
    if(1)
    {
        printf("Testing dbc_fft.\n");
        printf("Elements of X are in [-0.5;+0.5], uniformly distributed.\n");
        if(use_mflops) printf("Speed is in Cooley-Tukey gigaflops (CTGs): 5*N*log2(N)/(time in ns),\nsimilar to FFTW benchmarks.\n");
        else           printf("Time is in ns/(5*N*log2(N)), similar to FFTW benchmarks.\n");
        printf("SoA is separate real/imag, AoS is interleaved real/imag.\n");
        /*
            See the following benchmarks:
            http://www.fftw.org/speed/index.html
            https://github.com/project-gemmi/benchmarking-fft
        */
        printf("        %s:\n",types[0]);
        test_fft_f(MAXB/sizeof(float)/8);
        printf("        %s:\n",types[1]);
        test_fft_d(MAXB/sizeof(double)/8);
        printf("        %s:\n",types[2]);
        test_fft_l(MAXB/sizeof(long double)/8);
#if defined(__GNUC__) && (defined(__cplusplus) || (__STDC_VERSION__>=199901L)) && defined(USE_FLOAT128)
        printf("        %s:\n",types[3]);
        test_fft_q(2048);
#endif
#if (__cplusplus>=201103L) && defined(USE_FIXEDPOINT)
        printf("        %s:\n",types[4]);
        test_fft_x(2048);
#endif
        printf("\n");

    }
    if(1)
    {
        printf("Testing accuracy.\n");
        printf("Values reported are:\n");
        printf("    Err=RMS(error)/RMS(output),\n");
        printf("    Rel=Err/(E*log2(N)),\n");
        printf("where E=ULP(1), i.e. 2^{-23} for float and 2^{-52} for double.\n");
        if(sizeof(long double)<=sizeof(float))
            printf("WARNING: accuracy figures for float are likely wrong.\n");
        if(sizeof(long double)<=sizeof(double))
            printf("WARNING: accuracy figures for double are likely wrong.\n");
        printf("          |         float        |         double       \n");
        printf("        N |    Err     |   Rel   |    Err     |   Rel   \n");
        printf("----------+------------+---------+------------+---------\n");
        test_accuracy(MAXB/sizeof(float)/16/32);
        printf("\n");
    }
    return 0;
}
