/*
    dbc_fft.h - public domain single-file FFT library.
    Version: alpha.2.
    This is intended to be easy to use and reasonably efficient
    complex FFT implementation for arbitrary (e.g. non-power-of-2) input
    size.
    DBC stands for "Designed by Committee".

DEFINITIONS
    Below, the (forward) discrete Fourier transform (DFT) of a sequence
    X of N complex numbers (X[0] .. X[N-1]) is defined as:
        Y[j]=scale*sum(X[k]*exp(-2*pi*i*j*k)/N,0<=k<N),
    where i is the imaginary unit, and scale is a (real) scaling factor
    of choice.
    Similarly, the inverse discrete Fourier transform (IDFT) is defined as:
        X[k]=scale*sum(Y[j]*exp(+2*pi*i*j*k)/N,0<=j<N),
    with a possibly different scale.
    If scale=1 for both, then IDFT(DFT(X)))=DFT(IDFT(X))=N*X. Popular choices
    are (1:1/N) (Wikipedia's definition), (1/N:1) (DTFS), and
    (1/sqrt(N):1/sqrt(N)) (unitary DFT) for DFT and IDFT respectively.
    In dbc_fft.h you simply provide the desired scale explicitly.
    Fast Fourier transform (FFT) is a name of a class of efficient
    algorithms for computing DFT and IDFT.

USAGE
    The dbc_fft.h follows the general style of stb libraries
    ( https://github.com/nothings/stb ), which may help for those already
    familiar with them. Included normally, this file acts as a header.
    To create the implementation
        #define DBC_FFT_IMPLEMENTATION
    in *one* C/CPP file that includes this file. To make the implementation
    private to the file that generates it
        #define DBC_FFT_STATIC

    The function
        int dbc_fft_fc(
            dbcf_index num_elements,
            const float *src_real,const float *src_imag,
                  float *dst_real,      float *dst_imag,
            float scale);
    computes the FFT of num_elements complex numbers, whose real and imaginary
    parts are stored in contiguous arrays src_real and src_imag respectively.

    * src_real and/or src_imag can be NULL, which is treated as array of all
    zeros.
    * dst_real and dst_imag shall not be NULL.
    * src_real==dst_real and/or src_imag==dst_imag are allowed, to support
    in-place FFT's. Other overlap (e.g. partial overlap, or src_real==dst_imag)
    is not.
    * dbcf_index is simply an integer type the library uses for most of its
    size and index calculations. You can #define it to whatever you want
    (e.g. int or size_t; however, signed types are recommended - for one,
    they allow negative strides). If you don't, the library uses ptrdiff_t
    by default.
    * function returns zero on success and non-zero on error (e.g. size not
    being a power of 2).

    If you want the real and imaginary parts interleaved instead, the function
        int dbc_fft_fi(
            dbcf_index num_elements,
            const float *src,
                  float *dst,
            float scale);
    can be used.
    Note: the C++ Standard explicitly allows to reinterpret array of
    std::complex as an interleaved array twice the size of the underlying
    floating point type. The C Standard seems to suggest the same for
    its _Complex types.
    Again, src can be NULL, dst shall not, and src==dst is allowed,
    but other overlap is not.
    Note, that this function may be somewhat slower, than the contiguous
    version.

    More general, strided version (which the above are essentially wrappers of)
    of the function is:
        int dbc_fft_fs(
            dbcf_index num_elements,
            const float *src_real,const float *src_imag,
            dbcf_index src_real_stride,dbcf_index src_imag_stride,
                  float *dst_real,      float *dst_imag,
            dbcf_index dst_real_stride,dbcf_index dst_imag_stride,
            float scale);
    This function can be convenient, if e.g. performing 2D FFT, and
    want stride to be size of the row.
    Again, src can be NULL, dst shall not, and src==dst is allowed,
    but other overlap is not. src_*_stride can be 0, but not dst_*_stride.
    The function may be slower than contiguous version, especially for
    dst strides other than 1 (src strides have less impact).

    There are also corresponding IFFT versions (dbc_ifft_fc, dbc_ifft_fi,
    dbc_ifft_fs).

    These functions also exist in versions for double and long double.
    They have 'd' and 'l' instead of 'f' in their suffixes (e.g.
    dbc_fft_dc, dbc_fft_li).
    Generating code for specific type can be suppressed by defining
    DBC_FFT_NO_FLOAT, DBC_FFT_NO_DOUBLE, DBC_FFT_NO_LONGDOUBLE.

    For C++ all of the above (3 functions x 3 types) are available as
    overloads of dbc_fft and dbc_ifft, unless DBC_FFT_NO_CPP_OVERLOADS is
    defined.

ACCURACY
    Experimentally, the average error is estimated as:
        RMS(Error)<C*E*RMS(Output)*log2(N)
    where E is ULP(1) for the specific type (i.e. 2^{-23} for float,
    2^{-52} for double, etc.), and C=0.5 for power-of-2 sizes, and C=1 for
    non-power-of-2 sizes.

PERFORMANCE
    Performance of FFT is commonly measured in Cooley-Tukey gigaflops (CTGs),
    which are defined (for complex FFT) as CTG=5*N*log2(N)/(time in ns).
    On the test machine (Intel(R) Xeon(R) Platinum 8124M CPU @ 3.00GHz),
    the peak performance of bdc_fft.h is around 14 CTGs for float,
    and around 10 CTGs for double (both around N=4096), with N=128 being about
    twice slower. Interleaved output is about twice slower than contiguous
    for larger N in optimized (SIMD) case (not much difference otherwise).
    Non-power-of-2 case is about 10 times slower than power-of-2 of similar
    size.
    The benchmark was compiled as C++ using GCC 11.2 with -m64 -march=native -O3
    on Linux. -O2 is slightly (and -O1/-Os significantly) slower. -march=native
    helps somewhat. Observed slowdown if compiled as C in some situations.

THREAD SAFETY
    The library should be thread-safe in a sense that computing 2 distinct
    FFTs in different threads should cause no problems. WARNING: by default
    there is a performance penalty on x86 with SIMD enabled, unless the
    code is compiled as C++11 or newer. You can avoid it by
#define DBC_FFT_CACHE_CPU_DETECTION
    which makes the first call to dbc_fft not thread-safe. The first call
    may be e.g. dbc_fft_fi(0,0,0,0.0f);. Yes, this is lame.
    There is also
#define DBC_FFT_DONT_CACHE_CPU_DETECTION
    to disable it even in C++11.
    If you are using custom dbcf_detect_simd(), dbcf_malloc(), dbcf_free()
    you are responsible for their thread safety.
    The library itself uses no threading internally. 

CUSTOM TYPES
    dbc_fft.h uses #include __FILE__ trick to instantiate the FFT code for
    the concrete types. As a bonus, it is possible, with little effort,
    to instantiate an implementation for your own type, provided it has the
    usual arithmetic operations and a way to construct it from a
    floating-point-like literal. For example, the following provides (in
    addition to float/double/long double) an FFT implementation
    (dbc_fft_qc, etc.) for GCC's quadruple precision
    floats ( https://gcc.gnu.org/onlinedocs/gcc/Floating-Types.html ):
#define DBC_FFT_IMPLEMENTATION
#include "dbc_fft.h"
#define DBC_FFT_INSTANTIATION
#define DBCF_Type __float128
#define DBCF_Id q
// To compile as C++ may need -fext-numeric-literals switch.
// Also, produces "non-standard suffix on floating constant" warnings
// under -Wpedantic. Unfortunately, __extension__ only seems to
// help in C, but not in C++.
#define DBCF_LITERAL(x) (__extension__ x##q)
#include "dbc_fft.h"
#undef DBCF_Type 
#undef DBCF_Id 
#undef DBCF_LITERAL
#undef DBC_FFT_INSTANTIATION
    If you need ONLY the custom type, you may define DBC_FFT_NO_FLOAT,
    DBC_FFT_NO_DOUBLE, DBC_FFT_NO_LONGDOUBLE before the first #include.
    The internal complex exponent calculation routines are only accurate
    enough to about float128's precision. If the type is more precise,
    you may want to provide replacements:
#define DBCF_cexpm1      ...implementation...
#define DBCF_cexpm1_npot ...implementation...
    See code for their actual semantics and signatures. Another reason
    to provide them is to avoid dealing with literals altogether,
    as these are the only functions that use non-trivial literals (this
    may matter e.g. for fixed-point types). In this case you also need
    to provide the trivial constants:
#define DBCF_ZERO value
#define DBCF_ONE  value
    It is recommended to #undef those 4 after #include.
    If you really know what you are doing, you can also provide an optimized
    pass implementation (DBCF_butterfly_multipass_optimized).

MEMORY USAGE
    The heap ("dynamic") memory allocation only happens for non-power-of-2
    sizes. Memory is allocated/freed via the dbcf_malloc()/dbcf_free() calls,
    which you can #define to your own implementations. At most 10 times the
    size of output is allocated. You can
#define DBC_FFT_NO_NPOT
    to disable the non-power-of-2 code entirely (the call to fft functions
    with non-power-of-2 size will return DBCF_ERROR_INVALID_ARGUMENT in
    this case).
    For power-of-2 sizes no heap memory allocation whatsoever occurs.
    A few modest tables (less than 2KB in total for float+double+long double)
    are statically allocated. They can be disabled by
#define DBC_FFT_NO_BITREVERSE_TABLE // Saves 512 bytes.
#define DBCF_cexpm1 ...implementation... // Saves 34*sizeof(type) bytes per type.
    Note: DBCF_cexpm1 needs to be defined separately for each type. You may need
    to #include "dbc_fft.h" several times after #define DBC_FFT_IMPLEMENTATION,
    with different DBC_FFT_NO_FLOAT/etc. settings.
    A call to fft function allocates one modest temporary buffer on the stack.
    You can set its size by
#define DBCF_TMP_BUF_LOG2 value
    The default is 10 (resulting in 1024*sizeof(type) bytes). You can set it
    as low as 2, but that noticeably degrades performance compared
    to 4 (which itself is somewhat slower than the default). Increasing it
    may improve the performance slightly for large inputs.
    Also, a small amount of space (O(log(N))) on stack is used for recursion.
    Since dbc_fft can work inplace, the separate destination buffer might
    not be neccessary.

CODE FOOTPRINT
    At minimum (float only, no SIMD, no nPOT, no tables, compiled
    with -Os -s), dbc_fft adds around 4 KB to the binary size, as measured
    on x86. The size is around 128 KB on the default settings.

SIMD
    By default, dbc_fft.h tries to use SIMD where available. On x86/x64 it
    tries to automatically detect it at runtime. On other platforms SIMD is
    only available for GCC-style compilers (GCC, clang, ICC), and needs to
    be explicitly enabled at compile-time (no runtime detection is performed)
    for all desired type/width combinations, e.g.:
#define DBC_FFT_IMPLEMENTATION
#define DBC_FFT_FORCE_SIMD (DBCF_HAS_SIMD4F|DBCF_HAS_SIMD2D)
#inlcude "dbc_fft.h"
    For ARM NEON (if detected at compile time) the dbc_fft.h does this
    internally (with exactly DBCF_HAS_SIMD4F|DBCF_HAS_SIMD2D).
    Supported flags are DBCF_HAS_SIMD{4|8|16}F for float and
    DBCF_HAS_SIMD{2|4|8}D for double. You may also need to specify
    the complier options to actually enable the instructions in question.
    This can also be used on x86/x64 and/or NEON to override the default
    detection (you get exactly what you requested, and no runtime detection
    is performed).
    If you want to, you can provide a replacement CPU detection function:
#define dbcf_detect_simd() ...implementation...
    The expected signature is int dbcf_detect_simd(void), the return value is
    a bitmask of DBCF_HAS_SIMD* flags. It gets called several times per FFT,
    so you may want to cache the result.
    WARNING: Using dbcf_detect_simd to get the runtime SIMD detection
    on non-x86 platform (with GCC-style compiler) is NOT GUARANTEED TO WORK.
    If you want to try anyway, you probably need to provide the appropriate
    ABI declarations for all relevant type/width combinations, along
    the lines of:
// Targeting the hypothetical SuperLame CPU.
#define DBC_FFT_IMPLEMENTATION
#define DBC_FFT_FORCE_SIMD (DBCF_HAS_SIMD4F|DBCF_HAS_SIMD2D) // Enable generating code for SIMD functions.
#define dbcf_detect_simd() (__builtin_cpu_supports("super_lame_simd")?(DBCF_HAS_SIMD4F|DBCF_HAS_SIMD2D):0)
#define DBCF_DECL_SIMD4F __attribute__((target("super_lame_simd")))
#define DBCF_DECL_SIMD2D __attribute__((target("super_lame_simd")))
#inlcude "dbc_fft.h"
    WARNING: using SIMD in general, and runtime CPU detection specifically,
    can in some cases lead to problems, e.g.:
https://www.virtualdub.org/blog2/entry_363.html
https://gist.github.com/rygorous/f26f5f60284d9d9246f6
https://github.com/nothings/stb/issues/280
https://github.com/nothings/stb/issues/410
    To disable SIMD you can
#define DBC_FFT_NO_SIMD
    You can also selectively disable some instructions via
#define DBC_FFT_NO_AVX
#define DBC_FFT_NO_AVX512
    WARNING: By default SIMD is disabled for MinGW on 32-bit targets,
    because of the problems it causes:
https://www.peterstock.co.uk/games/mingw_sse/
https://github.com/nothings/stb/issues/81
    To enable it
#define DBC_FFT_ENABLE_MINGW_SIMD
    In this case you are expected to take measures to make it actually compile
    correctly (e.g. -mstackrealign). Note: DBC_FFT_FORCE_SIMD also overrides
    this, but you don't get runtime detection that way.
    WARNING: Some older OSes (Windows 95 and earlier, Linux kernel before
    something like 2.4) may not have the OS-level support for SSE (do not
    preserve the state on context switch). The dbc_fft.h does not check
    for it. If you need to support them, consider checking it yourself and
    providing a custom dbcf_detect_simd.

ALGORITHM
    The implementation details are documented below.
    For power-of-2 sizes a simple radix-2 decimation-in-time Cooley–Tukey FFT
    algorithm ( https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm )
    is used. The FFT function performs an explicit bit-reversal, followed by
    the butterflies. The bit-reversal is the only step that touches src.
    Both in-place and out-of-place cases are handled. Small (<=256 elements)
    inputs are reversed directly. Medium-sized inputs do recursion:
    reverse_and_swap(0XXX1,1XXX0), reverse(0XXX0), reverse(1XXX1). Large inputs
    are processed with the algorithm from "Towards an Optimal Bit-Reversal
    Permutation Program" by Larry Carter and Kang Su Gatlin. The out-of-place
    case transfers data to dst, and calls the in-place case, which (somewhat
    surprisingly) turns out to be faster (especially for smaller Q).
    The butterflies are done either as a number of passes over the entire input,
    or recursively on two halves, followed by a single butterfly pass (for larger
    inputs), to improve locality. The bottom passes are combined into
    hand-written FFT8.
    Interleaved dst is temporarily deinterleaved to make better use of SIMD
    (only when SIMD is enabled). The (de)interleave is simply bit-reversal
    permutation on the array and its halves.
    All twiddle factors are calculated from O(log(N)) values W^{2^k} (where
    W is the twiddle for the smallest angle), computed by complex exponent
    routines (which simply use precomputed tables or Taylor series), ensuring
    at most O(log(N)) arithmetic operations (and, therefore, roundoff errors)
    per twiddle. Actual calculations use exp(i*z)-1 rather than exp(i*z), for
    somewhat improved accuracy.
    Twiddles are written to a fixed-size buffer (allocated on the stack),
    which is then supplied to the butterfly passes. If the buffer is too
    small, the function recurses, and additional multipliers are supplied on top
    of buffer (so, roughly, twiddle[i]=multiplier*buffer[i%BUFFER_SIZE]).
    For all non-power-of-2 sizes Bluestein's algorithm is used, even though
    for some sizes different algorithms may be significantly more efficient.

LICENSE
    This software is dual-licensed to the public domain and under the following
    license: you are granted a perpetual, irrevocable license to copy, modify,
    publish, and distribute this file as you see fit.
*/

/*============================================================================*/
/* Header section. */
#ifndef DBC_FFT_H
#define DBC_FFT_H

#ifndef dbcf_index
#include <stddef.h>
typedef ptrdiff_t dbcf_index;
#endif

#ifdef DBC_FFT_STATIC
#define DBCF_DEF static
#else
#define DBCF_DEF extern
#endif

#define DBCF_ERROR_INVALID_ARGUMENT (-1)
#define DBCF_ERROR_OUT_OF_MEMORY    (-2)

#define DBCF_CONCAT1(x,y) x##y
#define DBCF_CONCAT(x,y) DBCF_CONCAT1(x,y)

#define DBCF_NAME2(nameL,nameR) DBCF_CONCAT(nameL,DBCF_CONCAT(_,DBCF_CONCAT(DBCF_Id,nameR)))
#define DBCF_NAME(name) DBCF_CONCAT(name,DBCF_CONCAT(_,DBCF_Id))/*DBCF_NAME2(name,)*/

/* Declarations */
#define DBC_FFT_DECLARATION

#ifndef DBC_FFT_NO_FLOAT

#define DBCF_Type float
#define DBCF_Id f
    #include __FILE__
#undef DBCF_Type
#undef DBCF_Id

#endif /* DBC_FFT_NO_FLOAT */

#ifndef DBC_FFT_NO_DOUBLE

#define DBCF_Type double
#define DBCF_Id d
    #include __FILE__
#undef DBCF_Type
#undef DBCF_Id

#endif /* DBC_FFT_NO_DOUBLE */

#ifndef DBC_FFT_NO_LONGDOUBLE

#define DBCF_Type long double
#define DBCF_Id l
    #include __FILE__
#undef DBCF_Type
#undef DBCF_Id

#endif /* DBC_FFT_NO_LONGDOUBLE */

#undef DBC_FFT_DECLARATION

#endif /* DBC_FFT_H */

/*============================================================================*/
/* Declaration section */
#ifdef DBC_FFT_DECLARATION

#ifdef __cplusplus
extern "C" {
#endif

DBCF_DEF int DBCF_NAME2(dbc_fft,c)(
    dbcf_index num_elements,
    const DBCF_Type *src_real,const DBCF_Type *src_imag,
          DBCF_Type *dst_real,      DBCF_Type *dst_imag,
    DBCF_Type scale);

DBCF_DEF int DBCF_NAME2(dbc_ifft,c)(
    dbcf_index num_elements,
    const DBCF_Type *src_real,const DBCF_Type *src_imag,
          DBCF_Type *dst_real,      DBCF_Type *dst_imag,
    DBCF_Type scale);

DBCF_DEF int DBCF_NAME2(dbc_fft,i)(
    dbcf_index num_elements,
    const DBCF_Type *src,
          DBCF_Type *dst,
    DBCF_Type scale);

DBCF_DEF int DBCF_NAME2(dbc_ifft,i)(
    dbcf_index num_elements,
    const DBCF_Type *src,
          DBCF_Type *dst,
    DBCF_Type scale);

DBCF_DEF int DBCF_NAME2(dbc_fft,s)(
    dbcf_index num_elements,
    const DBCF_Type *src_real,const DBCF_Type *src_imag,
    dbcf_index src_real_stride,dbcf_index src_imag_stride,
          DBCF_Type *dst_real,      DBCF_Type *dst_imag,
    dbcf_index dst_real_stride,dbcf_index dst_imag_stride,
    DBCF_Type scale);

DBCF_DEF int DBCF_NAME2(dbc_ifft,s)(
    dbcf_index num_elements,
    const DBCF_Type *src_real,const DBCF_Type *src_imag,
    dbcf_index src_real_stride,dbcf_index src_imag_stride,
          DBCF_Type *dst_real,      DBCF_Type *dst_imag,
    dbcf_index dst_real_stride,dbcf_index dst_imag_stride,
    DBCF_Type scale);

#ifdef __cplusplus
}
#endif

/* Overloaded versions. */
#if defined(__cplusplus) && !defined(DBC_FFT_NO_CPP_OVERLOADS)
DBCF_DEF int dbc_fft(
    dbcf_index num_elements,
    const DBCF_Type *src_real,const DBCF_Type *src_imag,
          DBCF_Type *dst_real,      DBCF_Type *dst_imag,
    DBCF_Type scale);
DBCF_DEF int dbc_ifft(
    dbcf_index num_elements,
    const DBCF_Type *src_real,const DBCF_Type *src_imag,
          DBCF_Type *dst_real,      DBCF_Type *dst_imag,
    DBCF_Type scale);
DBCF_DEF int dbc_fft(
    dbcf_index num_elements,
    const DBCF_Type *src,
          DBCF_Type *dst,
    DBCF_Type scale);
DBCF_DEF int dbc_ifft(
    dbcf_index num_elements,
    const DBCF_Type *src,
          DBCF_Type *dst,
    DBCF_Type scale);
DBCF_DEF int dbc_fft(
    dbcf_index num_elements,
    const DBCF_Type *src_real,const DBCF_Type *src_imag,
    dbcf_index src_real_stride,dbcf_index src_imag_stride,
          DBCF_Type *dst_real,      DBCF_Type *dst_imag,
    dbcf_index dst_real_stride,dbcf_index dst_imag_stride,
    DBCF_Type scale);
DBCF_DEF int dbc_ifft(
    dbcf_index num_elements,
    const DBCF_Type *src_real,const DBCF_Type *src_imag,
    dbcf_index src_real_stride,dbcf_index src_imag_stride,
          DBCF_Type *dst_real,      DBCF_Type *dst_imag,
    dbcf_index dst_real_stride,dbcf_index dst_imag_stride,
    DBCF_Type scale);
#endif /* defined(__cplusplus) && !defined(DBC_FFT_NO_CPP_OVERLOADS) */

#endif /* DBC_FFT_DECLARATION */

/*============================================================================*/
/* Implementation section. */
#if defined(DBC_FFT_IMPLEMENTATION) && !defined(DBC_FFT_DECLARATION) && !defined(DBC_FFT_INSTANTIATION)

#define DBCF_POW2(n) (((dbcf_index)1)<<(n))

#ifndef DBCF_TMP_BUF_LOG2
#define DBCF_TMP_BUF_LOG2 10
#endif
#define DBCF_TMP_BUF_SIZE DBCF_POW2(DBCF_TMP_BUF_LOG2)
#define DBCF_TWIDDLES_BUF_LOG2 ((DBCF_TMP_BUF_LOG2)-1)
#define DBCF_TWIDDLES_BUF_SIZE DBCF_POW2(DBCF_TWIDDLES_BUF_LOG2)
#ifndef DBCF_Q /* Parameter Q from  "Towards an Optimal Bit-Reversal Permutation Program". */
#define DBCF_Q (((DBCF_TMP_BUF_LOG2)>>1)<6?((DBCF_TMP_BUF_LOG2)>>1):6)
#endif

typedef int dbcF_static_assert_tmp_buf_size[(DBCF_TMP_BUF_LOG2)>=2                        ?1:-1];
typedef int dbcF_static_assert_Q           [(DBCF_Q)>=1&&(2*(DBCF_Q)<=(DBCF_TMP_BUF_LOG2))?1:-1];

/* malloc replacement. */
#ifndef dbcf_malloc
#include <stddef.h>
#include <stdlib.h>
#define dbcf_malloc(n) malloc((size_t)(n))
#define dbcf_free(p)   free(p)
#endif

/* Architecture detection. */
#if defined(__i386__) || defined(_M_IX86) || defined(_X86_) || defined(__i386) || defined(__THW_INTEL)
#define DBCF_X86
#endif

#if defined(__x86_64__) || defined(_M_X64) || defined(__amd64)
#define DBCF_X64
#endif

#if defined(DBCF_X86) || defined(DBCF_X64)
#define DBCF_X86_OR_X64
#endif

/* 8-bit bitreverse table. */
#ifndef DBC_FFT_NO_BITREVERSE_TABLE
static unsigned char dbcF_bitreverse_table[512]=
{
/* Adapted from https://graphics.stanford.edu/~seander/bithacks.html#BitReverseTable . */
#define DBCF_R1(n,c) n,n+c
#define DBCF_R2(n,c) DBCF_R1(n,2*c),DBCF_R1(n+c,2*c)
#define DBCF_R3(n,c) DBCF_R2(n,2*c),DBCF_R2(n+c,2*c)
#define DBCF_R4(n,c) DBCF_R3(n,2*c),DBCF_R3(n+c,2*c)
#define DBCF_R5(n,c) DBCF_R4(n,2*c),DBCF_R4(n+c,2*c)
#define DBCF_R6(n,c) DBCF_R5(n,2*c),DBCF_R5(n+c,2*c)
#define DBCF_R7(n,c) DBCF_R6(n,2*c),DBCF_R6(n+c,2*c)
#define DBCF_R8(n,c) DBCF_R7(n,2*c),DBCF_R7(n+c,2*c)
0,0,DBCF_R1(0,1),DBCF_R2(0,1),DBCF_R3(0,1),DBCF_R4(0,1),DBCF_R5(0,1),DBCF_R6(0,1),DBCF_R7(0,1),DBCF_R8(0,1)
#undef DBCF_R1
#undef DBCF_R2
#undef DBCF_R3
#undef DBCF_R4
#undef DBCF_R5
#undef DBCF_R6
#undef DBCF_R7
#undef DBCF_R8
};
#endif /* DBC_FFT_NO_BITREVERSE_TABLE */

static dbcf_index dbcF_bitreverse(dbcf_index i,dbcf_index bits)
{
#ifndef DBC_FFT_NO_BITREVERSE_TABLE
    if(bits<=8) return (dbcF_bitreverse_table+DBCF_POW2(bits))[i];
    return ((dbcf_index)((dbcF_bitreverse_table+256)[i&255])<<(bits-8))^dbcF_bitreverse(i>>8,bits-8);
#else
    dbcf_index j=i&255,b=(bits<8?bits:8);
    j=(j>>4)^(j<<4);
    j=((j&0xCC)>>2)^((j&0x33)<<2);
    j=((j&0xAA)>>1)^((j&0x55)<<1);
    j=j>>(8-b);
    if(bits>8) j=(j<<(bits-8))^dbcF_bitreverse(i>>8,bits-8);
    return j;
#endif
}

/* SIMD feature detection flags. */
#define DBCF_HAS_SIMD4F   1
#define DBCF_HAS_SIMD2D   2 
#define DBCF_HAS_SIMD8F   4 
#define DBCF_HAS_SIMD4D   8
#define DBCF_HAS_SIMD16F 16 
#define DBCF_HAS_SIMD8D  32

#if (defined(__MINGW32__)||defined(__MINGW64__))&&!defined(DBCF_X64)
#if !defined(DBC_FFT_ENABLE_MINGW_SIMD) && !defined(DBC_FFT_NO_SIMD) && !defined(DBC_FFT_FORCE_SIMD)
#define DBC_FFT_NO_SIMD
#endif
#endif

#if defined(DBC_FFT_FORCE_SIMD) && defined(DBC_FFT_NO_SIMD)
#error Both DBC_FFT_FORCE_SIMD and DBC_FFT_NO_SIMD are specified
#endif

#if defined(__GNUC__) && (defined(__ARM_NEON)||defined(__ARM_NEON__)) && !defined(DBC_FFT_NO_SIMD) && !defined(DBC_FFT_FORCE_SIMD)
#define DBC_FFT_FORCE_SIMD (DBCF_HAS_SIMD4F|DBCF_HAS_SIMD2D)
#endif

#if !defined(DBCF_X86_OR_X64) && !defined(DBC_FFT_NO_SIMD) && !defined(DBC_FFT_FORCE_SIMD)
#define DBC_FFT_NO_SIMD
#endif

#if !defined(DBCF_X86_OR_X64) && !defined(__GNUC__) && defined(DBC_FFT_FORCE_SIMD)
/* DBC_FFT_FORCE_SIMD is only supported on x86/x64 or for __GNUC__ compilers. */
#error DBC_FFT_FORCE_SIMD is not supported for this platform
#endif

#if !defined(DBCF_X86_OR_X64)
#if defined(DBC_FFT_NO_AVX)
#undef DBC_FFT_NO_AVX
#endif
#if defined(DBC_FFT_NO_AVX512)
#undef DBC_FFT_NO_AVX512
#endif
#endif

#if defined(DBCF_X86_OR_X64)
#if defined(DBC_FFT_NO_AVX)
#define DBCF_NO_SIMD8F
#define DBCF_NO_SIMD4D
#endif
#if defined(DBC_FFT_NO_AVX512)
#define DBCF_NO_SIMD16F
#define DBCF_NO_SIMD8D
#endif
#endif /* defined(DBCF_X86_OR_X64) */

#if defined(DBC_FFT_FORCE_SIMD)
#if !((DBC_FFT_FORCE_SIMD)&DBCF_HAS_SIMD4F)
#define DBCF_NO_SIMD4F
#endif
#if !((DBC_FFT_FORCE_SIMD)&DBCF_HAS_SIMD8F)
#define DBCF_NO_SIMD8F
#endif
#if !((DBC_FFT_FORCE_SIMD)&DBCF_HAS_SIMD16F)
#define DBCF_NO_SIMD16F
#endif
#if !((DBC_FFT_FORCE_SIMD)&DBCF_HAS_SIMD2D)
#define DBCF_NO_SIMD2D
#endif
#if !((DBC_FFT_FORCE_SIMD)&DBCF_HAS_SIMD4D)
#define DBCF_NO_SIMD4D
#endif
#if !((DBC_FFT_FORCE_SIMD)&DBCF_HAS_SIMD8D)
#define DBCF_NO_SIMD8D
#endif
#endif /* defined(DBC_FFT_FORCE_SIMD) */

#if defined(DBC_FFT_USE_VECTOR_EXTENSIONS) && defined(DBC_FFT_USE_INTRINSICS)
#error Both DBC_FFT_USE_VECTOR_EXTENSIONS and DBC_FFT_USE_INTRINSICS are specified
#endif

#if !defined(DBC_FFT_USE_VECTOR_EXTENSIONS) && !defined(DBC_FFT_USE_INTRINSICS)
#if defined(__GNUC__) || !defined(DBCF_X86_OR_X64)
#define DBC_FFT_USE_VECTOR_EXTENSIONS
#else
#define DBC_FFT_USE_INTRINSICS
#endif
#endif

#if defined(DBC_FFT_CACHE_CPU_DETECTION) && defined(DBC_FFT_DONT_CACHE_CPU_DETECTION)
#error Both DBC_FFT_CACHE_CPU_DETECTION and DBC_FFT_DONT_CACHE_CPU_DETECTION are specified
#endif

#if !defined(DBC_FFT_CACHE_CPU_DETECTION) && !defined(DBC_FFT_DONT_CACHE_CPU_DETECTION)
#if !(__cplusplus>=201103L) 
#define DBC_FFT_CACHE_CPU_DETECTION
#else
#define DBC_FFT_DONT_CACHE_CPU_DETECTION
#endif
#endif

/* Alignment specifiers. */
#if defined(__GNUC__)
#define DBCF_ALIGNED(n) __attribute__((aligned(n)))
#elif defined(_MSC_VER)
#define DBCF_ALIGNED(n) __declspec(align(n))
#else
#define DBCF_NO_ALIGNMENT
#define DBCF_ALIGNED(n)
#endif

/* Is this reliable? */
#define DBCF_IS_ALIGNED(p,n) ((((dbcf_index)(p))&((dbcf_index)(n)-1))==0)

/* SIMD stuff. */
#if !defined(DBC_FFT_NO_SIMD)

#if !defined(DBC_FFT_FORCE_SIMD) && !defined(dbcf_detect_simd) && defined(DBCF_X86_OR_X64) && !defined(_WIN16)
/*
    x86 runtime CPU detection.

    We need to check both CPU-level support and OS-level support. This
    is because the OS needs know to preserve SIMD state on context switch.
    Generally, the corresponding instructions are disabled (cause exception
    if used) on OSes that do not support them. XGETBV is used to check
    OS-level support for AVX and AVX512. Some older OSes (Windows 95 and
    earlier, Linux kernel before something like 2.4) may not have the
    OS-level support for SSE (do not preserve the state on context switch).
    The dbc_fft.h does not check for it. If you need to support them,
    consider checking it yourself and providing a custom dbcf_detect_simd.

    SSE2 is used as a minimum SIMD target on x86, even if SSE could
    suffice for float. This is because SSE can result in slow code
    with x86 GCC ( https://godbolt.org/z/4T7fx5b7o ), as well as in x64.
*/
#if !defined(__GNUC__) && (!defined(_MSC_VER) || (_MSC_VER>=1400))
#include <intrin.h>
#endif

/* All Pentiums and later should have CPUID. Some later 486s have it as well. */
static int dbcF_has_cpuid(void)
{
#if defined(DBCF_X64)
    return 1; /* Always have CPUID on x64. */
#elif defined(__GNUC__)
    /* Check by trying to set the 21st bit in EFLAGS. */
    unsigned a,b;
    __asm__ __volatile__(
        "pushfl\n\t"
        "pushfl\n\t"
        "popl %0\n\t"
        "movl %0,%1\n\t"
        "xorl $0x00200000,%0\n\t"
        "pushl %0\n\t"
        "popfl\n\t"
        "pushfl\n\t"
        "popl %0\n\t"
        "popfl\n\t"
	:"=&r"(a),"=&r"(b));
    return !!((a^b)&0x00200000u);
#else
    return 1; /* Assume the CPU has CPUID. */
#endif
}

static void dbcF_cpuid(int level,int sublevel,unsigned *eax,unsigned *ebx,unsigned *ecx,unsigned *edx)
{
#if defined(__GNUC__)
    unsigned a,b,c,d;
    __asm__ __volatile__(
        "cpuid\n\t"
        :"=a"(a),"=b"(b),"=c"(c),"=d"(d)
        :"0"(level),"2"(sublevel));
    *eax=a;
    *ebx=b;
    *ecx=c;
    *edx=d;
#elif defined(_MSC_VER) && !(_MSC_VER>=1400) /* VC6 */
    int res;
    __asm {
        mov  eax,1
        cpuid
        mov  res,edx
    }
    /* Only check edx. */
    *eax=0;
    *ebx=0;
    *ecx=0;
    *edx=(unsigned)res;
#else
    int regs[4];
    __cpuidex(regs,level,sublevel);
    *eax=(unsigned)regs[0];
    *ebx=(unsigned)regs[1];
    *ecx=(unsigned)regs[2];
    *edx=(unsigned)regs[3];
#endif
}

/* Only return lower 32 bits of xgetbv. */
static unsigned dbcF_xgetbv(unsigned level)
{
    unsigned ret=0;
#if defined(__GNUC__)
    /*
        Using __builtin_ia32_xgetbv (or _xgetbv intrinsic) needs the function
        to be compiled as AVX or higher, but using xgetbv from within the
        inline assembly does not. Seems to cause no problems (as long as
        the availability of XGETBV is checked at runtime).
    */
    unsigned eax,edx;
    __asm__ __volatile__("xgetbv\n\t":"=a"(eax),"=d"(edx):"c"(level));
    ret=eax; /* Full result: eax^((unsigned long long)edx<<32). */
#elif defined(_MSC_VER) && !(_MSC_VER>=1400) /* VC6 */
    /* Pretend that xgetbv always returns 0. */
#else
    ret=(unsigned)_xgetbv(level);
#endif
    return ret;
}
#endif /* !defined(DBC_FFT_FORCE_SIMD) && !defined(dbcf_detect_simd) && defined(DBCF_X86_OR_X64) && !defined(_WIN16) */

#ifndef dbcf_detect_simd
static int dbcF_detect_simd(void)
{
    int ret=0;
#if defined(DBC_FFT_FORCE_SIMD)
    ret=(DBC_FFT_FORCE_SIMD);
#elif defined(DBCF_X86_OR_X64)
    /* Compile-time CPU detection. */
#if defined(__SSE2__) || (_M_IX86_FP>=2) /* MSVC does not have __SSE2__ macro. */
    ret|=DBCF_HAS_SIMD4F |DBCF_HAS_SIMD2D;
#endif
#if !defined(DBC_FFT_NO_AVX) && defined(__AVX__)
    ret|=DBCF_HAS_SIMD8F |DBCF_HAS_SIMD4D;
#endif
#if !defined(DBC_FFT_NO_AVX512) && defined(__AVX512F__)
    ret|=DBCF_HAS_SIMD16F|DBCF_HAS_SIMD8D;
#endif
#if !defined(_WIN16)
    /* Runtime CPU detection. */
    if(dbcF_has_cpuid())
    {
        unsigned eax,ebx,ecx,edx;
        unsigned maxlevel;
        dbcF_cpuid(0,0,&maxlevel,&ebx,&ecx,&edx);
        if(maxlevel>0)
        {
            dbcF_cpuid(1,0,&eax,&ebx,&ecx,&edx);
            if(edx&0x04000000u) ret|=DBCF_HAS_SIMD4F |DBCF_HAS_SIMD2D; /* SSE2 */
#if !defined(DBC_FFT_NO_AVX)
            if(ecx&0x18000000u) /* CPU has AVX & XGETBV. */
            {
                unsigned xcr0=dbcF_xgetbv(0);
                if((xcr0&0x6)==0x6) /* OS-level support for AVX. */
                {
                    ret|=DBCF_HAS_SIMD8F |DBCF_HAS_SIMD4D; /* AVX */
#if !defined(DBC_FFT_NO_AVX512)
                    if(maxlevel>=7)
                    {
                        dbcF_cpuid(7,0,&eax,&ebx,&ecx,&edx);
                        if((xcr0&0xE6)==0xE6) /* OS-level support for AVX512. */
                        {
                            if(ebx&0x00010000u) ret|=DBCF_HAS_SIMD16F|DBCF_HAS_SIMD8D; /* AVX512 */
                        }
                    }
#endif
                }
            }
#endif
        }

    }
#endif /* !defined(_WIN16) */
#endif /* defined(DBC_FFT_FORCE_SIMD) */
    return ret;
}

static int dbcf_detect_simd(void)
{
#if defined(DBC_FFT_CACHE_CPU_DETECTION)
    static int cached=
#ifndef __cplusplus
    -1;
    if(cached==-1) cached=
#endif
    dbcF_detect_simd();
    return cached;
#else
    return dbcF_detect_simd();
#endif /* DBC_FFT_CACHE_CPU_DETECTION */
}
#endif /* dbcf_detect_simd */

#if defined(DBC_FFT_USE_VECTOR_EXTENSIONS)
#if !defined(DBC_FFT_NO_FLOAT) && !defined(DBCF_NO_SIMD4F)
typedef float  dbcf_simd4f  __attribute__((vector_size(16)));
#endif
#if !defined(DBC_FFT_NO_DOUBLE) && !defined(DBCF_NO_SIMD2D)
typedef double dbcf_simd2d  __attribute__((vector_size(16)));
#endif
#if !defined(DBC_FFT_NO_FLOAT) && !defined(DBCF_NO_SIMD8F)
typedef float  dbcf_simd8f  __attribute__((vector_size(32)));
#endif
#if !defined(DBC_FFT_NO_DOUBLE) && !defined(DBCF_NO_SIMD4D)
typedef double dbcf_simd4d  __attribute__((vector_size(32)));
#endif
#if !defined(DBC_FFT_NO_FLOAT) && !defined(DBCF_NO_SIMD16F)
typedef float  dbcf_simd16f __attribute__((vector_size(64)));
#endif
#if !defined(DBC_FFT_NO_DOUBLE) && !defined(DBCF_NO_SIMD8D)
typedef double dbcf_simd8d  __attribute__((vector_size(64)));
#endif
#elif defined(DBCF_X86_OR_X64)
#include <emmintrin.h>
#if !defined(DBC_FFT_NO_FLOAT) && !defined(DBCF_NO_SIMD4F)
typedef __m128  dbcf_simd4f;
#endif
#if !defined(DBC_FFT_NO_DOUBLE) && !defined(DBCF_NO_SIMD2D)
typedef __m128d dbcf_simd2d;
#endif
#ifndef DBC_FFT_NO_AVX
#include <immintrin.h>
#if !defined(DBC_FFT_NO_FLOAT) && !defined(DBCF_NO_SIMD8F)
typedef __m256  dbcf_simd8f;
#endif
#if !defined(DBC_FFT_NO_DOUBLE) && !defined(DBCF_NO_SIMD4D)
typedef __m256d dbcf_simd4d;
#endif
#endif /* DBC_FFT_NO_AVX */
#ifndef DBC_FFT_NO_AVX512
#include <immintrin.h>
#if !defined(DBC_FFT_NO_FLOAT) && !defined(DBCF_NO_SIMD16F)
typedef __m512  dbcf_simd16f;
#endif
#if !defined(DBC_FFT_NO_DOUBLE) && !defined(DBCF_NO_SIMD8D)
typedef __m512d dbcf_simd8d;
#endif
#endif /* DBC_FFT_NO_AVX512 */
#endif /* defined(DBC_FFT_USE_VECTOR_EXTENSIONS) */

/* SIMD functions attributes. */
#ifdef DBCF_X86_OR_X64
#if defined(__GNUC__)
#if defined(__SSE2__)
/*
    SSE2 is already globally enabled. No need to specify it,
    and possibly constrain the compiler from using better
    instructions.
*/
#define DBCF_AUTODECL_SIMD4F
#define DBCF_AUTODECL_SIMD2D
#else
/*
    Note: stdcall is there to work around the problem, pointed by
    "error: calling 'function' with SSE calling convention without SSE/SSE2 enabled"
    diagnostic. Short version: GCC seems to see fit to use SSE
    calling convention on static SSE functions, breaking the calls from
    FPU functions. So does clang, except with no error message.
    Using stdcall overrides that. Possibly unnecessary in this case,
    since the behavior seems to be triggered only if the function
    uses floating point arguments/returns.
    Note: fpmath=sse is not used, as it is not needed (SSE is only used for SIMD),
    and introduces additional hassle.
*/
#define DBCF_AUTODECL_SIMD4F __attribute__((target("sse2"),stdcall))
#define DBCF_AUTODECL_SIMD2D __attribute__((target("sse2"),stdcall))
#endif /* defined(__SSE2__) */
#else
#define DBCF_AUTODECL_SIMD4F
#define DBCF_AUTODECL_SIMD2D
#endif /* defined(__GNUC__) */

#if defined(__GNUC__)
#if defined(__AVX2__)
/*
    AVX2 is already globally enabled. No need to specify it,
    and possibly constrain the compiler from using better
    instructions.
*/
#define DBCF_AUTODECL_SIMD8F
#define DBCF_AUTODECL_SIMD4D
#else
#ifdef DBCF_X64
#define DBCF_AUTODECL_SIMD8F __attribute__((target("avx"))) /* No stdcall in x64. */
#define DBCF_AUTODECL_SIMD4D __attribute__((target("avx"))) /* No stdcall in x64. */
#else
#define DBCF_AUTODECL_SIMD8F __attribute__((target("avx"),stdcall))
#define DBCF_AUTODECL_SIMD4D __attribute__((target("avx"),stdcall))
#endif /* DBCF_X64 */
#endif /* defined(__AVX2__) */
#else
#define DBCF_AUTODECL_SIMD8F
#define DBCF_AUTODECL_SIMD4D
#endif /* defined(__GNUC__) */

#if defined(__GNUC__)
#if defined(__AVX512F__)
/*
    AVX512F is already globally enabled. No need to specify it,
    and possibly constrain the compiler from using better
    instructions.
*/
#define DBCF_AUTODECL_SIMD16F
#define DBCF_AUTODECL_SIMD8D
#else
#ifdef DBCF_X64
#define DBCF_AUTODECL_SIMD16F __attribute__((target("avx512f"))) /* No stdcall in x64. */
#define DBCF_AUTODECL_SIMD8D  __attribute__((target("avx512f"))) /* No stdcall in x64. */
#else
#define DBCF_AUTODECL_SIMD16F __attribute__((target("avx512f"),stdcall))
#define DBCF_AUTODECL_SIMD8D  __attribute__((target("avx512f"),stdcall))
#endif /* DBCF_X64 */
#endif /* defined(__AVX512F__) */
#else
#define DBCF_AUTODECL_SIMD16F
#define DBCF_AUTODECL_SIMD8D
#endif /* defined(__GNUC__) */

#else /* DBCF_X86_OR_X64 */
#define DBCF_AUTODECL_SIMD4F
#define DBCF_AUTODECL_SIMD2D
#define DBCF_AUTODECL_SIMD8F
#define DBCF_AUTODECL_SIMD4D
#define DBCF_AUTODECL_SIMD16F
#define DBCF_AUTODECL_SIMD8D
#endif /* DBCF_X86_OR_X64 */

/* Allow custom declarations. */
#ifndef DBCF_DECL_SIMD4F
#define DBCF_DECL_SIMD4F DBCF_AUTODECL_SIMD4F
#endif
#ifndef DBCF_DECL_SIMD2D
#define DBCF_DECL_SIMD2D DBCF_AUTODECL_SIMD2D
#endif
#ifndef DBCF_DECL_SIMD8F
#define DBCF_DECL_SIMD8F DBCF_AUTODECL_SIMD8F
#endif
#ifndef DBCF_DECL_SIMD4D
#define DBCF_DECL_SIMD4D DBCF_AUTODECL_SIMD4D
#endif
#ifndef DBCF_DECL_SIMD16F
#define DBCF_DECL_SIMD16F DBCF_AUTODECL_SIMD16F
#endif
#ifndef DBCF_DECL_SIMD8D
#define DBCF_DECL_SIMD8D DBCF_AUTODECL_SIMD8D
#endif

#if defined(DBC_FFT_USE_VECTOR_EXTENSIONS)
#define DBCF_DEF_SIMD_LOAD( name,type,alignment) static type name(const void *p) {type ret;__builtin_memcpy(&ret,__builtin_assume_aligned(p,alignment),sizeof(ret));return ret;}
#define DBCF_DEF_SIMD_STORE(name,type,alignment) static void name(type v,void *p) {__builtin_memcpy(__builtin_assume_aligned(p,alignment),&v,sizeof(v));}

#if !defined(DBC_FFT_NO_FLOAT) && !defined(DBCF_NO_SIMD4F)
DBCF_DECL_SIMD4F DBCF_DEF_SIMD_LOAD (dbcF_load4f ,dbcf_simd4f,1)
DBCF_DECL_SIMD4F DBCF_DEF_SIMD_STORE(dbcF_store4f,dbcf_simd4f,1)
DBCF_DECL_SIMD4F DBCF_DEF_SIMD_LOAD (dbcF_load4f_aligned ,dbcf_simd4f,16)
DBCF_DECL_SIMD4F DBCF_DEF_SIMD_STORE(dbcF_store4f_aligned,dbcf_simd4f,16)
DBCF_DECL_SIMD4F static dbcf_simd4f dbcF_set4f(float v0,float v1,float v2,float v3) {return (__extension__ (dbcf_simd4f){v0,v1,v2,v3});}
DBCF_DECL_SIMD4F static dbcf_simd4f dbcF_fill4f(float v) {return dbcF_set4f(v,v,v,v);}
DBCF_DECL_SIMD4F static dbcf_simd4f dbcF_add4f(dbcf_simd4f l,dbcf_simd4f r) {return l+r;}
DBCF_DECL_SIMD4F static dbcf_simd4f dbcF_sub4f(dbcf_simd4f l,dbcf_simd4f r) {return l-r;}
DBCF_DECL_SIMD4F static dbcf_simd4f dbcF_mul4f(dbcf_simd4f l,dbcf_simd4f r) {return l*r;}
#endif
#if !defined(DBC_FFT_NO_DOUBLE) && !defined(DBCF_NO_SIMD2D)
DBCF_DECL_SIMD2D DBCF_DEF_SIMD_LOAD (dbcF_load2d ,dbcf_simd2d,1)
DBCF_DECL_SIMD2D DBCF_DEF_SIMD_STORE(dbcF_store2d,dbcf_simd2d,1)
DBCF_DECL_SIMD2D DBCF_DEF_SIMD_LOAD (dbcF_load2d_aligned ,dbcf_simd2d,16)
DBCF_DECL_SIMD2D DBCF_DEF_SIMD_STORE(dbcF_store2d_aligned,dbcf_simd2d,16)
DBCF_DECL_SIMD2D static dbcf_simd2d dbcF_set2d(double v0,double v1) {return (__extension__ (dbcf_simd2d){v0,v1});}
DBCF_DECL_SIMD2D static dbcf_simd2d dbcF_fill2d(double v) {return dbcF_set2d(v,v);}
DBCF_DECL_SIMD2D static dbcf_simd2d dbcF_add2d(dbcf_simd2d l,dbcf_simd2d r) {return l+r;}
DBCF_DECL_SIMD2D static dbcf_simd2d dbcF_sub2d(dbcf_simd2d l,dbcf_simd2d r) {return l-r;}
DBCF_DECL_SIMD2D static dbcf_simd2d dbcF_mul2d(dbcf_simd2d l,dbcf_simd2d r) {return l*r;}
#endif
#if !defined(DBC_FFT_NO_FLOAT) && !defined(DBCF_NO_SIMD8F)
DBCF_DECL_SIMD8F DBCF_DEF_SIMD_LOAD (dbcF_load8f ,dbcf_simd8f,1)
DBCF_DECL_SIMD8F DBCF_DEF_SIMD_STORE(dbcF_store8f,dbcf_simd8f,1)
DBCF_DECL_SIMD8F DBCF_DEF_SIMD_LOAD (dbcF_load8f_aligned ,dbcf_simd8f,32)
DBCF_DECL_SIMD8F DBCF_DEF_SIMD_STORE(dbcF_store8f_aligned,dbcf_simd8f,32)
DBCF_DECL_SIMD8F static dbcf_simd8f dbcF_set8f(float v0,float v1,float v2,float v3,float v4,float v5,float v6,float v7) {return (__extension__ (dbcf_simd8f){v0,v1,v2,v3,v4,v5,v6,v7});}
DBCF_DECL_SIMD8F static dbcf_simd8f dbcF_fill8f(float v) {return dbcF_set8f(v,v,v,v,v,v,v,v);}
DBCF_DECL_SIMD8F static dbcf_simd8f dbcF_add8f(dbcf_simd8f l,dbcf_simd8f r) {return l+r;}
DBCF_DECL_SIMD8F static dbcf_simd8f dbcF_sub8f(dbcf_simd8f l,dbcf_simd8f r) {return l-r;}
DBCF_DECL_SIMD8F static dbcf_simd8f dbcF_mul8f(dbcf_simd8f l,dbcf_simd8f r) {return l*r;}
#endif
#if !defined(DBC_FFT_NO_DOUBLE) && !defined(DBCF_NO_SIMD4D)
DBCF_DECL_SIMD4D DBCF_DEF_SIMD_LOAD (dbcF_load4d ,dbcf_simd4d,1)
DBCF_DECL_SIMD4D DBCF_DEF_SIMD_STORE(dbcF_store4d,dbcf_simd4d,1)
DBCF_DECL_SIMD4D DBCF_DEF_SIMD_LOAD (dbcF_load4d_aligned ,dbcf_simd4d,32)
DBCF_DECL_SIMD4D DBCF_DEF_SIMD_STORE(dbcF_store4d_aligned,dbcf_simd4d,32)
DBCF_DECL_SIMD4D static dbcf_simd4d dbcF_set4d(double v0,double v1,double v2,double v3) {return (__extension__ (dbcf_simd4d){v0,v1,v2,v3});}
DBCF_DECL_SIMD4D static dbcf_simd4d dbcF_fill4d(double v) {return dbcF_set4d(v,v,v,v);}
DBCF_DECL_SIMD4D static dbcf_simd4d dbcF_add4d(dbcf_simd4d l,dbcf_simd4d r) {return l+r;}
DBCF_DECL_SIMD4D static dbcf_simd4d dbcF_sub4d(dbcf_simd4d l,dbcf_simd4d r) {return l-r;}
DBCF_DECL_SIMD4D static dbcf_simd4d dbcF_mul4d(dbcf_simd4d l,dbcf_simd4d r) {return l*r;}
#endif
#if !defined(DBC_FFT_NO_FLOAT) && !defined(DBCF_NO_SIMD16F)
DBCF_DECL_SIMD16F DBCF_DEF_SIMD_LOAD (dbcF_load16f ,dbcf_simd16f,1)
DBCF_DECL_SIMD16F DBCF_DEF_SIMD_STORE(dbcF_store16f,dbcf_simd16f,1)
DBCF_DECL_SIMD16F DBCF_DEF_SIMD_LOAD (dbcF_load16f_aligned ,dbcf_simd16f,64)
DBCF_DECL_SIMD16F DBCF_DEF_SIMD_STORE(dbcF_store16f_aligned,dbcf_simd16f,64)
DBCF_DECL_SIMD16F static dbcf_simd16f dbcF_set16f(float v0,float v1,float v2,float v3,float v4,float v5,float v6,float v7,float v8,float v9,float v10,float v11,float v12,float v13,float v14,float v15) {return (__extension__ (dbcf_simd16f){v0,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15});}
DBCF_DECL_SIMD16F static dbcf_simd16f dbcF_fill16f(float v) {return dbcF_set16f(v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v);}
DBCF_DECL_SIMD16F static dbcf_simd16f dbcF_add16f(dbcf_simd16f l,dbcf_simd16f r) {return l+r;}
DBCF_DECL_SIMD16F static dbcf_simd16f dbcF_sub16f(dbcf_simd16f l,dbcf_simd16f r) {return l-r;}
DBCF_DECL_SIMD16F static dbcf_simd16f dbcF_mul16f(dbcf_simd16f l,dbcf_simd16f r) {return l*r;}
#endif
#if !defined(DBC_FFT_NO_DOUBLE) && !defined(DBCF_NO_SIMD8D)
DBCF_DECL_SIMD8D DBCF_DEF_SIMD_LOAD (dbcF_load8d ,dbcf_simd8d,1)
DBCF_DECL_SIMD8D DBCF_DEF_SIMD_STORE(dbcF_store8d,dbcf_simd8d,1)
DBCF_DECL_SIMD8D DBCF_DEF_SIMD_LOAD (dbcF_load8d_aligned ,dbcf_simd8d,64)
DBCF_DECL_SIMD8D DBCF_DEF_SIMD_STORE(dbcF_store8d_aligned,dbcf_simd8d,64)
DBCF_DECL_SIMD8D static dbcf_simd8d dbcF_set8d(double v0,double v1,double v2,double v3,double v4,double v5,double v6,double v7) {return (__extension__ (dbcf_simd8d){v0,v1,v2,v3,v4,v5,v6,v7});}
DBCF_DECL_SIMD8D static dbcf_simd8d dbcF_fill8d(double v) {return dbcF_set8d(v,v,v,v,v,v,v,v);}
DBCF_DECL_SIMD8D static dbcf_simd8d dbcF_add8d(dbcf_simd8d l,dbcf_simd8d r) {return l+r;}
DBCF_DECL_SIMD8D static dbcf_simd8d dbcF_sub8d(dbcf_simd8d l,dbcf_simd8d r) {return l-r;}
DBCF_DECL_SIMD8D static dbcf_simd8d dbcF_mul8d(dbcf_simd8d l,dbcf_simd8d r) {return l*r;}
#endif

#elif defined(DBCF_X86_OR_X64)

#if !defined(DBC_FFT_NO_FLOAT) && !defined(DBCF_NO_SIMD4F)
DBCF_DECL_SIMD4F static dbcf_simd4f dbcF_load4f (const void *p) {return _mm_loadu_ps((const float*)p);}
DBCF_DECL_SIMD4F static void        dbcF_store4f(dbcf_simd4f v,void *p) {_mm_storeu_ps((float*)p,v);}
DBCF_DECL_SIMD4F static dbcf_simd4f dbcF_load4f_aligned (const void *p) {return _mm_load_ps((const float*)p);}
DBCF_DECL_SIMD4F static void        dbcF_store4f_aligned(dbcf_simd4f v,void *p) {_mm_store_ps((float*)p,v);}
DBCF_DECL_SIMD4F static dbcf_simd4f dbcF_set4f(float v0,float v1,float v2,float v3) {return _mm_setr_ps(v0,v1,v2,v3);}
DBCF_DECL_SIMD4F static dbcf_simd4f dbcF_fill4f(float v) {return dbcF_set4f(v,v,v,v);}
DBCF_DECL_SIMD4F static dbcf_simd4f dbcF_add4f(dbcf_simd4f l,dbcf_simd4f r) {return _mm_add_ps(l,r);}
DBCF_DECL_SIMD4F static dbcf_simd4f dbcF_sub4f(dbcf_simd4f l,dbcf_simd4f r) {return _mm_sub_ps(l,r);}
DBCF_DECL_SIMD4F static dbcf_simd4f dbcF_mul4f(dbcf_simd4f l,dbcf_simd4f r) {return _mm_mul_ps(l,r);}
#endif
#if !defined(DBC_FFT_NO_DOUBLE) && !defined(DBCF_NO_SIMD2D)
DBCF_DECL_SIMD2D static dbcf_simd2d dbcF_load2d (const void *p) {return _mm_loadu_pd((const double*)p);}
DBCF_DECL_SIMD2D static void        dbcF_store2d(dbcf_simd2d v,void *p) {_mm_storeu_pd((double*)p,v);}
DBCF_DECL_SIMD2D static dbcf_simd2d dbcF_load2d_aligned (const void *p) {return _mm_load_pd((const double*)p);}
DBCF_DECL_SIMD2D static void        dbcF_store2d_aligned(dbcf_simd2d v,void *p) {_mm_store_pd((double*)p,v);}
DBCF_DECL_SIMD2D static dbcf_simd2d dbcF_set2d(double v0,double v1) {return _mm_setr_pd(v0,v1);}
DBCF_DECL_SIMD2D static dbcf_simd2d dbcF_fill2d(double v) {return dbcF_set2d(v,v);}
DBCF_DECL_SIMD2D static dbcf_simd2d dbcF_add2d(dbcf_simd2d l,dbcf_simd2d r) {return _mm_add_pd(l,r);}
DBCF_DECL_SIMD2D static dbcf_simd2d dbcF_sub2d(dbcf_simd2d l,dbcf_simd2d r) {return _mm_sub_pd(l,r);}
DBCF_DECL_SIMD2D static dbcf_simd2d dbcF_mul2d(dbcf_simd2d l,dbcf_simd2d r) {return _mm_mul_pd(l,r);}
#endif
#if !defined(DBC_FFT_NO_FLOAT) && !defined(DBCF_NO_SIMD8F)
DBCF_DECL_SIMD8F static dbcf_simd8f dbcF_load8f (const void *p) {return _mm256_loadu_ps((const float*)p);}
DBCF_DECL_SIMD8F static void        dbcF_store8f(dbcf_simd8f v,void *p) {_mm256_storeu_ps((float*)p,v);}
DBCF_DECL_SIMD8F static dbcf_simd8f dbcF_load8f_aligned (const void *p) {return _mm256_load_ps((const float*)p);}
DBCF_DECL_SIMD8F static void        dbcF_store8f_aligned(dbcf_simd8f v,void *p) {_mm256_store_ps((float*)p,v);}
DBCF_DECL_SIMD8F static dbcf_simd8f dbcF_set8f(float v0,float v1,float v2,float v3,float v4,float v5,float v6,float v7) {return _mm256_setr_ps(v0,v1,v2,v3,v4,v5,v6,v7);}
DBCF_DECL_SIMD8F static dbcf_simd8f dbcF_fill8f(float v) {return dbcF_set8f(v,v,v,v,v,v,v,v);}
DBCF_DECL_SIMD8F static dbcf_simd8f dbcF_add8f(dbcf_simd8f l,dbcf_simd8f r) {return _mm256_add_ps(l,r);}
DBCF_DECL_SIMD8F static dbcf_simd8f dbcF_sub8f(dbcf_simd8f l,dbcf_simd8f r) {return _mm256_sub_ps(l,r);}
DBCF_DECL_SIMD8F static dbcf_simd8f dbcF_mul8f(dbcf_simd8f l,dbcf_simd8f r) {return _mm256_mul_ps(l,r);}
#endif
#if !defined(DBC_FFT_NO_DOUBLE) && !defined(DBCF_NO_SIMD4D)
DBCF_DECL_SIMD4D static dbcf_simd4d dbcF_load4d (const void *p) {return _mm256_loadu_pd((const double*)p);}
DBCF_DECL_SIMD4D static void        dbcF_store4d(dbcf_simd4d v,void *p) {_mm256_storeu_pd((double*)p,v);}
DBCF_DECL_SIMD4D static dbcf_simd4d dbcF_load4d_aligned (const void *p) {return _mm256_load_pd((const double*)p);}
DBCF_DECL_SIMD4D static void        dbcF_store4d_aligned(dbcf_simd4d v,void *p) {_mm256_store_pd((double*)p,v);}
DBCF_DECL_SIMD4D static dbcf_simd4d dbcF_set4d(double v0,double v1,double v2,double v3) {return _mm256_setr_pd(v0,v1,v2,v3);}
DBCF_DECL_SIMD4D static dbcf_simd4d dbcF_fill4d(double v) {return dbcF_set4d(v,v,v,v);}
DBCF_DECL_SIMD4D static dbcf_simd4d dbcF_add4d(dbcf_simd4d l,dbcf_simd4d r) {return _mm256_add_pd(l,r);}
DBCF_DECL_SIMD4D static dbcf_simd4d dbcF_sub4d(dbcf_simd4d l,dbcf_simd4d r) {return _mm256_sub_pd(l,r);}
DBCF_DECL_SIMD4D static dbcf_simd4d dbcF_mul4d(dbcf_simd4d l,dbcf_simd4d r) {return _mm256_mul_pd(l,r);}
#endif
#if !defined(DBC_FFT_NO_FLOAT) && !defined(DBCF_NO_SIMD16F)
DBCF_DECL_SIMD16F static dbcf_simd16f dbcF_load16f (const void *p) {return _mm512_loadu_ps((const float*)p);}
DBCF_DECL_SIMD16F static void        dbcF_store16f(dbcf_simd16f v,void *p) {_mm512_storeu_ps((float*)p,v);}
DBCF_DECL_SIMD16F static dbcf_simd16f dbcF_load16f_aligned (const void *p) {return _mm512_load_ps((const float*)p);}
DBCF_DECL_SIMD16F static void        dbcF_store16f_aligned(dbcf_simd16f v,void *p) {_mm512_store_ps((float*)p,v);}
DBCF_DECL_SIMD16F static dbcf_simd16f dbcF_set16f(float v0,float v1,float v2,float v3,float v4,float v5,float v6,float v7,float v8,float v9,float v10,float v11,float v12,float v13,float v14,float v15) {return _mm512_setr_ps(v0,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15);}
DBCF_DECL_SIMD16F static dbcf_simd16f dbcF_fill16f(float v) {return dbcF_set16f(v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v);}
DBCF_DECL_SIMD16F static dbcf_simd16f dbcF_add16f(dbcf_simd16f l,dbcf_simd16f r) {return _mm512_add_ps(l,r);}
DBCF_DECL_SIMD16F static dbcf_simd16f dbcF_sub16f(dbcf_simd16f l,dbcf_simd16f r) {return _mm512_sub_ps(l,r);}
DBCF_DECL_SIMD16F static dbcf_simd16f dbcF_mul16f(dbcf_simd16f l,dbcf_simd16f r) {return _mm512_mul_ps(l,r);}
#endif
#if !defined(DBC_FFT_NO_DOUBLE) && !defined(DBCF_NO_SIMD8D)
DBCF_DECL_SIMD8D static dbcf_simd8d dbcF_load8d (const void *p) {return _mm512_loadu_pd((const double*)p);}
DBCF_DECL_SIMD8D static void        dbcF_store8d(dbcf_simd8d v,void *p) {_mm512_storeu_pd((double*)p,v);}
DBCF_DECL_SIMD8D static dbcf_simd8d dbcF_load8d_aligned (const void *p) {return _mm512_load_pd((const double*)p);}
DBCF_DECL_SIMD8D static void        dbcF_store8d_aligned(dbcf_simd8d v,void *p) {_mm512_store_pd((double*)p,v);}
DBCF_DECL_SIMD8D static dbcf_simd8d dbcF_set8d(double v0,double v1,double v2,double v3,double v4,double v5,double v6,double v7) {return _mm512_setr_pd(v0,v1,v2,v3,v4,v5,v6,v7);}
DBCF_DECL_SIMD8D static dbcf_simd8d dbcF_fill8d(double v) {return dbcF_set8d(v,v,v,v,v,v,v,v);}
DBCF_DECL_SIMD8D static dbcf_simd8d dbcF_add8d(dbcf_simd8d l,dbcf_simd8d r) {return _mm512_add_pd(l,r);}
DBCF_DECL_SIMD8D static dbcf_simd8d dbcF_sub8d(dbcf_simd8d l,dbcf_simd8d r) {return _mm512_sub_pd(l,r);}
DBCF_DECL_SIMD8D static dbcf_simd8d dbcF_mul8d(dbcf_simd8d l,dbcf_simd8d r) {return _mm512_mul_pd(l,r);}
#endif

#endif /* defined(DBC_FFT_USE_VECTOR_EXTENSIONS) */

#define DBCF_DEF_SIMD_BLOCK(name,type,size,simd,load_t,load_d,store_d,fill,add,sub,mul,cexp)\
static void name(                                                                             \
    dbcf_index log2n,                                                                         \
    dbcf_index log2b,                                                                         \
    type *LR,type *LI,                                                                        \
    type *HR,type *HI,                                                                        \
    type C,type S,                                                                            \
    int inverse,                                                                              \
    const type *tr,const type *ti)                                                            \
{                                                                                             \
    dbcf_index b=DBCF_POW2(log2b),h=b>>1;                                                     \
    if(log2b<=DBCF_TWIDDLES_BUF_LOG2)                                                         \
    {                                                                                         \
        /* The block is small, we have enough precomputed twiddles. */                        \
        dbcf_index i;                                                                         \
        simd CC=fill(C),SS=fill(S);                                                           \
        for(i=0;i<b;i+=size)                                                                  \
        {                                                                                     \
            simd TR=load_t(tr+i),TI=load_t(ti+i);                                             \
            simd c=sub(mul(CC,TR),mul(SS,TI)),s=add(mul(SS,TR),mul(CC,TI));                   \
            simd xl=load_d(LR+i),yl=load_d(LI+i);                                             \
            simd xr=load_d(HR+i),yr=load_d(HI+i);                                             \
            simd x=sub(mul(c,xr),mul(s,yr)),y=add(mul(s,xr),mul(c,yr));                       \
            store_d(add(xl,x),LR+i);store_d(add(yl,y),LI+i);                                  \
            store_d(sub(xl,x),HR+i);store_d(sub(yl,y),HI+i);                                  \
        }                                                                                     \
    }                                                                                         \
    else                                                                                      \
    {                                                                                         \
        type X,Y;                                                                             \
        /* The block is large, we process it's halves recursively. */                         \
        cexp(log2n-log2b+1,&X,&Y);                                                            \
        if(!inverse) Y=-Y;                                                                    \
        name(log2n,log2b-1,LR  ,LI  ,HR  ,HI  ,C      ,S      ,inverse,tr,ti);                \
        name(log2n,log2b-1,LR+h,LI+h,HR+h,HI+h,C*X-S*Y,S*X+C*Y,inverse,tr,ti);                \
    }                                                                                         \
}

#define DBCF_DEF_SIMD_PASS(name,type,size,simd,load_t,load_d,store_d,add,sub,mul,block)\
static void name(                                                                             \
    dbcf_index log2n,                                                                         \
    dbcf_index log2c,                                                                         \
    type *real,type *imag,                                                                    \
    int inverse,                                                                              \
	dbcf_index log2t,                                                                         \
    const type *tr,const type *ti)                                                            \
{                                                                                             \
    dbcf_index n=DBCF_POW2(log2n),h=n>>1;                                                     \
    dbcf_index c=DBCF_POW2(log2c);                                                            \
    dbcf_index i;                                                                             \
    type *LR=real;                                                                            \
    type *HR=real+h;                                                                          \
    type *LI=imag;                                                                            \
    type *HI=imag+h;                                                                          \
                                                                                              \
    if(log2n-1<=log2t)                                                                        \
    {                                                                                         \
        /*  We have as much precomputed twiddles            */                                \
        /*  as the block needs, so we supply them directly. */                                \
        if(h>size) for(i=0;i<c;++i)                                                           \
        {                                                                                     \
            /* Unroll x2. Slightly faster. */                                                 \
            dbcf_index d;                                                                     \
            for(d=0;d<h;)                                                                     \
            {                                                                                 \
                simd c=load_t(tr+d),s=load_t(ti+d);                                           \
                simd xl=load_d(LR+d),yl=load_d(LI+d);                                         \
                simd xr=load_d(HR+d),yr=load_d(HI+d);                                         \
                simd x=sub(mul(c,xr),mul(s,yr)),y=add(mul(s,xr),mul(c,yr));                   \
                store_d(add(xl,x),LR+d);store_d(add(yl,y),LI+d);                              \
                store_d(sub(xl,x),HR+d);store_d(sub(yl,y),HI+d);                              \
                d+=size;                                                                      \
                c=load_t(tr+d);s=load_t(ti+d);                                                \
                xl=load_d(LR+d);yl=load_d(LI+d);                                              \
                xr=load_d(HR+d);yr=load_d(HI+d);                                              \
                x=sub(mul(c,xr),mul(s,yr));y=add(mul(s,xr),mul(c,yr));                        \
                store_d(add(xl,x),LR+d);store_d(add(yl,y),LI+d);                              \
                store_d(sub(xl,x),HR+d);store_d(sub(yl,y),HI+d);                              \
                d+=size;                                                                      \
            }                                                                                 \
            LR+=n;LI+=n;                                                                      \
            HR+=n;HI+=n;                                                                      \
        }                                                                                     \
        else for(i=0;i<c;++i)                                                                 \
        {                                                                                     \
            simd C=load_t(tr),S=load_t(ti);                                                   \
            simd xl=load_d(LR),yl=load_d(LI);                                                 \
            simd xr=load_d(HR),yr=load_d(HI);                                                 \
            simd x=sub(mul(C,xr),mul(S,yr)),y=add(mul(S,xr),mul(C,yr));                       \
            store_d(add(xl,x),LR);store_d(add(yl,y),LI);                                      \
            store_d(sub(xl,x),HR);store_d(sub(yl,y),HI);                                      \
            LR+=n;LI+=n;                                                                      \
            HR+=n;HI+=n;                                                                      \
        }                                                                                     \
    }                                                                                         \
    else                                                                                      \
    {                                                                                         \
        for(i=0;i<c;++i)                                                                      \
        {                                                                                     \
            block(log2n,log2n-1,LR,LI,HR,HI,1.0f,0.0f,inverse,tr,ti);                         \
            LR+=n;LI+=n;                                                                      \
            HR+=n;HI+=n;                                                                      \
        }                                                                                     \
    }                                                                                         \
}

/* Not actually SIMDified, but at least uses the right instruction level. */
#define DBCF_DEF_SIMD_FFT8(name,type,suffix,lsuffix)\
static void name(type *real,type *imag,int inverse) {dbcF_fft8_##suffix(real,imag,1,1,inverse,0.70710678118654752438##lsuffix);}

#define DBCF_DEF_SIMD_COMPUTE_TWIDDLES(name,type,size,simd,lsuffix,load,store,fill,add,sub,mul,cexpm1)\
static void name(dbcf_index log2n,dbcf_index log2b,type *real,type *imag,int inverse)         \
{                                                                                             \
    dbcf_index i;                                                                             \
    real[0]=0.0##lsuffix;                                                                     \
    imag[0]=0.0##lsuffix;                                                                     \
    for(i=0;i<log2b;++i)                                                                      \
    {                                                                                         \
        dbcf_index j,k=DBCF_POW2(i);                                                          \
        type x,y;                                                                             \
        cexpm1(log2n-i,&x,&y);                                                                \
        if(!inverse) y=-y;                                                                    \
        if(k>=size)                                                                           \
        {                                                                                     \
            simd X=fill(x);                                                                   \
            simd Y=fill(y);                                                                   \
            for(j=0;j<k;j+=size)                                                              \
            {                                                                                 \
                simd R=load(real+j);                                                          \
                simd I=load(imag+j);                                                          \
                store(add(sub(mul(X,R),mul(Y,I)),add(X,R)),real+k+j);                         \
                store(add(add(mul(Y,R),mul(X,I)),add(Y,I)),imag+k+j);                         \
            }                                                                                 \
        }                                                                                     \
        else for(j=0;j<k;++j)                                                                 \
        {                                                                                     \
            real[k+j]=(x*real[j]-y*imag[j])+(x+real[j]);                                      \
            imag[k+j]=(y*real[j]+x*imag[j])+(y+imag[j]);                                      \
        }                                                                                     \
    }                                                                                         \
    for(i=0;i<DBCF_POW2(log2b);++i)                                                           \
        real[i]=1.0##lsuffix+real[i];                                                         \
}

#define DBCF_DEF_SIMD_FUNCTIONS(decl,type,size,suffix,lsuffix)\
    decl DBCF_DEF_SIMD_BLOCK(dbcF_butterfly_block_##size##suffix##_uu,type,size,dbcf_simd##size##suffix,dbcF_load##size##suffix          ,dbcF_load##size##suffix          ,dbcF_store##size##suffix          ,dbcF_fill##size##suffix,dbcF_add##size##suffix,dbcF_sub##size##suffix,dbcF_mul##size##suffix,dbcF_cexp_##suffix)\
    decl DBCF_DEF_SIMD_BLOCK(dbcF_butterfly_block_##size##suffix##_au,type,size,dbcf_simd##size##suffix,dbcF_load##size##suffix##_aligned,dbcF_load##size##suffix          ,dbcF_store##size##suffix          ,dbcF_fill##size##suffix,dbcF_add##size##suffix,dbcF_sub##size##suffix,dbcF_mul##size##suffix,dbcF_cexp_##suffix)\
    decl DBCF_DEF_SIMD_BLOCK(dbcF_butterfly_block_##size##suffix##_ua,type,size,dbcf_simd##size##suffix,dbcF_load##size##suffix          ,dbcF_load##size##suffix##_aligned,dbcF_store##size##suffix##_aligned,dbcF_fill##size##suffix,dbcF_add##size##suffix,dbcF_sub##size##suffix,dbcF_mul##size##suffix,dbcF_cexp_##suffix)\
    decl DBCF_DEF_SIMD_BLOCK(dbcF_butterfly_block_##size##suffix##_aa,type,size,dbcf_simd##size##suffix,dbcF_load##size##suffix##_aligned,dbcF_load##size##suffix##_aligned,dbcF_store##size##suffix##_aligned,dbcF_fill##size##suffix,dbcF_add##size##suffix,dbcF_sub##size##suffix,dbcF_mul##size##suffix,dbcF_cexp_##suffix)\
    decl DBCF_DEF_SIMD_PASS(dbcF_butterfly_pass_##size##suffix##_uu,type,size,dbcf_simd##size##suffix,dbcF_load##size##suffix          ,dbcF_load##size##suffix          ,dbcF_store##size##suffix          ,dbcF_add##size##suffix,dbcF_sub##size##suffix,dbcF_mul##size##suffix,dbcF_butterfly_block_##size##suffix##_uu)\
    decl DBCF_DEF_SIMD_PASS(dbcF_butterfly_pass_##size##suffix##_au,type,size,dbcf_simd##size##suffix,dbcF_load##size##suffix##_aligned,dbcF_load##size##suffix          ,dbcF_store##size##suffix          ,dbcF_add##size##suffix,dbcF_sub##size##suffix,dbcF_mul##size##suffix,dbcF_butterfly_block_##size##suffix##_au)\
    decl DBCF_DEF_SIMD_PASS(dbcF_butterfly_pass_##size##suffix##_ua,type,size,dbcf_simd##size##suffix,dbcF_load##size##suffix          ,dbcF_load##size##suffix##_aligned,dbcF_store##size##suffix##_aligned,dbcF_add##size##suffix,dbcF_sub##size##suffix,dbcF_mul##size##suffix,dbcF_butterfly_block_##size##suffix##_ua)\
    decl DBCF_DEF_SIMD_PASS(dbcF_butterfly_pass_##size##suffix##_aa,type,size,dbcf_simd##size##suffix,dbcF_load##size##suffix##_aligned,dbcF_load##size##suffix##_aligned,dbcF_store##size##suffix##_aligned,dbcF_add##size##suffix,dbcF_sub##size##suffix,dbcF_mul##size##suffix,dbcF_butterfly_block_##size##suffix##_aa)\
    decl DBCF_DEF_SIMD_COMPUTE_TWIDDLES(dbcF_compute_twiddles_##size##suffix##_u,type,size,dbcf_simd##size##suffix,lsuffix,dbcF_load##size##suffix          ,dbcF_store##size##suffix          ,dbcF_fill##size##suffix,dbcF_add##size##suffix,dbcF_sub##size##suffix,dbcF_mul##size##suffix,dbcF_cexpm1_##suffix)\
    decl DBCF_DEF_SIMD_COMPUTE_TWIDDLES(dbcF_compute_twiddles_##size##suffix##_a,type,size,dbcf_simd##size##suffix,lsuffix,dbcF_load##size##suffix##_aligned,dbcF_store##size##suffix##_aligned,dbcF_fill##size##suffix,dbcF_add##size##suffix,dbcF_sub##size##suffix,dbcF_mul##size##suffix,dbcF_cexpm1_##suffix)\
    decl DBCF_DEF_SIMD_FFT8(dbcF_fft8_##size##suffix,type,suffix,lsuffix)

#ifndef DBC_FFT_NO_FLOAT
static void dbcF_cexpm1_f(dbcf_index log2n,float  *real,float  *imag);
static void dbcF_cexp_f(dbcf_index log2n,float  *real,float  *imag);
static void dbcF_fft8_f(float *real,float *imag,dbcf_index real_stride,dbcf_index imag_stride,int inverse,float c);
#endif
#ifndef DBC_FFT_NO_DOUBLE
static void dbcF_cexpm1_d(dbcf_index log2n,double *real,double *imag);
static void dbcF_cexp_d(dbcf_index log2n,double *real,double *imag);
static void dbcF_fft8_d(double *real,double *imag,dbcf_index real_stride,dbcf_index imag_stride,int inverse,double c);
#endif

#if !defined(DBC_FFT_NO_FLOAT) && !defined(DBCF_NO_SIMD4F)
DBCF_DEF_SIMD_FUNCTIONS(DBCF_DECL_SIMD4F ,float , 4,f,f)
#endif
#if !defined(DBC_FFT_NO_DOUBLE) && !defined(DBCF_NO_SIMD2D)
DBCF_DEF_SIMD_FUNCTIONS(DBCF_DECL_SIMD2D ,double, 2,d,e0)
#endif
#if !defined(DBC_FFT_NO_FLOAT) && !defined(DBCF_NO_SIMD8F)
DBCF_DEF_SIMD_FUNCTIONS(DBCF_DECL_SIMD8F ,float , 8,f,f)
#endif
#if !defined(DBC_FFT_NO_DOUBLE) && !defined(DBCF_NO_SIMD4D)
DBCF_DEF_SIMD_FUNCTIONS(DBCF_DECL_SIMD4D ,double, 4,d,e0)
#endif
#if !defined(DBC_FFT_NO_FLOAT) && !defined(DBCF_NO_SIMD16F)
DBCF_DEF_SIMD_FUNCTIONS(DBCF_DECL_SIMD16F,float ,16,f,f)
#endif
#if !defined(DBC_FFT_NO_DOUBLE) && !defined(DBCF_NO_SIMD8D)
DBCF_DEF_SIMD_FUNCTIONS(DBCF_DECL_SIMD8D ,double, 8,d,e0)
#endif

#define DBCF_TRY_SIMD_PASS(type,size,suffix,SUFFIX)\
    if(simd_flags&DBCF_HAS_SIMD##size##SUFFIX)                                                                        \
    {                                                                                                                 \
        int alignt=DBCF_IS_ALIGNED(tr  ,size*sizeof(type))&&DBCF_IS_ALIGNED(ti  ,size*sizeof(type));                  \
        int alignd=DBCF_IS_ALIGNED(real,size*sizeof(type))&&DBCF_IS_ALIGNED(imag,size*sizeof(type));                  \
        if(((size<<2)>>log2n)<=1&&((size<<1)>>log2t)<=1)                                                              \
        {                                                                                                             \
            if(alignt) dbcF_compute_twiddles_##size##suffix##_a(log2n,log2t,tr,ti,inverse);                           \
            else       dbcF_compute_twiddles_##size##suffix##_u(log2n,log2t,tr,ti,inverse);                           \
            switch(2*alignd+alignt)                                                                                   \
            {                                                                                                         \
                case 0: dbcF_butterfly_pass_##size##suffix##_uu(log2n,log2c,real,imag,inverse,log2t,tr,ti); break;    \
                case 1: dbcF_butterfly_pass_##size##suffix##_au(log2n,log2c,real,imag,inverse,log2t,tr,ti); break;    \
                case 2: dbcF_butterfly_pass_##size##suffix##_ua(log2n,log2c,real,imag,inverse,log2t,tr,ti); break;    \
                case 3: dbcF_butterfly_pass_##size##suffix##_aa(log2n,log2c,real,imag,inverse,log2t,tr,ti); break;    \
            }                                                                                                         \
            return 1;                                                                                                 \
        }                                                                                                             \
    }

#ifndef DBC_FFT_NO_FLOAT
static dbcf_index dbcF_butterfly_pass_optimized_float(
    dbcf_index log2n,
    dbcf_index log2c,
    float *real,float *imag,
    int inverse,
	dbcf_index log2t,
    float *tr,float *ti,
    int simd_flags)
{
#ifndef DBCF_NO_SIMD16F
    DBCF_TRY_SIMD_PASS(float,16,f,F);
#endif
#ifndef DBCF_NO_SIMD8F
    DBCF_TRY_SIMD_PASS(float,8,f,F);
#endif
#ifndef DBCF_NO_SIMD4F
    DBCF_TRY_SIMD_PASS(float,4,f,F);
#endif
    return 0;
}

/* Returns the number of passes actually performed (always contiguous, starting from log2n-depth+1). */
static dbcf_index dbcF_butterfly_multipass_optimized_float(
    dbcf_index log2n,
    dbcf_index log2c,
    dbcf_index depth,
    float *real,float *imag,
    dbcf_index real_stride,dbcf_index imag_stride,
    int inverse,
    float *tr,float *ti)
{
    dbcf_index ret=0;
    int simd_flags;
    if(real_stride!=1||imag_stride!=1) return 0;
    simd_flags=dbcf_detect_simd();
    if(!(simd_flags&(DBCF_HAS_SIMD4F|DBCF_HAS_SIMD8F|DBCF_HAS_SIMD16F))) return 0;
    if(DBCF_TWIDDLES_BUF_LOG2<3) return 0;
    if(depth==log2n&&depth>=3)
    {
        dbcf_index j,m=DBCF_POW2(log2n+log2c-3);
#ifndef DBCF_NO_SIMD16F
        if(simd_flags&DBCF_HAS_SIMD16F) {for(j=0;j<m;++j) {dbcF_fft8_16f(real+8*j,imag+8*j,inverse);} goto ok;}
#endif
#ifndef DBCF_NO_SIMD8F
        if(simd_flags&DBCF_HAS_SIMD8F)  {for(j=0;j<m;++j) {dbcF_fft8_8f (real+8*j,imag+8*j,inverse);} goto ok;}
#endif
#ifndef DBCF_NO_SIMD4F
        if(simd_flags&DBCF_HAS_SIMD4F)  {for(j=0;j<m;++j) {dbcF_fft8_4f (real+8*j,imag+8*j,inverse);} goto ok;}
#endif
        return 0;
        ok: depth-=3; ret=3;
    }
    if(log2n-depth+1>3)
    {
        dbcf_index log2d;
        for(log2d=log2n-depth+1;log2d<=log2n;++log2d)
        {
            dbcf_index log2t=(log2d-1<DBCF_TWIDDLES_BUF_LOG2?log2d-1:DBCF_TWIDDLES_BUF_LOG2);
            if(dbcF_butterfly_pass_optimized_float(log2d,log2c+log2n-log2d,real,imag,inverse,log2t,tr,ti,simd_flags)) ++ret;
            else break;
        }
        return ret;
    }
    return 0;
}
#endif /* DBC_FFT_NO_FLOAT */

#ifndef DBC_FFT_NO_DOUBLE
static dbcf_index dbcF_butterfly_pass_optimized_double(
    dbcf_index log2n,
    dbcf_index log2c,
    double *real,double *imag,
    int inverse,
	dbcf_index log2t,
    double *tr,double *ti,
    int simd_flags)
{
#ifndef DBCF_NO_SIMD8D
    DBCF_TRY_SIMD_PASS(double,8,d,D);
#endif
#ifndef DBCF_NO_SIMD4D
    DBCF_TRY_SIMD_PASS(double,4,d,D);
#endif
#ifndef DBCF_NO_SIMD2D
    DBCF_TRY_SIMD_PASS(double,2,d,D);
#endif
    return 0;
}

/* Returns the number of passes actually performed (always contiguous, starting from log2n-depth+1). */
static dbcf_index dbcF_butterfly_multipass_optimized_double(
    dbcf_index log2n,
    dbcf_index log2c,
    dbcf_index depth,
    double *real,double *imag,
    dbcf_index real_stride,dbcf_index imag_stride,
    int inverse,
    double *tr,double *ti)
{
    dbcf_index ret=0;
    int simd_flags;
    if(real_stride!=1||imag_stride!=1) return 0;
    simd_flags=dbcf_detect_simd();
    if(!(simd_flags&(DBCF_HAS_SIMD2D|DBCF_HAS_SIMD4D))) return 0;
    if(DBCF_TWIDDLES_BUF_LOG2<2) return 0;
    if(depth==log2n&&depth>=3)
    {
        dbcf_index j,m=DBCF_POW2(log2n+log2c-3);
#ifndef DBCF_NO_SIMD8D
        if(simd_flags&DBCF_HAS_SIMD8D) {for(j=0;j<m;++j) {dbcF_fft8_8d(real+8*j,imag+8*j,inverse);} goto ok;}
#endif
#ifndef DBCF_NO_SIMD4D
        if(simd_flags&DBCF_HAS_SIMD4D) {for(j=0;j<m;++j) {dbcF_fft8_4d(real+8*j,imag+8*j,inverse);} goto ok;}
#endif
#ifndef DBCF_NO_SIMD2D
        if(simd_flags&DBCF_HAS_SIMD2D) {for(j=0;j<m;++j) {dbcF_fft8_2d(real+8*j,imag+8*j,inverse);} goto ok;}
#endif
        ok: depth-=3; ret=3;
    }
    if(log2n-depth+1>3)
    {
        dbcf_index log2d;
        for(log2d=log2n-depth+1;log2d<=log2n;++log2d)
        {
            dbcf_index log2t=(log2d-1<DBCF_TWIDDLES_BUF_LOG2?log2d-1:DBCF_TWIDDLES_BUF_LOG2);
            if(dbcF_butterfly_pass_optimized_double(log2d,log2c+log2n-log2d,real,imag,inverse,log2t,tr,ti,simd_flags)) ++ret;
            else break;
        }
        return ret;
    }
    return 0;
}
#endif /* DBC_FFT_NO_DOUBLE */

#endif /* DBC_FFT_NO_SIMD */

static void dbcF_init()
{
#if !defined(DBC_FFT_NO_SIMD) && !defined(dbcf_detect_simd) && defined(DBC_FFT_CACHE_CPU_DETECTION)
    (void)dbcf_detect_simd();
#endif
}

/* Instantiations */
#define DBC_FFT_INSTANTIATION

#ifndef DBC_FFT_NO_FLOAT

#define DBCF_Type float
#define DBCF_Id f
#define DBCF_LITERAL(x) DBCF_CONCAT(x,f)
#ifndef DBC_FFT_NO_SIMD
#define DBCF_butterfly_multipass_optimized dbcF_butterfly_multipass_optimized_float
#endif
#include __FILE__
#undef DBCF_Type
#undef DBCF_Id
#undef DBCF_LITERAL
#ifndef DBC_FFT_NO_SIMD
#undef DBCF_butterfly_multipass_optimized
#endif

#endif /* DBC_FFT_NO_FLOAT */

#ifndef DBC_FFT_NO_DOUBLE

#define DBCF_Type double
#define DBCF_Id d
#define DBCF_LITERAL(x) (x)
#ifndef DBC_FFT_NO_SIMD
#define DBCF_butterfly_multipass_optimized dbcF_butterfly_multipass_optimized_double
#endif
#include __FILE__
#undef DBCF_Type
#undef DBCF_Id
#undef DBCF_LITERAL
#ifndef DBC_FFT_NO_SIMD
#undef DBCF_butterfly_multipass_optimized
#endif

#endif /* DBC_FFT_NO_DOUBLE */

#ifndef DBC_FFT_NO_LONGDOUBLE

#define DBCF_Type long double
#define DBCF_Id l
#define DBCF_LITERAL(x) DBCF_CONCAT(x,l)
#include __FILE__
#undef DBCF_Type
#undef DBCF_Id
#undef DBCF_LITERAL

#endif /* DBC_FFT_NO_LONGDOUBLE */

#undef DBC_FFT_INSTANTIATION

#endif /* defined(DBC_FFT_IMPLEMENTATION) && !defined(DBC_FFT_DECLARATION) && !defined(DBC_FFT_INSTANTIATION) */

/*============================================================================*/
/* Instantiation section. */

#ifdef DBC_FFT_INSTANTIATION

#ifndef DBCF_ZERO
#define DBCF_ZERO DBCF_LITERAL(0.0)
#endif
#ifndef DBCF_ONE
#define DBCF_ONE  DBCF_LITERAL(1.0)
#endif

/*
    Compute exp(2*pi*i/n)-1.
    We don't actually need to depend on <math.h>.
*/
static void DBCF_NAME(dbcF_cexpm1)(dbcf_index log2n,DBCF_Type *real,DBCF_Type *imag)
{
#ifdef DBCF_cexpm1
    DBCF_cexpm1(log2n,real,imag);
#else
    /* All coefficients are entered with sub-ULP precision, even for float128. */
    static DBCF_Type table[][2]={
        { DBCF_LITERAL(0.0e0)                                       ,DBCF_LITERAL(0.0e0)                                       },
        {-DBCF_LITERAL(2.0e0)                                       ,DBCF_LITERAL(0.0e0)                                       },
        {-DBCF_LITERAL(1.0e0)                                       ,DBCF_LITERAL(1.0e0)                                       },
        {-DBCF_LITERAL(2.928932188134524755991556378951509607151e-1),DBCF_LITERAL(7.071067811865475244008443621048490392848e-1)},
        {-DBCF_LITERAL(7.612046748871324387181681060321171317758e-2),DBCF_LITERAL(3.826834323650897717284599840303988667613e-1)},
        {-DBCF_LITERAL(1.921471959676955087381776386576096302606e-2),DBCF_LITERAL(1.950903220161282678482848684770222409276e-1)},
        {-DBCF_LITERAL(4.815273327803113755163046890520078424525e-3),DBCF_LITERAL(9.801714032956060199419556388864184586113e-2)},
        {-DBCF_LITERAL(1.204543794827607285228395240899305556796e-3),DBCF_LITERAL(4.906767432741801425495497694268265831474e-2)},
        {-DBCF_LITERAL(3.011813037957798842343503338278031499389e-4),DBCF_LITERAL(2.454122852291228803173452945928292506546e-2)},
        {-DBCF_LITERAL(7.529816085545907835350880361677564939353e-5),DBCF_LITERAL(1.227153828571992607940826195100321214037e-2)},
        {-DBCF_LITERAL(1.882471739885734300956227143228382608274e-5),DBCF_LITERAL(6.135884649154475359640234590372580917057e-3)},
        {-DBCF_LITERAL(4.706190423828488419874299880100447012366e-6),DBCF_LITERAL(3.067956762965976270145365490919842518944e-3)},
        {-DBCF_LITERAL(1.176548298090070974289828473980951732077e-6),DBCF_LITERAL(1.533980186284765612303697150264079079954e-3)},
        {-DBCF_LITERAL(2.941371177808397717822612343228837361006e-7),DBCF_LITERAL(7.669903187427045269385683579485766431409e-4)},
        {-DBCF_LITERAL(7.353428214885526851929261214305179884431e-8),DBCF_LITERAL(3.834951875713955890724616811813812633950e-4)},
        {-DBCF_LITERAL(1.838357070619165308459709028549492394875e-8),DBCF_LITERAL(1.917475973107033074399095619890009334688e-4)},
        {-DBCF_LITERAL(4.595892687109028066860393851041105696810e-9),DBCF_LITERAL(9.587379909597734587051721097647635118706e-5)}
    };
    if(log2n<(dbcf_index)(sizeof(table)/(sizeof(table[0]))))
    {
        *real=table[log2n][0];
        *imag=table[log2n][1];
    }
    else
    {
        /* For small x Taylor series is accurate to couple of ULPs, even for float128. */
        dbcf_index n=((dbcf_index)1)<<log2n;
        const DBCF_Type C1=DBCF_LITERAL(1.0e0);
        const DBCF_Type C2=DBCF_LITERAL(5.0e-1);
        const DBCF_Type C3=DBCF_LITERAL(1.666666666666666666666666666666666666666e-1);
        const DBCF_Type C4=DBCF_LITERAL(4.166666666666666666666666666666666666666e-2);
        const DBCF_Type C5=DBCF_LITERAL(8.333333333333333333333333333333333333333e-3);
        const DBCF_Type C6=DBCF_LITERAL(1.388888888888888888888888888888888888888e-3);
        const DBCF_Type C7=DBCF_LITERAL(1.984126984126984126984126984126984126984e-4);
        const DBCF_Type C8=DBCF_LITERAL(2.480158730158730158730158730158730158730e-5);
        DBCF_Type x=DBCF_LITERAL(6.283185307179586476925286766559005768)/(DBCF_Type)n;
        DBCF_Type x2=x*x;
        *real=-x2*(C2-x2*(C4-x2*(C6-x2*C8)));
        *imag=x*(C1-x2*(C3-x2*(C5-x2*C7)));
    }
#endif
}

static void DBCF_NAME(dbcF_cexp)(dbcf_index log2n,DBCF_Type *real,DBCF_Type *imag)
{
    DBCF_NAME(dbcF_cexpm1)(log2n,real,imag);
    *real=DBCF_ONE+*real;
}

/*
    Compute exp(2*pi*i*k/n), for 0<=k<b.
    At most O(log(n)) arithmetic operations (and, therefore, roundoff errors)
    are involved in calculating any given twiddle.
*/
static void DBCF_NAME(dbcF_compute_twiddles)(dbcf_index log2n,dbcf_index log2b,DBCF_Type *real,DBCF_Type *imag,int inverse)
{
    dbcf_index i;
    real[0]=DBCF_ZERO;
    imag[0]=DBCF_ZERO;
    for(i=0;i<log2b;++i)
    {
        dbcf_index j,k=DBCF_POW2(i);
        /*
            Note: accuracy seems to be slightly better when working
            with (cos-1,sin), rather than (cos,sin).
        */
        DBCF_Type x,y;
        DBCF_NAME(dbcF_cexpm1)(log2n-i,&x,&y);
        if(!inverse) y=-y;
        for(j=0;j<k;++j)
        {
            real[k+j]=(x*real[j]-y*imag[j])+(x+real[j]);
            imag[k+j]=(y*real[j]+x*imag[j])+(y+imag[j]);
        }
    }
    for(i=0;i<DBCF_POW2(log2b);++i)
        real[i]=DBCF_ONE+real[i];
}

#ifndef DBC_FFT_NO_NPOT
/*
    Compute exp(2*pi*i*(p/q))-1.
    We don't actually need to depend on <math.h>.
*/
static void DBCF_NAME(dbcF_cexpm1_npot)(dbcf_index p,dbcf_index q,DBCF_Type *real,DBCF_Type *imag)
{
#ifdef DBCF_cexpm1_npot
    DBCF_cexpm1_npot(p,q,real,imag);
#else
    dbcf_index i;
    DBCF_Type C=DBCF_ONE,S=DBCF_ONE;
    DBCF_Type x=DBCF_LITERAL(6.283185307179586476925286766559005768)*(DBCF_Type)p/(DBCF_Type)q,x2=x*x;
    DBCF_Type I=DBCF_LITERAL(32.0);
    for(i=32;i>=0;--i)
    {
        DBCF_Type J=(DBCF_LITERAL(2.0)*I+DBCF_LITERAL(3.0)),K=I+I+DBCF_LITERAL(3.0);
        J=J*J;
        C=DBCF_LITERAL(1.0)-x2*C/(J+K);
        S=DBCF_LITERAL(1.0)-x2*S/(J-K);
        I=I-DBCF_LITERAL(1.0);
    }
    C=-C*DBCF_LITERAL(0.5)*x2;
    S=S*x;
    *real=C;
    *imag=S;
#endif
}

static void DBCF_NAME(dbcF_compute_twiddles_npot)(dbcf_index n,DBCF_Type *real,DBCF_Type *imag,int inverse)
{
    /* Note: always gets called with even n. */
    dbcf_index i,j,k,m=n>>1,h=(m+2)>>1;
    if(n<1) return;
    real[0]=DBCF_ZERO;
    imag[0]=DBCF_ZERO;
    /*
        Note: accuracy is somewhat better when working with (cos-1,sin),
        rather than (cos,sin).
    */
    for(i=1;i<h;i*=2)
    {
        DBCF_Type X,Y;
        DBCF_NAME(dbcF_cexpm1_npot)(i,n,&X,&Y);
        if(!inverse) Y=-Y;
        j=(h<i*2?h-i:i);
        for(k=0;k<j;++k)
        {
            real[i+k]=(X*real[k]-Y*imag[k])+(X+real[k]);
            imag[i+k]=(Y*real[k]+X*imag[k])+(Y+imag[k]);
        }
    }
    for(i=0;i<h;++i) real[i]=DBCF_ONE+real[i];
    for(i=h;i<m;++i)
    {
        real[i]=-real[m-i];
        imag[i]= imag[m-i];
    }
    for(i=0;i<m;++i)
    {
        real[m+i]=-real[i];
        imag[m+i]=-imag[i];
    }
}
#endif

/* Bit-reversal permutations. */
static void DBCF_NAME(dbcF_bitreversal_swap)(dbcf_index log2n,DBCF_Type *src,dbcf_index src_stride,DBCF_Type *dst,dbcf_index dst_stride)
{
    dbcf_index i,n=DBCF_POW2(log2n),h=n>>1;
    if(log2n<=8)
    {
#ifndef DBC_FFT_NO_BITREVERSE_TABLE
        const unsigned char *idx=dbcF_bitreverse_table+DBCF_POW2(log2n);
#endif
        for(i=0;i<n;++i)
        {
#ifndef DBC_FFT_NO_BITREVERSE_TABLE
            dbcf_index j=(dbcf_index)(idx[i]);
#else
            dbcf_index j=dbcF_bitreverse(i,log2n);
#endif
            DBCF_Type x=src[i*src_stride];
            DBCF_Type y=dst[j*dst_stride];
            src[i*src_stride]=y;
            dst[j*dst_stride]=x;
        }
    }
    else
    {
        DBCF_NAME(dbcF_bitreversal_swap)(log2n-1,src           ,2*src_stride,dst             ,dst_stride);
        DBCF_NAME(dbcF_bitreversal_swap)(log2n-1,src+src_stride,2*src_stride,dst+h*dst_stride,dst_stride);
    }
}

static void DBCF_NAME(dbcF_bitreversal_permutation)(dbcf_index log2n,const DBCF_Type *src,dbcf_index src_stride,DBCF_Type *dst,dbcf_index dst_stride,DBCF_Type *tmp)
{
    dbcf_index i,n=DBCF_POW2(log2n),h=n>>1;
    if(src_stride==0)
    {
        DBCF_Type x=src[0];
        for(i=0;i<n;++i)
            dst[i*dst_stride]=x;
    }
    else if(src==dst)
    {
        /* In-place case. */
        if(log2n<=8)
        {
#ifndef DBC_FFT_NO_BITREVERSE_TABLE
            const unsigned char *idx=dbcF_bitreverse_table+DBCF_POW2(log2n);
#endif
            for(i=0;i<n;++i)
            {
#ifndef DBC_FFT_NO_BITREVERSE_TABLE
                dbcf_index j=(dbcf_index)(idx[i]);
#else
                dbcf_index j=dbcF_bitreverse(i,log2n);
#endif
                if(i<j)
                {
                    DBCF_Type x=dst[i*dst_stride];
                    DBCF_Type y=dst[j*dst_stride];
                    dst[i*dst_stride]=y;
                    dst[j*dst_stride]=x;
                }
            }
        }
        else if(log2n<=2*(DBCF_Q)+2||log2n<=16)
        {
            /* Exchange 0X...X1's and 1X...X0's */
            DBCF_NAME(dbcF_bitreversal_swap)(log2n-2,dst+dst_stride,2*dst_stride,dst+h*dst_stride,2*dst_stride);
            /* Reverse 0X...X0's */
            DBCF_NAME(dbcF_bitreversal_permutation)(log2n-2,dst                 ,2*dst_stride,dst                 ,2*dst_stride,tmp);
            /* Reverse 1X...X1's */
            DBCF_NAME(dbcF_bitreversal_permutation)(log2n-2,dst+(h+1)*dst_stride,2*dst_stride,dst+(h+1)*dst_stride,2*dst_stride,tmp);
        }
        else
        {
            /*
                The algorithm is based on
                "Towards an Optimal Bit-Reversal Permutation Program"
                by Larry Carter and Kang Su Gatlin.
            */
            dbcf_index a,b,c,log2m=log2n-2*(DBCF_Q);
            dbcf_index m=DBCF_POW2(log2m);
            dbcf_index pow2q=DBCF_POW2(DBCF_Q);
            for(b=0;b<m;++b)
            {
                dbcf_index ib=dbcF_bitreverse(b,log2m);
                if(ib<b) continue;
                for(a=0;a<pow2q;++a)
                    for(c=0;c<pow2q;++c)
                        tmp[(a<<(DBCF_Q))^c]=dst[((a<<(log2n-(DBCF_Q)))^(b<<(DBCF_Q))^c)*dst_stride];
                for(c=0;c<pow2q;++c)
                {
                    dbcf_index ic=dbcF_bitreverse(c,(DBCF_Q));
                    for(a=0;a<pow2q;++a)
                    {

                        dbcf_index ia=dbcF_bitreverse(a,(DBCF_Q));
                        DBCF_Type t;
                        i=(ic<<(log2n-(DBCF_Q)))^(ib<<(DBCF_Q))^ia;
                        t=dst[i*dst_stride];
                        dst[i*dst_stride]=tmp[(a<<(DBCF_Q))^c];
                        tmp[(a<<(DBCF_Q))^c]=t;
                    }
                }
                if(b!=ib)
                    for(a=0;a<pow2q;++a)
                        for(c=0;c<pow2q;++c)
                            dst[((a<<(log2n-(DBCF_Q)))^(b<<(DBCF_Q))^c)*dst_stride]=tmp[(a<<(DBCF_Q))^c];
            }
        }
    }
    else
    {
        if(log2n<=8)
        {
#ifndef DBC_FFT_NO_BITREVERSE_TABLE
            const unsigned char *idx=dbcF_bitreverse_table+DBCF_POW2(log2n);
#endif
            for(i=0;i<n;++i)
            {
#ifndef DBC_FFT_NO_BITREVERSE_TABLE
                dbcf_index j=(dbcf_index)(idx[i]);
#else
                dbcf_index j=dbcF_bitreverse(i,log2n);
#endif
                dst[j*dst_stride]=src[i*src_stride];
            }
        }
        else if(log2n<=16)
        {
            DBCF_NAME(dbcF_bitreversal_permutation)(log2n-1,src           ,2*src_stride,dst             ,dst_stride,tmp);
            DBCF_NAME(dbcF_bitreversal_permutation)(log2n-1,src+src_stride,2*src_stride,dst+h*dst_stride,dst_stride,tmp);
        }
        else
        {
            /*
                Do one pass, and call the in-place case. This turns out to be faster,
                especially for smaller Q.
            */
            for(i=0;i<h;++i)
            {
                dst[ i   *dst_stride]=src[(2*i  )*src_stride];
                dst[(i+h)*dst_stride]=src[(2*i+1)*src_stride];
            }
            DBCF_NAME(dbcF_bitreversal_permutation)(log2n-1,dst             ,dst_stride,dst             ,dst_stride,tmp);
            DBCF_NAME(dbcF_bitreversal_permutation)(log2n-1,dst+h*dst_stride,dst_stride,dst+h*dst_stride,dst_stride,tmp);
        }
    }
}

/* Hand-coded (I)FFT for size 8. */
static void DBCF_NAME(dbcF_fft8)(
    DBCF_Type *real,DBCF_Type *imag,
    dbcf_index real_stride,dbcf_index imag_stride,
    int inverse,DBCF_Type c)
{
    DBCF_Type r0,r1,r2,r3,r4,r5,r6,r7;
    DBCF_Type i0,i1,i2,i3,i4,i5,i6,i7;
    DBCF_Type R0,R1,R2,R3,R4,R5,R6,R7;
    DBCF_Type I0,I1,I2,I3,I4,I5,I6,I7;
    DBCF_Type p5,m5,p7,m7;
    r0=real[0*real_stride];i0=imag[0*imag_stride];
    r1=real[1*real_stride];i1=imag[1*imag_stride];
    r2=real[2*real_stride];i2=imag[2*imag_stride];
    r3=real[3*real_stride];i3=imag[3*imag_stride];
    r4=real[4*real_stride];i4=imag[4*imag_stride];
    r5=real[5*real_stride];i5=imag[5*imag_stride];
    r6=real[6*real_stride];i6=imag[6*imag_stride];
    r7=real[7*real_stride];i7=imag[7*imag_stride];
    R0=r0+r1;R1=r0-r1;I0=i0+i1;I1=i0-i1;
    R2=r2+r3;R3=r2-r3;I2=i2+i3;I3=i2-i3;
    R4=r4+r5;R5=r4-r5;I4=i4+i5;I5=i4-i5;
    R6=r6+r7;R7=r6-r7;I6=i6+i7;I7=i6-i7;
    if(!inverse)
    {
        r0=R0+R2;i0=I0+I2;
        r1=R1+I3;i1=I1-R3;
        r2=R0-R2;i2=I0-I2;
        r3=R1-I3;i3=I1+R3;
        r4=R4+R6;i4=I4+I6;
        r5=R5+I7;i5=I5-R7;
        r6=R4-R6;i6=I4-I6;
        r7=R5-I7;i7=I5+R7;
        p5=c*(r5+i5);m5=c*(r5-i5);
        p7=c*(r7+i7);m7=c*(r7-i7);
        real[0*real_stride]=r0+r4;imag[0*imag_stride]=i0+i4;
        real[1*real_stride]=r1+p5;imag[1*imag_stride]=i1-m5;
        real[2*real_stride]=r2+i6;imag[2*imag_stride]=i2-r6;
        real[3*real_stride]=r3-m7;imag[3*imag_stride]=i3-p7;
        real[4*real_stride]=r0-r4;imag[4*imag_stride]=i0-i4;
        real[5*real_stride]=r1-p5;imag[5*imag_stride]=i1+m5;
        real[6*real_stride]=r2-i6;imag[6*imag_stride]=i2+r6;
        real[7*real_stride]=r3+m7;imag[7*imag_stride]=i3+p7;
    }
    else
    {
        r0=R0+R2;i0=I0+I2;
        r1=R1-I3;i1=I1+R3;
        r2=R0-R2;i2=I0-I2;
        r3=R1+I3;i3=I1-R3;
        r4=R4+R6;i4=I4+I6;
        r5=R5-I7;i5=I5+R7;
        r6=R4-R6;i6=I4-I6;
        r7=R5+I7;i7=I5-R7;
        p5=c*(r5+i5);m5=c*(r5-i5);
        p7=c*(r7+i7);m7=c*(r7-i7);
        real[0*real_stride]=r0+r4;imag[0*imag_stride]=i0+i4;
        real[1*real_stride]=r1+m5;imag[1*imag_stride]=i1+p5;
        real[2*real_stride]=r2-i6;imag[2*imag_stride]=i2+r6;
        real[3*real_stride]=r3-p7;imag[3*imag_stride]=i3+m7;
        real[4*real_stride]=r0-r4;imag[4*imag_stride]=i0-i4;
        real[5*real_stride]=r1-m5;imag[5*imag_stride]=i1-p5;
        real[6*real_stride]=r2+i6;imag[6*imag_stride]=i2-r6;
        real[7*real_stride]=r3+p7;imag[7*imag_stride]=i3-m7;
    }
}

/*
    Compute a part of a butterfly of size n, on a block of size b,
    ensuring that each computed twiddle is touched by
    at most O(log(n)) arithmetic operations, to improve accuracy.
    The individual computed twiddles are products of precomputed
    twiddles (tr, ti) with recursively passed multipliers (C, S).
*/
static void DBCF_NAME(dbcF_butterfly_block)(
    dbcf_index log2n,
    dbcf_index log2b,
    DBCF_Type *LR,DBCF_Type *LI,
    DBCF_Type *HR,DBCF_Type *HI,
    dbcf_index real_stride,dbcf_index imag_stride,
    DBCF_Type C,DBCF_Type S,
    int inverse,
    const DBCF_Type *tr,const DBCF_Type *ti)
{
    dbcf_index b=DBCF_POW2(log2b),h=b>>1;
    DBCF_Type X,Y;
    if(log2b<=DBCF_TWIDDLES_BUF_LOG2)
    {
        /* The block is small, we have enough precomputed twiddles. */
        dbcf_index i,j=0,k=0;
        for(i=0;i<b;++i)
        {
            DBCF_Type c=C*tr[i]-S*ti[i];
            DBCF_Type s=S*tr[i]+C*ti[i];
            DBCF_Type xl=LR[j],yl=LI[k];
            DBCF_Type xr=HR[j],yr=HI[k];
            DBCF_Type x=c*xr-s*yr;
            DBCF_Type y=s*xr+c*yr;
            LR[j]=xl+x;
            LI[k]=yl+y;
            HR[j]=xl-x;
            HI[k]=yl-y;
            j+=real_stride;
            k+=imag_stride;
        }
    }
    else
    {
        /* The block is large, we process it's halves recursively. */
        DBCF_NAME(dbcF_cexp)(log2n-log2b+1,&X,&Y);
        if(!inverse) Y=-Y;
        DBCF_NAME(dbcF_butterfly_block)(log2n,log2b-1,LR              ,LI              ,HR              ,HI              ,real_stride,imag_stride,C      ,S      ,inverse,tr,ti);
        DBCF_NAME(dbcF_butterfly_block)(log2n,log2b-1,LR+h*real_stride,LI+h*imag_stride,HR+h*real_stride,HI+h*imag_stride,real_stride,imag_stride,C*X-S*Y,S*X+C*Y,inverse,tr,ti);
    }
}

/* Compute a butterfly pass over c blocks of size n. */
static void DBCF_NAME(dbcF_butterfly_pass)(
    dbcf_index log2n,
    dbcf_index log2c,
    DBCF_Type *real,DBCF_Type *imag,
    dbcf_index real_stride,dbcf_index imag_stride,
    int inverse,
	dbcf_index log2t,
    const DBCF_Type *tr,const DBCF_Type *ti)
{
    dbcf_index n=DBCF_POW2(log2n),h=n>>1;
    dbcf_index c=DBCF_POW2(log2c);
    dbcf_index i;
    DBCF_Type *LR=real,*HR=real+h*real_stride;
    DBCF_Type *LI=imag,*HI=imag+h*imag_stride;
    if(log2n==0) return;
    if(log2n-1<=log2t)
    {
        /*
            We have as much precomputed twiddles
            as the block needs, so we supply them directly.
        */
        if(h>1)
        {
            for(i=0;i<c;++i)
            {
                dbcf_index d,j=0,k=0;
                for(d=0;d<h;d+=2)
                {
                    /* Unroll x2. Slightly faster. */
                    DBCF_Type C,S,xl,yl,xr,yr,x,y;
                    C=tr[d];S=ti[d];
                    xl=LR[j];yl=LI[k];
                    xr=HR[j];yr=HI[k];
                    x=C*xr-S*yr;y=S*xr+C*yr;
                    LR[j]=xl+x;LI[k]=yl+y;
                    HR[j]=xl-x;HI[k]=yl-y;
                    j+=real_stride;k+=imag_stride;
                    C=tr[d+1];S=ti[d+1];
                    xl=LR[j];yl=LI[k];
                    xr=HR[j];yr=HI[k];
                    x=C*xr-S*yr;y=S*xr+C*yr;
                    LR[j]=xl+x;LI[k]=yl+y;
                    HR[j]=xl-x;HI[k]=yl-y;
                    j+=real_stride;k+=imag_stride;
                }
                LR+=n*real_stride;LI+=n*imag_stride;
                HR+=n*real_stride;HI+=n*imag_stride;
            }
        }
        else
        {
            for(i=0;i<c;++i)
            {
                DBCF_Type xl,yl,xr,yr;
                xl=LR[0];yl=LI[0];
                xr=HR[0];yr=HI[0];
                LR[0]=xl+xr;LI[0]=yl+yr;
                HR[0]=xl-xr;HI[0]=yl-yr;
                LR+=n*real_stride;LI+=n*imag_stride;
                HR+=n*real_stride;HI+=n*imag_stride;
            }
        }
    }
    else
    {
        /* We compute butterfly recursively, supplying only t precomputed twiddles. */
        for(i=0;i<c;++i)
        {
            DBCF_NAME(dbcF_butterfly_block)(log2n,log2n-1,LR,LI,HR,HI,real_stride,imag_stride,DBCF_ONE,DBCF_ZERO,inverse,tr,ti);
            LR+=n*real_stride;LI+=n*imag_stride;
            HR+=n*real_stride;HI+=n*imag_stride;
        }
    }
}

/* Compute a series of butterfly passes from (log2n-depth+1,log2c+depth+1) to (log2n,log2c). */
static void DBCF_NAME(dbcF_butterfly_multipass)(
    dbcf_index log2n,
    dbcf_index log2c,
    dbcf_index depth,
    DBCF_Type *real,DBCF_Type *imag,
    dbcf_index real_stride,dbcf_index imag_stride,
    int inverse,
    DBCF_Type *tr,DBCF_Type *ti)
{
    while(depth>0)
    {
        dbcf_index log2d,log2t;
#if defined(DBCF_butterfly_multipass_optimized)
        dbcf_index d=DBCF_butterfly_multipass_optimized(
            log2n,log2c,depth,
            real,imag,
            real_stride,imag_stride,
            inverse,
            tr,ti);
        if(d>0) {depth-=d;continue;}
#endif
        if(depth==log2n&&depth>=3)
        {
            dbcf_index j,m=DBCF_POW2(log2n+log2c-3);
            DBCF_NAME(dbcF_cexp)(3,tr,ti);
            for(j=0;j<m;++j)
                DBCF_NAME(dbcF_fft8)(real+8*real_stride*j,imag+8*imag_stride*j,real_stride,imag_stride,inverse,tr[0]);
            depth-=3;
            continue;
        }
        log2d=log2n-depth+1;
        log2t=(log2d-1<DBCF_TWIDDLES_BUF_LOG2?log2d-1:DBCF_TWIDDLES_BUF_LOG2);
        DBCF_NAME(dbcF_compute_twiddles)(log2d,log2t,tr,ti,inverse);
        DBCF_NAME(dbcF_butterfly_pass)(
            log2d,
            log2c+log2n-log2d,
            real,imag,
            real_stride,imag_stride,
            inverse,
            log2t,
            tr,ti);
        depth-=1;
    }
}

static void DBCF_NAME(dbcF_butterfly)(
    dbcf_index log2n,
    DBCF_Type *real,DBCF_Type *imag,
    dbcf_index real_stride,dbcf_index imag_stride,
    int inverse,
    DBCF_Type *tmp)
{
    DBCF_Type *tr=tmp;
    DBCF_Type *ti=tmp+DBCF_TWIDDLES_BUF_SIZE;
    if(log2n>12)
    {
        DBCF_NAME(dbcF_butterfly)(
            log2n-1,
            real,imag,
            real_stride,imag_stride,
            inverse,
            tmp);
        DBCF_NAME(dbcF_butterfly)(
            log2n-1,
            real+DBCF_POW2(log2n-1)*real_stride,imag+DBCF_POW2(log2n-1)*imag_stride,
            real_stride,imag_stride,
            inverse,
            tmp);
        DBCF_NAME(dbcF_butterfly_multipass)(
            log2n,0,1,
            real,imag,
            real_stride,imag_stride,
            inverse,
            tr,ti);
    }
    else
    {
        DBCF_NAME(dbcF_butterfly_multipass)(
            log2n,0,log2n,
            real,imag,
            real_stride,imag_stride,
            inverse,
            tr,ti);
    }
}

#ifdef DBCF_butterfly_multipass_optimized
/*
    Only provide (de)interleave if we have an optimized (SIMD)
    implementation, which actually cares about it.
*/
static void DBCF_NAME(dbcF_deinterleave)(DBCF_Type *dst,dbcf_index log2n,DBCF_Type *tmp)
{
    dbcf_index n=DBCF_POW2(log2n),h=n>>1;
    if(n<=2) return;
    if(n<=DBCF_TMP_BUF_SIZE)
    {
        dbcf_index i,h=n>>1;
        DBCF_Type *real=tmp,*imag=tmp+h;
        for(i=0;i<h;++i)
        {
            real[i]=dst[2*i+0];
            imag[i]=dst[2*i+1];
        }
        for(i=0;i<n;++i) dst[i]=tmp[i];
        return;
    }
    DBCF_NAME(dbcF_bitreversal_permutation)(log2n  ,dst  ,1,dst  ,1,tmp);
    DBCF_NAME(dbcF_bitreversal_permutation)(log2n-1,dst  ,1,dst  ,1,tmp);
    DBCF_NAME(dbcF_bitreversal_permutation)(log2n-1,dst+h,1,dst+h,1,tmp);
}

static void DBCF_NAME(dbcF_interleave)(DBCF_Type *dst,dbcf_index log2n,DBCF_Type *tmp)
{
    dbcf_index n=DBCF_POW2(log2n),h=n>>1;
    if(n<=2) return;
    if(n<=DBCF_TMP_BUF_SIZE)
    {
        dbcf_index i,h=n>>1;
        const DBCF_Type *real=dst,*imag=dst+h;
        for(i=0;i<h;++i)
        {
            tmp[2*i+0]=real[i];
            tmp[2*i+1]=imag[i];
        }
        for(i=0;i<n;++i) dst[i]=tmp[i];
        return;
    }
    DBCF_NAME(dbcF_bitreversal_permutation)(log2n-1,dst  ,1,dst  ,1,tmp);
    DBCF_NAME(dbcF_bitreversal_permutation)(log2n-1,dst+h,1,dst+h,1,tmp);
    DBCF_NAME(dbcF_bitreversal_permutation)(log2n  ,dst  ,1,dst  ,1,tmp);
}
#endif /* DBCF_butterfly_multipass_optimized */

/* Power-of-2 case. */
static int DBCF_NAME(dbcF_fft_pot)(
    dbcf_index num_elements,
    const DBCF_Type *src_real,const DBCF_Type *src_imag,
    dbcf_index src_real_stride,dbcf_index src_imag_stride,
          DBCF_Type *dst_real,      DBCF_Type *dst_imag,
    dbcf_index dst_real_stride,dbcf_index dst_imag_stride,
    int inverse,
    DBCF_Type scale)
{
#ifndef DBC_FFT_NO_SIMD
    DBCF_ALIGNED(32) DBCF_Type dummy[2]={DBCF_ZERO,DBCF_ZERO};
    DBCF_ALIGNED(32) DBCF_Type tmp[DBCF_TMP_BUF_SIZE];
#else
    DBCF_Type dummy[2]={DBCF_ZERO,DBCF_ZERO};
    DBCF_Type tmp[DBCF_TMP_BUF_SIZE];
#endif
    dbcf_index n=num_elements;
    dbcf_index log2n=(dbcf_index)-1;
    dbcf_index i;
#ifdef DBCF_butterfly_multipass_optimized
    /* Deinterleave dst to make better use of SIMD. */
    int needs_deinterleave=(dst_real_stride==2&&dst_imag_stride==2&&dst_imag==dst_real+1&&num_elements>16);
#endif
    while(n) {n>>=1;++log2n;}
    if(!src_real) {src_real=dummy  ;src_real_stride=0;}
    if(!src_imag) {src_imag=dummy+1;src_imag_stride=0;}
    DBCF_NAME(dbcF_bitreversal_permutation)(log2n,src_real,src_real_stride,dst_real,dst_real_stride,tmp);
    DBCF_NAME(dbcF_bitreversal_permutation)(log2n,src_imag,src_imag_stride,dst_imag,dst_imag_stride,tmp);
#ifdef DBCF_butterfly_multipass_optimized
    if(needs_deinterleave)
    {
        DBCF_NAME(dbcF_deinterleave)(dst_real,log2n+1,tmp);
        DBCF_NAME(dbcF_butterfly)(
            log2n,
            dst_real,dst_real+num_elements,
            1,1,
            inverse,tmp);
    }
    else
#endif
    DBCF_NAME(dbcF_butterfly)(
        log2n,
        dst_real,dst_imag,
        dst_real_stride,dst_imag_stride,
        inverse,tmp);
#ifdef DBCF_butterfly_multipass_optimized
    if(needs_deinterleave) DBCF_NAME(dbcF_interleave)(dst_real,log2n+1,tmp);
#endif
    if(scale!=DBCF_ONE) for(i=0;i<num_elements;++i)
    {
        dst_real[i]=dst_real[i]*scale;
        dst_imag[i]=dst_imag[i]*scale;
    }
    return 0;
}

#ifndef DBC_FFT_NO_NPOT
/* Non-power-of-2 case. */
static int DBCF_NAME(dbcF_fft_npot)(
    dbcf_index n,
    const DBCF_Type *src_real,const DBCF_Type *src_imag,
    dbcf_index src_real_stride,dbcf_index src_imag_stride,
          DBCF_Type *dst_real,      DBCF_Type *dst_imag,
    dbcf_index dst_real_stride,dbcf_index dst_imag_stride,
    int inverse,
    DBCF_Type scale)
{
    unsigned char *buf,*mem=0;
    DBCF_Type *ar,*ai,*br,*bi,*tr,*ti;
    /*
        M has the same value as m, but different type. This avoids
        dbcf_index->DBCF_Type cast, in case the custom type does
        not provide it.
    */
    DBCF_Type M=DBCF_ONE;
    dbcf_index i,j,log2m=0,m;
#ifdef DBCF_butterfly_multipass_optimized
    dbcf_index alignment=64;
#else
    dbcf_index alignment=0;
#endif
    while(DBCF_POW2(log2m)<2*n-1) {++log2m;M=M+M;}
    m=DBCF_POW2(log2m);
    if(!(mem=(unsigned char*)dbcf_malloc((4*m+4*n)*(dbcf_index)sizeof(DBCF_Type)+alignment))) return DBCF_ERROR_OUT_OF_MEMORY;
    buf=mem;
    if(alignment)
    {
        dbcf_index offset=((dbcf_index)buf)&(alignment-1);
        if(offset) buf+=alignment-offset;
    }
    ar=(DBCF_Type*)buf+0*m;
    ai=(DBCF_Type*)buf+1*m;
    br=(DBCF_Type*)buf+2*m;
    bi=(DBCF_Type*)buf+3*m;
    tr=(DBCF_Type*)buf+4*m+0*n;
    ti=(DBCF_Type*)buf+4*m+2*n;
    DBCF_NAME(dbcF_compute_twiddles_npot)(2*n,tr,ti,inverse);
    for(i=0,j=0;i<n;++i)
    {
        DBCF_Type c=tr[j],s=ti[j];
        DBCF_Type x=src_real[i*src_real_stride],y=src_imag[i*src_imag_stride];
        ar[i]=x*c-y*s;
        ai[i]=x*s+y*c;
        br[i]= c;
        bi[i]=-s;
        if(i>0)
        {
            br[m-i]= c;
            bi[m-i]=-s;
        }
        j+=(2*i+1);
        if(j>=2*n) j-=2*n;
    }
    for(i=n;i<m;++i)
    {
        ar[i]=DBCF_ZERO;
        ai[i]=DBCF_ZERO;
    }
    for(i=n;i<=m-n;++i)
    {
        br[i]=DBCF_ZERO;
        bi[i]=DBCF_ZERO;
    }
    /*
        Note: the scale factors for FFTs are (1/M,1,scale), rather than, say,
        (1,1,scale/M). This helps to keep intermediate results from
        overflowing/underflowing, when the range is limited (fixed-point,
        maybe half-floats).
    */
    DBCF_NAME(dbcF_fft_pot)(m,ar,ai,1,1,ar,ai,1,1,0,DBCF_ONE/M);
    DBCF_NAME(dbcF_fft_pot)(m,br,bi,1,1,br,bi,1,1,0,DBCF_ONE);
    for(i=0;i<m;++i)
    {
        DBCF_Type c=br[i],s=bi[i],x=ar[i],y=ai[i];
        ar[i]=c*x-s*y;
        ai[i]=c*y+s*x;
    }
    DBCF_NAME(dbcF_fft_pot)(m,ar,ai,1,1,ar,ai,1,1,1,scale);
    for(i=0,j=0;i<n;++i)
    {
        DBCF_Type c=tr[j],s=ti[j],x=ar[i],y=ai[i];
        dst_real[i*dst_real_stride]=c*x-s*y;
        dst_imag[i*dst_imag_stride]=c*y+s*x;
        j+=(2*i+1);
        if(j>=2*n) j-=2*n;
    }
    dbcf_free(mem);
    return 0;
}
#endif /* DBC_FFT_NO_NPOT */

static int DBCF_NAME(dbcF_fft)(
    dbcf_index num_elements,
    const DBCF_Type *src_real,const DBCF_Type *src_imag,
    dbcf_index src_real_stride,dbcf_index src_imag_stride,
          DBCF_Type *dst_real,      DBCF_Type *dst_imag,
    dbcf_index dst_real_stride,dbcf_index dst_imag_stride,
    int inverse,
    DBCF_Type scale)
{
    dbcF_init();
    if(num_elements<1) return 0;
    if(src_real==dst_real&&src_real_stride!=dst_real_stride) return DBCF_ERROR_INVALID_ARGUMENT;
    if(src_imag==dst_imag&&src_imag_stride!=dst_imag_stride) return DBCF_ERROR_INVALID_ARGUMENT;
    if(src_imag==dst_real) return DBCF_ERROR_INVALID_ARGUMENT;
    if(src_real==dst_imag) return DBCF_ERROR_INVALID_ARGUMENT;
    if(num_elements&(num_elements-1))
#ifndef DBC_FFT_NO_NPOT
        return DBCF_NAME(dbcF_fft_npot)(
            num_elements,
            src_real,src_imag,
            src_real_stride,src_imag_stride,
            dst_real,dst_imag,
            dst_real_stride,dst_imag_stride,
            inverse,
            scale);
#else
        return DBCF_ERROR_INVALID_ARGUMENT;
#endif /* DBC_FFT_NO_NPOT*/
    return DBCF_NAME(dbcF_fft_pot)(
        num_elements,
        src_real,src_imag,
        src_real_stride,src_imag_stride,
        dst_real,dst_imag,
        dst_real_stride,dst_imag_stride,
        inverse,
        scale);
}

DBCF_DEF int DBCF_NAME2(dbc_fft,c)(
    dbcf_index num_elements,
    const DBCF_Type *src_real,const DBCF_Type *src_imag,
          DBCF_Type *dst_real,      DBCF_Type *dst_imag,
    DBCF_Type scale)
{
    return DBCF_NAME(dbcF_fft)(num_elements,
        src_real,src_imag,
        1,1,
        dst_real,dst_imag,
        1,1,
        0,
        scale);
}

DBCF_DEF int DBCF_NAME2(dbc_ifft,c)(
    dbcf_index num_elements,
    const DBCF_Type *src_real,const DBCF_Type *src_imag,
          DBCF_Type *dst_real,      DBCF_Type *dst_imag,
    DBCF_Type scale)
{
    return DBCF_NAME(dbcF_fft)(num_elements,
        src_real,src_imag,
        1,1,
        dst_real,dst_imag,
        1,1,
        1,
        scale);
}

DBCF_DEF int DBCF_NAME2(dbc_fft,i)(
    dbcf_index num_elements,
    const DBCF_Type *src,
          DBCF_Type *dst,
    DBCF_Type scale)
{
    return DBCF_NAME(dbcF_fft)(num_elements,
        src,(src?src+1:src),
        (src?2:0),(src?2:0),
        dst,dst+1,
        (src?2:0),(src?2:0),
        0,
        scale);
}

DBCF_DEF int DBCF_NAME2(dbc_ifft,i)(
    dbcf_index num_elements,
    const DBCF_Type *src,
          DBCF_Type *dst,
    DBCF_Type scale)
{
    return DBCF_NAME(dbcF_fft)(num_elements,
        src,(src?src+1:src),
        (src?2:0),(src?2:0),
        dst,dst+1,
        (src?2:0),(src?2:0),
        1,
        scale);
}

DBCF_DEF int DBCF_NAME2(dbc_fft,s)(
    dbcf_index num_elements,
    const DBCF_Type *src_real,const DBCF_Type *src_imag,
    dbcf_index src_real_stride,dbcf_index src_imag_stride,
          DBCF_Type *dst_real,      DBCF_Type *dst_imag,
    dbcf_index dst_real_stride,dbcf_index dst_imag_stride,
    DBCF_Type scale)
{
    return DBCF_NAME(dbcF_fft)(num_elements,
        src_real,src_imag,
        src_real_stride,src_imag_stride,
        dst_real,dst_imag,
        dst_real_stride,dst_imag_stride,
        0,
        scale);
}

DBCF_DEF int DBCF_NAME2(dbc_ifft,s)(
    dbcf_index num_elements,
    const DBCF_Type *src_real,const DBCF_Type *src_imag,
    dbcf_index src_real_stride,dbcf_index src_imag_stride,
          DBCF_Type *dst_real,      DBCF_Type *dst_imag,
    dbcf_index dst_real_stride,dbcf_index dst_imag_stride,
    DBCF_Type scale)
{
    return DBCF_NAME(dbcF_fft)(num_elements,
        src_real,src_imag,
        src_real_stride,src_imag_stride,
        dst_real,dst_imag,
        dst_real_stride,dst_imag_stride,
        1,
        scale);
}

#if defined(__cplusplus) && !defined(DBC_FFT_NO_CPP_OVERLOADS)
DBCF_DEF int dbc_fft(
    dbcf_index num_elements,
    const DBCF_Type *src_real,const DBCF_Type *src_imag,
          DBCF_Type *dst_real,      DBCF_Type *dst_imag,
    DBCF_Type scale)
{
    return DBCF_NAME(dbcF_fft)(num_elements,
        src_real,src_imag,
        1,1,
        dst_real,dst_imag,
        1,1,
        0,
        scale);
}

DBCF_DEF int dbc_ifft(
    dbcf_index num_elements,
    const DBCF_Type *src_real,const DBCF_Type *src_imag,
          DBCF_Type *dst_real,      DBCF_Type *dst_imag,
    DBCF_Type scale)
{
    return DBCF_NAME(dbcF_fft)(num_elements,
        src_real,src_imag,
        1,1,
        dst_real,dst_imag,
        1,1,
        1,
        scale);
}

DBCF_DEF int dbc_fft(
    dbcf_index num_elements,
    const DBCF_Type *src,
          DBCF_Type *dst,
    DBCF_Type scale)
{
    return DBCF_NAME(dbcF_fft)(num_elements,
        src,(src?src+1:src),
        (src?2:0),(src?2:0),
        dst,dst+1,
        (src?2:0),(src?2:0),
        0,
        scale);
}

DBCF_DEF int dbc_ifft(
    dbcf_index num_elements,
    const DBCF_Type *src,
          DBCF_Type *dst,
    DBCF_Type scale)
{
    return DBCF_NAME(dbcF_fft)(num_elements,
        src,(src?src+1:src),
        (src?2:0),(src?2:0),
        dst,dst+1,
        (src?2:0),(src?2:0),
        1,
        scale);
}

DBCF_DEF int dbc_fft(
    dbcf_index num_elements,
    const DBCF_Type *src_real,const DBCF_Type *src_imag,
    dbcf_index src_real_stride,dbcf_index src_imag_stride,
          DBCF_Type *dst_real,      DBCF_Type *dst_imag,
    dbcf_index dst_real_stride,dbcf_index dst_imag_stride,
    DBCF_Type scale)
{
    return DBCF_NAME(dbcF_fft)(num_elements,
        src_real,src_imag,
        src_real_stride,src_imag_stride,
        dst_real,dst_imag,
        dst_real_stride,dst_imag_stride,
        0,
        scale);
}

DBCF_DEF int dbc_ifft(
    dbcf_index num_elements,
    const DBCF_Type *src_real,const DBCF_Type *src_imag,
    dbcf_index src_real_stride,dbcf_index src_imag_stride,
          DBCF_Type *dst_real,      DBCF_Type *dst_imag,
    dbcf_index dst_real_stride,dbcf_index dst_imag_stride,
    DBCF_Type scale)
{
    return DBCF_NAME(dbcF_fft)(num_elements,
        src_real,src_imag,
        src_real_stride,src_imag_stride,
        dst_real,dst_imag,
        dst_real_stride,dst_imag_stride,
        1,
        scale);
}
#endif /* defined(__cplusplus) && !defined(DBC_FFT_NO_CPP_OVERLOADS) */

#endif /* DBC_FFT_INSTANTIATION */
