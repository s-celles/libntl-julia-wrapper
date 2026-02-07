/**
 * LibNTL Julia Wrapper
 *
 * CxxWrap-based wrapper for NTL (Number Theory Library) providing
 * Julia bindings for ZZ, ZZ_p, ZZX, ZZ_pX, Vec<ZZ>, Vec<ZZ_p>, Mat<ZZ>,
 * GF2, GF2X, Vec<GF2>, Mat<GF2>, zz_p, zz_pX, ZZ_pE, ZZ_pEX, RR,
 * PrimeSeq, and number theory functions (PowerMod, ProbPrime, etc.).
 */

#include <jlcxx/jlcxx.hpp>
#include <jlcxx/stl.hpp>

// Core NTL types
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/mat_ZZ.h>

// Binary field types (GF2)
#include <NTL/GF2.h>
#include <NTL/GF2X.h>
#include <NTL/vec_GF2.h>
#include <NTL/mat_GF2.h>

// Small prime types (single-precision)
#include <NTL/lzz_p.h>
#include <NTL/lzz_pX.h>

// Extension field types
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>

// Arbitrary precision floating point
#include <NTL/RR.h>

// Factoring (for DetIrredTest, etc.)
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/lzz_pXFactoring.h>

#include <sstream>
#include <string>
#include <stdexcept>

using namespace NTL;

/**
 * Module entry point for CxxWrap
 *
 * This function defines all Julia bindings for NTL types.
 * Types and functions are registered here and become available
 * when the Julia module loads the shared library.
 */
JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    // =========================================================================
    // ZZ Type (Arbitrary-Precision Integers)
    // =========================================================================

    mod.add_type<ZZ>("ZZ")
        .constructor<>()
        .constructor<long>()
        .method("__copy__", [](const ZZ& a) { return ZZ(a); });

    // ZZ String Conversion
    mod.method("ZZ_from_string", [](const std::string& s) {
        ZZ z;
        std::istringstream iss(s);
        iss >> z;
        if (iss.fail()) {
            throw std::invalid_argument("Invalid integer string: " + s);
        }
        return z;
    });

    mod.method("ZZ_to_string", [](const ZZ& z) {
        std::ostringstream oss;
        oss << z;
        return oss.str();
    });

    // ZZ Arithmetic Operations
    mod.method("ZZ_add", [](const ZZ& a, const ZZ& b) { return a + b; });
    mod.method("ZZ_sub", [](const ZZ& a, const ZZ& b) { return a - b; });
    mod.method("ZZ_mul", [](const ZZ& a, const ZZ& b) { return a * b; });

    mod.method("ZZ_div", [](const ZZ& a, const ZZ& b) {
        if (IsZero(b)) throw std::domain_error("Division by zero");
        return a / b;
    });

    mod.method("ZZ_rem", [](const ZZ& a, const ZZ& b) {
        if (IsZero(b)) throw std::domain_error("Division by zero");
        return a % b;
    });

    mod.method("ZZ_divrem", [](const ZZ& a, const ZZ& b) {
        if (IsZero(b)) throw std::domain_error("Division by zero");
        ZZ q, r;
        DivRem(q, r, a, b);
        return std::make_tuple(q, r);
    });

    mod.method("ZZ_power", [](const ZZ& a, long e) {
        if (e < 0) throw std::domain_error("Negative exponent");
        return power(a, e);
    });

    mod.method("ZZ_negate", [](const ZZ& a) { return -a; });
    mod.method("ZZ_abs", [](const ZZ& a) { return abs(a); });

    // ZZ GCD Operations
    mod.method("ZZ_gcd", [](const ZZ& a, const ZZ& b) {
        return GCD(a, b);
    });

    mod.method("ZZ_gcdx", [](const ZZ& a, const ZZ& b) {
        ZZ d, s, t;
        XGCD(d, s, t, a, b);
        return std::make_tuple(d, s, t);
    });

    // ZZ Comparison Operations
    mod.method("ZZ_equal", [](const ZZ& a, const ZZ& b) { return a == b; });
    mod.method("ZZ_less", [](const ZZ& a, const ZZ& b) { return a < b; });
    mod.method("ZZ_lesseq", [](const ZZ& a, const ZZ& b) { return a <= b; });

    // ZZ Predicates
    mod.method("ZZ_iszero", [](const ZZ& a) { return IsZero(a); });
    mod.method("ZZ_isone", [](const ZZ& a) { return IsOne(a); });
    mod.method("ZZ_sign", [](const ZZ& a) { return sign(a); });
    mod.method("ZZ_isodd", [](const ZZ& a) { return IsOdd(a); });

    // ZZ Size Queries
    mod.method("ZZ_numbits", [](const ZZ& a) { return NumBits(a); });
    mod.method("ZZ_numbytes", [](const ZZ& a) { return NumBytes(a); });

    // =========================================================================
    // ZZ_p Type (Integers Modulo p)
    // =========================================================================

    mod.add_type<ZZ_p>("ZZ_p")
        .constructor<>()
        .constructor<long>()
        .method("__copy__", [](const ZZ_p& a) { return ZZ_p(a); });

    // ZZ_p Modulus Management
    mod.method("ZZ_p_init", [](const ZZ& p) {
        if (p <= 1) throw std::domain_error("Modulus must be > 1");
        ZZ_p::init(p);
    });

    mod.method("ZZ_p_modulus", []() {
        return ZZ_p::modulus();
    });

    // ZZ_pContext Type
    mod.add_type<ZZ_pContext>("ZZ_pContext")
        .constructor<>()
        .constructor<const ZZ&>();

    mod.method("ZZ_pContext_save", [](ZZ_pContext& ctx) {
        ctx.save();
    });

    mod.method("ZZ_pContext_restore", [](const ZZ_pContext& ctx) {
        ctx.restore();
    });

    // ZZ_p Representation
    mod.method("ZZ_p_rep", [](const ZZ_p& a) {
        return rep(a);
    });

    // ZZ_p Arithmetic Operations
    mod.method("ZZ_p_add", [](const ZZ_p& a, const ZZ_p& b) { return a + b; });
    mod.method("ZZ_p_sub", [](const ZZ_p& a, const ZZ_p& b) { return a - b; });
    mod.method("ZZ_p_mul", [](const ZZ_p& a, const ZZ_p& b) { return a * b; });
    mod.method("ZZ_p_negate", [](const ZZ_p& a) { return -a; });

    mod.method("ZZ_p_inv", [](const ZZ_p& a) {
        if (IsZero(a)) {
            throw std::domain_error("Inverse of zero");
        }
        return inv(a);
    });

    mod.method("ZZ_p_div", [](const ZZ_p& a, const ZZ_p& b) {
        if (IsZero(b)) {
            throw std::domain_error("Division by zero");
        }
        return a / b;
    });

    mod.method("ZZ_p_power", [](const ZZ_p& a, long e) {
        return power(a, e);
    });

    mod.method("ZZ_p_power_ZZ", [](const ZZ_p& a, const ZZ& e) {
        return power(a, e);
    });

    // ZZ_p Predicates
    mod.method("ZZ_p_iszero", [](const ZZ_p& a) { return IsZero(a); });
    mod.method("ZZ_p_isone", [](const ZZ_p& a) { return IsOne(a); });

    // =========================================================================
    // ZZX Type (Polynomials over Z)
    // =========================================================================

    mod.add_type<ZZX>("ZZX")
        .constructor<>()
        .constructor<long>()
        .constructor<const ZZ&>()
        .method("__copy__", [](const ZZX& f) { return ZZX(f); });

    // ZZX Coefficient Access
    mod.method("ZZX_deg", [](const ZZX& f) { return deg(f); });

    mod.method("ZZX_coeff", [](const ZZX& f, long i) {
        return coeff(f, i);
    });

    mod.method("ZZX_setcoeff", [](ZZX& f, long i, const ZZ& c) {
        SetCoeff(f, i, c);
    });

    mod.method("ZZX_leadcoeff", [](const ZZX& f) { return LeadCoeff(f); });
    mod.method("ZZX_constterm", [](const ZZX& f) { return ConstTerm(f); });

    // ZZX Arithmetic Operations
    mod.method("ZZX_add", [](const ZZX& f, const ZZX& g) { return f + g; });
    mod.method("ZZX_sub", [](const ZZX& f, const ZZX& g) { return f - g; });
    mod.method("ZZX_mul", [](const ZZX& f, const ZZX& g) { return f * g; });
    mod.method("ZZX_mul_scalar", [](const ZZ& c, const ZZX& f) { return c * f; });
    mod.method("ZZX_negate", [](const ZZX& f) { return -f; });

    mod.method("ZZX_div", [](const ZZX& f, const ZZX& g) {
        if (IsZero(g)) throw std::domain_error("Division by zero polynomial");
        return f / g;
    });

    mod.method("ZZX_rem", [](const ZZX& f, const ZZX& g) {
        if (IsZero(g)) throw std::domain_error("Division by zero polynomial");
        return f % g;
    });

    mod.method("ZZX_divrem", [](const ZZX& f, const ZZX& g) {
        if (IsZero(g)) throw std::domain_error("Division by zero polynomial");
        ZZX q, r;
        DivRem(q, r, f, g);
        return std::make_tuple(q, r);
    });

    mod.method("ZZX_gcd", [](const ZZX& f, const ZZX& g) {
        return GCD(f, g);
    });

    // ZZX Polynomial Operations
    mod.method("ZZX_diff", [](const ZZX& f) {
        ZZX df;
        diff(df, f);
        return df;
    });

    mod.method("ZZX_content", [](const ZZX& f) {
        ZZ c;
        content(c, f);
        return c;
    });

    mod.method("ZZX_primpart", [](const ZZX& f) {
        ZZX pp;
        PrimitivePart(pp, f);
        return pp;
    });

    mod.method("ZZX_eval", [](const ZZX& f, const ZZ& x) {
        // Horner's method for polynomial evaluation
        if (IsZero(f)) return ZZ(0);
        long d = deg(f);
        ZZ result = coeff(f, d);
        for (long i = d - 1; i >= 0; i--) {
            result = result * x + coeff(f, i);
        }
        return result;
    });

    // ZZX Predicates and Conversion
    mod.method("ZZX_iszero", [](const ZZX& f) { return IsZero(f); });

    mod.method("ZZX_to_string", [](const ZZX& f) {
        std::ostringstream oss;
        oss << f;
        return oss.str();
    });

    // =========================================================================
    // Vec<ZZ> Type (Vector of Integers)
    // =========================================================================

    mod.add_type<Vec<ZZ>>("VecZZ")
        .constructor<>()
        .method("__copy__", [](const Vec<ZZ>& v) { return Vec<ZZ>(v); });

    mod.method("VecZZ_length", [](const Vec<ZZ>& v) { return v.length(); });

    mod.method("VecZZ_getindex", [](const Vec<ZZ>& v, long i) {
        // Julia is 1-indexed, NTL operator() is 1-indexed
        if (i < 1 || i > v.length()) {
            throw std::out_of_range("Index out of bounds");
        }
        return v(i);
    });

    mod.method("VecZZ_setindex!", [](Vec<ZZ>& v, long i, const ZZ& x) {
        if (i < 1 || i > v.length()) {
            throw std::out_of_range("Index out of bounds");
        }
        v(i) = x;
    });

    mod.method("VecZZ_setlength!", [](Vec<ZZ>& v, long n) {
        v.SetLength(n);
    });

    mod.method("VecZZ_append!", [](Vec<ZZ>& v, const ZZ& x) {
        append(v, x);
    });

    mod.method("VecZZ_to_string", [](const Vec<ZZ>& v) {
        std::ostringstream oss;
        oss << v;
        return oss.str();
    });

    // =========================================================================
    // Mat<ZZ> Type (Matrix of Integers)
    // =========================================================================

    mod.add_type<Mat<ZZ>>("MatZZ")
        .constructor<>()
        .method("__copy__", [](const Mat<ZZ>& m) { return Mat<ZZ>(m); });

    mod.method("MatZZ_nrows", [](const Mat<ZZ>& m) { return m.NumRows(); });
    mod.method("MatZZ_ncols", [](const Mat<ZZ>& m) { return m.NumCols(); });

    mod.method("MatZZ_getindex", [](const Mat<ZZ>& m, long i, long j) {
        // Julia is 1-indexed, NTL operator() is 1-indexed
        if (i < 1 || i > m.NumRows() || j < 1 || j > m.NumCols()) {
            throw std::out_of_range("Index out of bounds");
        }
        return m(i, j);
    });

    mod.method("MatZZ_setindex!", [](Mat<ZZ>& m, long i, long j, const ZZ& x) {
        if (i < 1 || i > m.NumRows() || j < 1 || j > m.NumCols()) {
            throw std::out_of_range("Index out of bounds");
        }
        m(i, j) = x;
    });

    mod.method("MatZZ_setdims!", [](Mat<ZZ>& m, long n, long k) {
        m.SetDims(n, k);
    });

    mod.method("MatZZ_mul", [](const Mat<ZZ>& a, const Mat<ZZ>& b) {
        Mat<ZZ> c;
        mul(c, a, b);
        return c;
    });

    mod.method("MatZZ_add", [](const Mat<ZZ>& a, const Mat<ZZ>& b) {
        Mat<ZZ> c;
        add(c, a, b);
        return c;
    });

    mod.method("MatZZ_sub", [](const Mat<ZZ>& a, const Mat<ZZ>& b) {
        Mat<ZZ> c;
        sub(c, a, b);
        return c;
    });

    mod.method("MatZZ_to_string", [](const Mat<ZZ>& m) {
        std::ostringstream oss;
        oss << m;
        return oss.str();
    });

    // =========================================================================
    // PrimeSeq (Prime Number Iterator)
    // =========================================================================

    mod.add_type<PrimeSeq>("PrimeSeq")
        .constructor<>();

    mod.method("PrimeSeq_next", [](PrimeSeq& ps) {
        return ps.next();
    });

    mod.method("PrimeSeq_reset", [](PrimeSeq& ps, long start) {
        ps.reset(start);
    });

    // =========================================================================
    // Number Theory Functions
    // =========================================================================

    // PowerMod: a^e mod n
    mod.method("ZZ_PowerMod", [](const ZZ& a, const ZZ& e, const ZZ& n) {
        if (n <= 1) throw std::domain_error("Modulus must be > 1");
        return PowerMod(a, e, n);
    });

    // PowerMod with long exponent
    mod.method("ZZ_PowerMod_long", [](const ZZ& a, long e, const ZZ& n) {
        if (n <= 1) throw std::domain_error("Modulus must be > 1");
        return PowerMod(a, e, n);
    });

    // Bit operations
    mod.method("ZZ_bit", [](const ZZ& a, long i) {
        return bit(a, i);
    });

    // RandomBnd: random in [0, n-1]
    mod.method("ZZ_RandomBnd", [](const ZZ& n) {
        if (n <= 0) throw std::domain_error("Bound must be > 0");
        return RandomBnd(n);
    });

    // RandomBits: random n-bit number
    mod.method("ZZ_RandomBits", [](long n) {
        if (n < 0) throw std::domain_error("Number of bits must be >= 0");
        return RandomBits_ZZ(n);
    });

    // ProbPrime: Miller-Rabin primality test
    mod.method("ZZ_ProbPrime", [](const ZZ& n, long num_trials) {
        return ProbPrime(n, num_trials);
    });

    // =========================================================================
    // ZZ_pX Type (Polynomials over Z/pZ)
    // =========================================================================

    mod.add_type<ZZ_pX>("ZZ_pX")
        .constructor<>()
        .constructor<long>()
        .constructor<const ZZ_p&>()
        .method("__copy__", [](const ZZ_pX& f) { return ZZ_pX(f); });

    // ZZ_pX Coefficient Access
    mod.method("ZZ_pX_deg", [](const ZZ_pX& f) { return deg(f); });

    mod.method("ZZ_pX_coeff", [](const ZZ_pX& f, long i) {
        return coeff(f, i);
    });

    mod.method("ZZ_pX_setcoeff", [](ZZ_pX& f, long i, const ZZ_p& c) {
        SetCoeff(f, i, c);
    });

    mod.method("ZZ_pX_setcoeff_long", [](ZZ_pX& f, long i, long c) {
        SetCoeff(f, i, c);
    });

    mod.method("ZZ_pX_leadcoeff", [](const ZZ_pX& f) { return LeadCoeff(f); });
    mod.method("ZZ_pX_constterm", [](const ZZ_pX& f) { return ConstTerm(f); });

    // ZZ_pX Arithmetic Operations
    mod.method("ZZ_pX_add", [](const ZZ_pX& f, const ZZ_pX& g) { return f + g; });
    mod.method("ZZ_pX_sub", [](const ZZ_pX& f, const ZZ_pX& g) { return f - g; });
    mod.method("ZZ_pX_mul", [](const ZZ_pX& f, const ZZ_pX& g) { return f * g; });
    mod.method("ZZ_pX_mul_scalar", [](const ZZ_p& c, const ZZ_pX& f) { return c * f; });
    mod.method("ZZ_pX_negate", [](const ZZ_pX& f) { return -f; });

    mod.method("ZZ_pX_div", [](const ZZ_pX& f, const ZZ_pX& g) {
        if (IsZero(g)) throw std::domain_error("Division by zero polynomial");
        return f / g;
    });

    mod.method("ZZ_pX_rem", [](const ZZ_pX& f, const ZZ_pX& g) {
        if (IsZero(g)) throw std::domain_error("Division by zero polynomial");
        return f % g;
    });

    mod.method("ZZ_pX_divrem", [](const ZZ_pX& f, const ZZ_pX& g) {
        if (IsZero(g)) throw std::domain_error("Division by zero polynomial");
        ZZ_pX q, r;
        DivRem(q, r, f, g);
        return std::make_tuple(q, r);
    });

    mod.method("ZZ_pX_gcd", [](const ZZ_pX& f, const ZZ_pX& g) {
        return GCD(f, g);
    });

    // ZZ_pX Polynomial Operations
    mod.method("ZZ_pX_diff", [](const ZZ_pX& f) {
        ZZ_pX df;
        diff(df, f);
        return df;
    });

    mod.method("ZZ_pX_eval", [](const ZZ_pX& f, const ZZ_p& x) {
        return eval(f, x);
    });

    // ZZ_pX Irreducibility Test
    mod.method("ZZ_pX_DetIrredTest", [](const ZZ_pX& f) {
        return DetIrredTest(f);
    });

    // ZZ_pX Predicates and Conversion
    mod.method("ZZ_pX_iszero", [](const ZZ_pX& f) { return IsZero(f); });
    mod.method("ZZ_pX_isone", [](const ZZ_pX& f) { return IsOne(f); });

    mod.method("ZZ_pX_to_string", [](const ZZ_pX& f) {
        std::ostringstream oss;
        oss << f;
        return oss.str();
    });

    // Build irreducible polynomial of given degree
    mod.method("ZZ_pX_BuildIrred", [](long n) {
        ZZ_pX f;
        BuildIrred(f, n);
        return f;
    });

    // =========================================================================
    // Vec<ZZ_p> Type (Vector over Z/pZ)
    // =========================================================================

    mod.add_type<Vec<ZZ_p>>("VecZZ_p")
        .constructor<>()
        .method("__copy__", [](const Vec<ZZ_p>& v) { return Vec<ZZ_p>(v); });

    mod.method("VecZZ_p_length", [](const Vec<ZZ_p>& v) { return v.length(); });

    mod.method("VecZZ_p_getindex", [](const Vec<ZZ_p>& v, long i) {
        if (i < 1 || i > v.length()) {
            throw std::out_of_range("Index out of bounds");
        }
        return v(i);
    });

    mod.method("VecZZ_p_setindex!", [](Vec<ZZ_p>& v, long i, const ZZ_p& x) {
        if (i < 1 || i > v.length()) {
            throw std::out_of_range("Index out of bounds");
        }
        v(i) = x;
    });

    mod.method("VecZZ_p_setlength!", [](Vec<ZZ_p>& v, long n) {
        v.SetLength(n);
    });

    mod.method("VecZZ_p_append!", [](Vec<ZZ_p>& v, const ZZ_p& x) {
        append(v, x);
    });

    mod.method("VecZZ_p_add", [](const Vec<ZZ_p>& a, const Vec<ZZ_p>& b) {
        Vec<ZZ_p> c;
        add(c, a, b);
        return c;
    });

    mod.method("VecZZ_p_sub", [](const Vec<ZZ_p>& a, const Vec<ZZ_p>& b) {
        Vec<ZZ_p> c;
        sub(c, a, b);
        return c;
    });

    mod.method("VecZZ_p_negate", [](const Vec<ZZ_p>& v) {
        Vec<ZZ_p> result;
        negate(result, v);
        return result;
    });

    mod.method("VecZZ_p_InnerProduct", [](const Vec<ZZ_p>& a, const Vec<ZZ_p>& b) {
        ZZ_p result;
        InnerProduct(result, a, b);
        return result;
    });

    mod.method("VecZZ_p_to_string", [](const Vec<ZZ_p>& v) {
        std::ostringstream oss;
        oss << v;
        return oss.str();
    });

    // =========================================================================
    // GF2 Type (Binary Field Element)
    // =========================================================================

    mod.add_type<GF2>("GF2")
        .constructor<>()
        .constructor<long>()
        .method("__copy__", [](const GF2& a) { return GF2(a); });

    mod.method("GF2_IsZero", [](const GF2& a) { return IsZero(a); });
    mod.method("GF2_IsOne", [](const GF2& a) { return IsOne(a); });
    // GF2 is 0 or 1; 1 is odd, 0 is even
    mod.method("GF2_IsOdd", [](const GF2& a) { return IsOne(a); });

    mod.method("GF2_add", [](const GF2& a, const GF2& b) { return a + b; });
    mod.method("GF2_sub", [](const GF2& a, const GF2& b) { return a - b; });
    mod.method("GF2_mul", [](const GF2& a, const GF2& b) { return a * b; });
    mod.method("GF2_negate", [](const GF2& a) { return -a; });

    mod.method("GF2_IsEqual", [](const GF2& a, const GF2& b) { return a == b; });

    mod.method("GF2_rep", [](const GF2& a) {
        return IsOne(a) ? 1L : 0L;
    });

    // =========================================================================
    // GF2X Type (Polynomials over GF(2))
    // =========================================================================

    mod.add_type<GF2X>("GF2X")
        .constructor<>()
        .constructor<long>()
        .constructor<const GF2&>()
        .method("__copy__", [](const GF2X& f) { return GF2X(f); });

    mod.method("GF2X_deg", [](const GF2X& f) { return deg(f); });

    mod.method("GF2X_coeff", [](const GF2X& f, long i) {
        return coeff(f, i);
    });

    mod.method("GF2X_setcoeff", [](GF2X& f, long i, const GF2& c) {
        SetCoeff(f, i, IsOne(c) ? 1 : 0);
    });

    mod.method("GF2X_setcoeff_long", [](GF2X& f, long i, long c) {
        SetCoeff(f, i, c);
    });

    mod.method("GF2X_leadcoeff", [](const GF2X& f) { return LeadCoeff(f); });
    mod.method("GF2X_constterm", [](const GF2X& f) { return ConstTerm(f); });

    mod.method("GF2X_add", [](const GF2X& f, const GF2X& g) { return f + g; });
    mod.method("GF2X_sub", [](const GF2X& f, const GF2X& g) { return f - g; });
    mod.method("GF2X_mul", [](const GF2X& f, const GF2X& g) { return f * g; });
    mod.method("GF2X_negate", [](const GF2X& f) { return -f; });

    mod.method("GF2X_div", [](const GF2X& f, const GF2X& g) {
        if (IsZero(g)) throw std::domain_error("Division by zero polynomial");
        return f / g;
    });

    mod.method("GF2X_rem", [](const GF2X& f, const GF2X& g) {
        if (IsZero(g)) throw std::domain_error("Division by zero polynomial");
        return f % g;
    });

    mod.method("GF2X_divrem", [](const GF2X& f, const GF2X& g) {
        if (IsZero(g)) throw std::domain_error("Division by zero polynomial");
        GF2X q, r;
        DivRem(q, r, f, g);
        return std::make_tuple(q, r);
    });

    mod.method("GF2X_gcd", [](const GF2X& f, const GF2X& g) {
        return GCD(f, g);
    });

    mod.method("GF2X_diff", [](const GF2X& f) {
        GF2X df;
        diff(df, f);
        return df;
    });

    // Evaluate GF2X polynomial at a GF2 point
    // At x=0: return constant term; at x=1: return XOR of all coefficients
    mod.method("GF2X_eval", [](const GF2X& f, const GF2& x) {
        if (IsZero(x)) {
            return ConstTerm(f);
        } else {
            // x=1: sum all coefficients (XOR in GF(2))
            GF2 result;
            clear(result);
            for (long i = 0; i <= deg(f); i++) {
                result += coeff(f, i);
            }
            return result;
        }
    });

    // For GF2X, use IterIrredTest (DetIrredTest is for ZZ_pX/zz_pX)
    mod.method("GF2X_IterIrredTest", [](const GF2X& f) {
        return IterIrredTest(f);
    });

    mod.method("GF2X_iszero", [](const GF2X& f) { return IsZero(f); });
    mod.method("GF2X_isone", [](const GF2X& f) { return IsOne(f); });

    mod.method("GF2X_to_string", [](const GF2X& f) {
        std::ostringstream oss;
        oss << f;
        return oss.str();
    });

    // =========================================================================
    // vec_GF2 Type (Vector over GF(2))
    // =========================================================================

    mod.add_type<vec_GF2>("VecGF2")
        .constructor<>()
        .method("__copy__", [](const vec_GF2& v) { return vec_GF2(v); });

    mod.method("VecGF2_length", [](const vec_GF2& v) { return v.length(); });

    mod.method("VecGF2_getindex", [](const vec_GF2& v, long i) {
        if (i < 1 || i > v.length()) {
            throw std::out_of_range("Index out of bounds");
        }
        return IsOne(v(i)) ? 1L : 0L;
    });

    mod.method("VecGF2_setindex!", [](vec_GF2& v, long i, long x) {
        if (i < 1 || i > v.length()) {
            throw std::out_of_range("Index out of bounds");
        }
        v(i) = x;
    });

    mod.method("VecGF2_setlength!", [](vec_GF2& v, long n) {
        v.SetLength(n);
    });

    mod.method("VecGF2_add", [](const vec_GF2& a, const vec_GF2& b) {
        vec_GF2 c;
        add(c, a, b);
        return c;
    });

    mod.method("VecGF2_InnerProduct", [](const vec_GF2& a, const vec_GF2& b) {
        GF2 result;
        InnerProduct(result, a, b);
        return IsOne(result) ? 1L : 0L;
    });

    mod.method("VecGF2_to_string", [](const vec_GF2& v) {
        std::ostringstream oss;
        oss << v;
        return oss.str();
    });

    // =========================================================================
    // mat_GF2 Type (Matrix over GF(2))
    // =========================================================================

    mod.add_type<mat_GF2>("MatGF2")
        .constructor<>()
        .method("__copy__", [](const mat_GF2& m) { return mat_GF2(m); });

    mod.method("MatGF2_nrows", [](const mat_GF2& m) { return m.NumRows(); });
    mod.method("MatGF2_ncols", [](const mat_GF2& m) { return m.NumCols(); });

    mod.method("MatGF2_getindex", [](const mat_GF2& m, long i, long j) {
        if (i < 1 || i > m.NumRows() || j < 1 || j > m.NumCols()) {
            throw std::out_of_range("Index out of bounds");
        }
        return IsOne(m(i, j)) ? 1L : 0L;
    });

    mod.method("MatGF2_setindex!", [](mat_GF2& m, long i, long j, long x) {
        if (i < 1 || i > m.NumRows() || j < 1 || j > m.NumCols()) {
            throw std::out_of_range("Index out of bounds");
        }
        m(i, j) = x;
    });

    mod.method("MatGF2_setdims!", [](mat_GF2& m, long n, long k) {
        m.SetDims(n, k);
    });

    mod.method("MatGF2_mul", [](const mat_GF2& a, const mat_GF2& b) {
        mat_GF2 c;
        mul(c, a, b);
        return c;
    });

    mod.method("MatGF2_add", [](const mat_GF2& a, const mat_GF2& b) {
        mat_GF2 c;
        add(c, a, b);
        return c;
    });

    mod.method("MatGF2_gauss", [](mat_GF2& m) {
        return gauss(m);
    });

    mod.method("MatGF2_transpose", [](const mat_GF2& m) {
        mat_GF2 t;
        transpose(t, m);
        return t;
    });

    mod.method("MatGF2_to_string", [](const mat_GF2& m) {
        std::ostringstream oss;
        oss << m;
        return oss.str();
    });

    // =========================================================================
    // zz_p Type (Small Prime Modular Integers - Single Precision)
    // =========================================================================

    mod.add_type<zz_p>("zz_p")
        .constructor<>()
        .constructor<long>()
        .method("__copy__", [](const zz_p& a) { return zz_p(a); });

    mod.method("zz_p_init", [](long p) {
        if (p <= 1) throw std::domain_error("Modulus must be > 1");
        zz_p::init(p);
    });

    mod.method("zz_p_FFTInit", [](long i) {
        zz_p::FFTInit(i);
    });

    mod.method("zz_p_modulus", []() {
        return zz_p::modulus();
    });

    mod.method("zz_p_rep", [](const zz_p& a) {
        return rep(a);
    });

    mod.method("zz_p_add", [](const zz_p& a, const zz_p& b) { return a + b; });
    mod.method("zz_p_sub", [](const zz_p& a, const zz_p& b) { return a - b; });
    mod.method("zz_p_mul", [](const zz_p& a, const zz_p& b) { return a * b; });
    mod.method("zz_p_negate", [](const zz_p& a) { return -a; });

    mod.method("zz_p_inv", [](const zz_p& a) {
        if (IsZero(a)) throw std::domain_error("Inverse of zero");
        return inv(a);
    });

    mod.method("zz_p_div", [](const zz_p& a, const zz_p& b) {
        if (IsZero(b)) throw std::domain_error("Division by zero");
        return a / b;
    });

    mod.method("zz_p_power", [](const zz_p& a, long e) {
        return power(a, e);
    });

    mod.method("zz_p_iszero", [](const zz_p& a) { return IsZero(a); });
    mod.method("zz_p_isone", [](const zz_p& a) { return IsOne(a); });

    // zz_pContext
    mod.add_type<zz_pContext>("zz_pContext")
        .constructor<>()
        .constructor<long>();

    mod.method("zz_pContext_save", [](zz_pContext& ctx) {
        ctx.save();
    });

    mod.method("zz_pContext_restore", [](const zz_pContext& ctx) {
        ctx.restore();
    });

    // =========================================================================
    // zz_pX Type (Polynomials over zz_p)
    // =========================================================================

    mod.add_type<zz_pX>("zz_pX")
        .constructor<>()
        .constructor<long>()
        .constructor<const zz_p&>()
        .method("__copy__", [](const zz_pX& f) { return zz_pX(f); });

    mod.method("zz_pX_deg", [](const zz_pX& f) { return deg(f); });

    mod.method("zz_pX_coeff", [](const zz_pX& f, long i) {
        return coeff(f, i);
    });

    mod.method("zz_pX_setcoeff", [](zz_pX& f, long i, const zz_p& c) {
        SetCoeff(f, i, c);
    });

    mod.method("zz_pX_setcoeff_long", [](zz_pX& f, long i, long c) {
        SetCoeff(f, i, c);
    });

    mod.method("zz_pX_add", [](const zz_pX& f, const zz_pX& g) { return f + g; });
    mod.method("zz_pX_sub", [](const zz_pX& f, const zz_pX& g) { return f - g; });
    mod.method("zz_pX_mul", [](const zz_pX& f, const zz_pX& g) { return f * g; });
    mod.method("zz_pX_negate", [](const zz_pX& f) { return -f; });

    mod.method("zz_pX_div", [](const zz_pX& f, const zz_pX& g) {
        if (IsZero(g)) throw std::domain_error("Division by zero polynomial");
        return f / g;
    });

    mod.method("zz_pX_rem", [](const zz_pX& f, const zz_pX& g) {
        if (IsZero(g)) throw std::domain_error("Division by zero polynomial");
        return f % g;
    });

    mod.method("zz_pX_gcd", [](const zz_pX& f, const zz_pX& g) {
        return GCD(f, g);
    });

    mod.method("zz_pX_DetIrredTest", [](const zz_pX& f) {
        return DetIrredTest(f);
    });

    mod.method("zz_pX_iszero", [](const zz_pX& f) { return IsZero(f); });

    mod.method("zz_pX_to_string", [](const zz_pX& f) {
        std::ostringstream oss;
        oss << f;
        return oss.str();
    });

    // =========================================================================
    // ZZ_pE Type (Extension Field Elements)
    // =========================================================================

    mod.add_type<ZZ_pE>("ZZ_pE")
        .constructor<>()
        .constructor<long>()
        .constructor<const ZZ_p&>()
        .method("__copy__", [](const ZZ_pE& a) { return ZZ_pE(a); });

    mod.method("ZZ_pE_init", [](const ZZ_pX& P) {
        ZZ_pE::init(P);
    });

    mod.method("ZZ_pE_degree", []() {
        return ZZ_pE::degree();
    });

    mod.method("ZZ_pE_modulus", []() {
        return ZZ_pE::modulus();
    });

    mod.method("ZZ_pE_rep", [](const ZZ_pE& a) {
        return rep(a);
    });

    mod.method("ZZ_pE_add", [](const ZZ_pE& a, const ZZ_pE& b) { return a + b; });
    mod.method("ZZ_pE_sub", [](const ZZ_pE& a, const ZZ_pE& b) { return a - b; });
    mod.method("ZZ_pE_mul", [](const ZZ_pE& a, const ZZ_pE& b) { return a * b; });
    mod.method("ZZ_pE_negate", [](const ZZ_pE& a) { return -a; });

    mod.method("ZZ_pE_inv", [](const ZZ_pE& a) {
        if (IsZero(a)) throw std::domain_error("Inverse of zero");
        return inv(a);
    });

    mod.method("ZZ_pE_div", [](const ZZ_pE& a, const ZZ_pE& b) {
        if (IsZero(b)) throw std::domain_error("Division by zero");
        return a / b;
    });

    mod.method("ZZ_pE_power", [](const ZZ_pE& a, long e) {
        return power(a, e);
    });

    mod.method("ZZ_pE_power_ZZ", [](const ZZ_pE& a, const ZZ& e) {
        return power(a, e);
    });

    mod.method("ZZ_pE_iszero", [](const ZZ_pE& a) { return IsZero(a); });
    mod.method("ZZ_pE_isone", [](const ZZ_pE& a) { return IsOne(a); });

    mod.method("ZZ_pE_random", []() {
        return random_ZZ_pE();
    });

    // ZZ_pEContext
    mod.add_type<ZZ_pEContext>("ZZ_pEContext")
        .constructor<>();

    mod.method("ZZ_pEContext_save", [](ZZ_pEContext& ctx) {
        ctx.save();
    });

    mod.method("ZZ_pEContext_restore", [](const ZZ_pEContext& ctx) {
        ctx.restore();
    });

    // =========================================================================
    // ZZ_pEX Type (Polynomials over Extension Field)
    // =========================================================================

    mod.add_type<ZZ_pEX>("ZZ_pEX")
        .constructor<>()
        .constructor<long>()
        .constructor<const ZZ_pE&>()
        .method("__copy__", [](const ZZ_pEX& f) { return ZZ_pEX(f); });

    mod.method("ZZ_pEX_deg", [](const ZZ_pEX& f) { return deg(f); });

    mod.method("ZZ_pEX_coeff", [](const ZZ_pEX& f, long i) {
        return coeff(f, i);
    });

    mod.method("ZZ_pEX_setcoeff", [](ZZ_pEX& f, long i, const ZZ_pE& c) {
        SetCoeff(f, i, c);
    });

    mod.method("ZZ_pEX_add", [](const ZZ_pEX& f, const ZZ_pEX& g) { return f + g; });
    mod.method("ZZ_pEX_sub", [](const ZZ_pEX& f, const ZZ_pEX& g) { return f - g; });
    mod.method("ZZ_pEX_mul", [](const ZZ_pEX& f, const ZZ_pEX& g) { return f * g; });
    mod.method("ZZ_pEX_negate", [](const ZZ_pEX& f) { return -f; });

    mod.method("ZZ_pEX_div", [](const ZZ_pEX& f, const ZZ_pEX& g) {
        if (IsZero(g)) throw std::domain_error("Division by zero polynomial");
        return f / g;
    });

    mod.method("ZZ_pEX_rem", [](const ZZ_pEX& f, const ZZ_pEX& g) {
        if (IsZero(g)) throw std::domain_error("Division by zero polynomial");
        return f % g;
    });

    mod.method("ZZ_pEX_gcd", [](const ZZ_pEX& f, const ZZ_pEX& g) {
        return GCD(f, g);
    });

    mod.method("ZZ_pEX_iszero", [](const ZZ_pEX& f) { return IsZero(f); });

    mod.method("ZZ_pEX_random", [](long n) {
        return random_ZZ_pEX(n);
    });

    mod.method("ZZ_pEX_to_string", [](const ZZ_pEX& f) {
        std::ostringstream oss;
        oss << f;
        return oss.str();
    });

    // Extension field polynomial functions
    mod.method("ZZ_pEX_MinPolyMod", [](const ZZ_pEX& g, const ZZ_pEX& f) {
        ZZ_pEX h;
        MinPolyMod(h, g, f);
        return h;
    });

    mod.method("ZZ_pEX_CompMod", [](const ZZ_pEX& g, const ZZ_pEX& h, const ZZ_pEX& f) {
        ZZ_pEX result;
        CompMod(result, g, h, f);
        return result;
    });

    // =========================================================================
    // RR Type (Arbitrary Precision Floating Point)
    // =========================================================================

    mod.add_type<RR>("RR")
        .constructor<>()
        .constructor<double>()
        .method("__copy__", [](const RR& a) { return RR(a); });

    mod.method("RR_from_ZZ", [](const ZZ& z) {
        return to_RR(z);
    });

    mod.method("RR_from_string", [](const std::string& s) {
        RR r;
        std::istringstream iss(s);
        iss >> r;
        if (iss.fail()) {
            throw std::invalid_argument("Invalid RR string: " + s);
        }
        return r;
    });

    mod.method("RR_to_string", [](const RR& r) {
        std::ostringstream oss;
        oss << r;
        return oss.str();
    });

    mod.method("RR_SetPrecision", [](long p) {
        RR::SetPrecision(p);
    });

    mod.method("RR_precision", []() {
        return RR::precision();
    });

    mod.method("RR_SetOutputPrecision", [](long p) {
        RR::SetOutputPrecision(p);
    });

    mod.method("RR_OutputPrecision", []() {
        return RR::OutputPrecision();
    });

    mod.method("RR_add", [](const RR& a, const RR& b) { return a + b; });
    mod.method("RR_sub", [](const RR& a, const RR& b) { return a - b; });
    mod.method("RR_mul", [](const RR& a, const RR& b) { return a * b; });
    mod.method("RR_div", [](const RR& a, const RR& b) {
        if (IsZero(b)) throw std::domain_error("Division by zero");
        return a / b;
    });
    mod.method("RR_negate", [](const RR& a) { return -a; });
    mod.method("RR_abs", [](const RR& a) { return abs(a); });

    mod.method("RR_sqrt", [](const RR& a) {
        if (a < 0) throw std::domain_error("Square root of negative number");
        return sqrt(a);
    });

    mod.method("RR_exp", [](const RR& a) { return exp(a); });
    mod.method("RR_log", [](const RR& a) {
        if (a <= 0) throw std::domain_error("Logarithm of non-positive number");
        return log(a);
    });

    mod.method("RR_sin", [](const RR& a) { return sin(a); });
    mod.method("RR_cos", [](const RR& a) { return cos(a); });

    mod.method("RR_power", [](const RR& a, long e) {
        return power(a, e);
    });

    mod.method("RR_power_RR", [](const RR& a, const RR& e) {
        return pow(a, e);
    });

    mod.method("RR_iszero", [](const RR& a) { return IsZero(a); });
    mod.method("RR_isone", [](const RR& a) { return IsOne(a); });

    mod.method("RR_compare", [](const RR& a, const RR& b) {
        return compare(a, b);
    });

    mod.method("RR_less", [](const RR& a, const RR& b) { return a < b; });
    mod.method("RR_lesseq", [](const RR& a, const RR& b) { return a <= b; });
    mod.method("RR_equal", [](const RR& a, const RR& b) { return a == b; });

    // Pi constant
    mod.method("RR_ComputePi", []() {
        return ComputePi_RR();
    });
}
