# libntl-julia-wrapper

CxxWrap-based Julia bindings for [NTL](https://libntl.org/) (Number Theory Library).

## Features

- **ZZ**: Arbitrary-precision integers with arithmetic, GCD, and comparison operations
- **ZZ_p**: Integers modulo p with context management
- **ZZX**: Polynomials over Z with arithmetic and polynomial operations
- **Vec<ZZ>**: Vectors of arbitrary-precision integers
- **Mat<ZZ>**: Matrices of arbitrary-precision integers
- **PrimeSeq**: Prime number iterator
- **Number Theory Functions**: PowerMod, ProbPrime, RandomBnd, RandomBits, etc.

## Requirements

- CMake >= 3.10
- C++17 compiler
- [NTL](https://libntl.org/)
- [CxxWrap.jl](https://github.com/JuliaInterop/CxxWrap.jl) / libcxxwrap-julia

## Building

```bash
mkdir build && cd build
cmake ..
make
```

## License

MIT License - see [LICENSE](LICENSE) for details.
