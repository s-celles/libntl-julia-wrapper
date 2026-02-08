# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- CI build pipeline using GitHub Actions
  - Automatic build verification on push to any branch
  - Build status checks on pull requests
  - Multi-platform matrix: Ubuntu (GCC, Clang), macOS (Clang), Windows (MSVC)
  - Compiler caching via ccache for faster rebuilds
  - Julia package caching for faster dependency resolution
- Windows CI support
  - NTL and GMP installation via vcpkg
  - Visual Studio 2022 (MSVC) compiler
  - vcpkg toolchain integration for CMake
