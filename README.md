# Numerical

A numerical library for the GEOframe framework, supporting **hydrological modeling workflows**.

---

## Metadata
- **Language**: Java
- **Build system**: Gradle
- **@author**: Niccolo' Tubini
- **@license**: See the `LICENSE` file for license information (  GNU GENERAL PUBLIC LICENSE Version 3)

---

## Overview

This is a Gradle-based Java project that collects implementations of numerical algorithms used across the GEOframe framework.
The library is intended to provide reusable, general-purpose numerical methods to support scientific and hydrological modeling workflows.

It is designed as a **library module**, not as a standalone application, and is consumed by other GEOframe components.

---

## Numerical Methods

The `it.geoframe.blogspot.numerical` package includes implementations of:

- **Linear system solvers**
  - Thomas algorithm
  - Conjugate Gradient method

- **Nonlinear solvers**
  - Nested Newton methods

- **Root-finding algorithms**
  - Bisection method

- **ODE utilities**
  - Utilities and examples related to solid and liquid snow mass dynamics


---

## Project Structure

Source code are placed in the folder **src/main/java/it/geoframe/blogspot/numerical**.

The package structure reflects the organization of numerical methods by functionality (e.g., solvers, utilities).

---

## Build and Test

---

## Scope and Design Notes


---

