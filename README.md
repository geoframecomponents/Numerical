# Numerical

A numerical library for the GEOframe framework, supporting **hydrological modeling workflows**.

---
## Metadata
- **Language**: Java
- **Build system**: Maven
- **@author**: Niccolo' Tubini, Concetta D'Amato, Daniele Andreis, Giuseppe Formetta
- **@license**: See the `LICENSE` file for license information (  GNU GENERAL PUBLIC LICENSE Version 3)

---
## License

This project is licensed under the GNU General Public License v3.0. See LICENSE.

---

[!WARNING]
 **Repository Update**

 This repository has undergone several changes compared to the original version. In particular, the build system has been migrated from **Gradle** to **Maven**.

 If you are looking for the original Gradle-based implementation or wish to continue development from the previous codebase, you can find it under the Git tag **`v1.0.3`**.


 ## Installation

The library is available on **[Maven Central](https://central.sonatype.com/artifact/org.geoframe/numerical)**

Add the following dependency to your `pom.xml`:

```xml
<dependency>
    <groupId>org.geoframe</groupId>
    <artifactId>numerical</artifactId>
    <version>1.0.4</version>
</dependency>
```

---

## Numerical Methods


The `it.geoframe.blogspot.numerical` packages include:

- **Linear system solvers**
  - Thomas algorithm (tridiagonal systems)
  - Conjugate Gradient method (matrix-free via `Matop`; intended for SPD systems)

- **Nonlinear solvers**
  - Nested Newton methods (Thomas-based and CG-based variants)

- **Root finding**
  - Bisection method

- **ODE utilities**
  - Utilities and examples related to solid and liquid snow mass dynamics


---

## Project Structure

Source code is located under:

`src/main/java/it/geoframe/blogspot/numerical`

Main packages:

- `linearsystemsolver/` — Thomas and Conjugate Gradient solvers  
- `newtonalgorithm/` — Nested Newton implementations  
- `matop/` — matrix-free operator abstraction (`Matop`)  
- `rootfinding/` — bisection method  
- `ode/` — ODE utilities and examples

---
## Design notes

* This repository provides reusable numerical building blocks; domain-specific models live in other GEOframe components.

* Some solvers support matrix-free computations through the Matop abstraction.

* The Thomas algorithm modifies internal working arrays; pass copies if you need to reuse input diagonals.

* Nested Newton and bisection methods depend on GEOframe closure-equation abstractions (e.g., EquationState).


---

## Build and Test

---

## Scope and Design Notes


---

