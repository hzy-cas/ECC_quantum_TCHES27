# Binary-ECC quantum-resource comparison tables

This directory contains five PDF tables covering single-step in-place addition, [JSB+25] point addition, the complete controlled point-addition stage of Shor's algorithm, AC-based inversion resources, and the Windowed QROM architecture.

## 1. Comparison of Three Single-Step In-Place Implementations

File: [PDF](point_inplace_comparison.pdf);

This PDF compares three implementations of single-step binary elliptic curve in-place addition: FLT-in, AC-based multiplication with [JSB+25] inversion, and AC-based multiplication with AC-based inversion. The table shows the resource trade-offs of the three schemes in terms of gate count, width, Tofoli depth, NCT depth, DW-cost, and TDW-cost.

## 2. Single-Point Addition Resources for [JSB+25]

File: [PDF](JSB25_point_addition_comparison_table.pdf);

This PDF presents resource estimates for the FLT-in single-point in-place addition circuit.

## 3. Shor's Algorithm's Full Controlled Point Addition Stage

File: [PDF](<Costs for the controlled point addition stage in Shor’s algorithm.pdf>)

This PDF compares the [JSB+25] In-place and [JSB+25] Out-of-place strategies in Shor's algorithm's full controlled point addition stage with our Balanced, $T$-optimal, and AC-based Balanced strategies. Tables are used to observe the overall trade-offs in gate number, depth, width, and spatiotemporal costs of the full stage.

## 4. AC-Based Inversion Circuit Resources

File: [PDF](Table5_inversion_update.pdf)

## 5. Windowed QROM Architecture for the Double-Scalar-Multiplication Stage

File: [PDF](window.pdf)
