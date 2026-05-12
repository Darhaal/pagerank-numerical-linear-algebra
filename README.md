# PageRank and Numerical Linear Algebra

Numerical analysis project implementing the PageRank algorithm using sparse matrix methods, power iteration, dangling-node correction, and direct-solve verification in MATLAB.

This project was completed at Florida Institute of Technology for MTH 4311 – Numerical Analysis.

## Repository Description

Numerical analysis project implementing PageRank with sparse power iteration, dangling-node correction, and direct-solve verification in MATLAB.

## Project Overview

This project studies PageRank as an application of numerical linear algebra to search engine systems.

PageRank is usually introduced as a web ranking algorithm, but mathematically it becomes:

- a fixed-point problem
- an eigenvector problem
- a sparse linear algebra problem
- an iterative numerical method problem

The web is represented as a directed graph. Web pages are nodes, and hyperlinks are directed edges. From this graph, the project constructs an adjacency matrix, transition matrix, and PageRank vector.

The main numerical method is sparse power iteration.

## Main Question

How can PageRank be formulated and solved as a numerical linear algebra problem?

## Mathematical Idea

If page j links to page i, then the adjacency matrix entry is:

A(i,j) = 1

This column convention means each column describes outgoing links from one page.

For a page with outgoing degree d_j, the transition probability is:

S(i,j) = A(i,j) / d_j

The PageRank vector satisfies the fixed-point equation:

r = alpha * S_full * r + (1 - alpha) * v

Equivalently, it can be written as a linear system:

(I - alpha * S_full) r = (1 - alpha) v

This shows that PageRank is not only a ranking idea. It is also a numerical linear algebra problem.

## Numerical Method

The main method used in this project is power iteration.

Starting from a uniform probability vector, the algorithm repeatedly applies the PageRank update until the L1 change between two consecutive vectors is below the stopping tolerance.

Stopping criterion:

||r^(k+1) - r^(k)||_1 < 10^-10

The PageRank vector is interpreted as the dominant eigenvector of the Google matrix.

## Sparse Matrix Implementation

The implementation uses sparse matrix operations in MATLAB.

Instead of building the transition matrix by modifying sparse columns one by one, the transition matrix is built using sparse diagonal scaling:

S = A * D^-1

This is important because MATLAB sparse matrices use compressed sparse column storage. Index-vector assembly and diagonal scaling avoid unnecessary memory reallocation.

## Dangling Node Correction

A dangling node is a page with no outgoing links.

In this project, page 10 is a dangling node.

Instead of filling the dangling column with a dense uniform vector, the algorithm handles dangling nodes implicitly.

At each iteration, the total PageRank mass on dangling pages is computed and redistributed through the teleportation vector.

This preserves sparse structure and keeps the PageRank vector as a probability distribution.

## Experiment

The numerical experiment uses a small directed web graph.

Main setup:

- Number of pages: 10
- Number of directed links: 23
- Dangling node: page 10
- Damping factor: alpha = 0.85
- Stopping tolerance: 10^-10
- Maximum iterations: 1000

The graph is intentionally small so that the algorithm and numerical behavior can be inspected clearly.

## Results

Sparse power iteration converged in 27 iterations.

Main results:

- Final L1 change: 9.474 x 10^-11
- Sum of PageRank vector: 1.000000000000
- Maximum mass error before renormalization: approximately machine precision
- Direct solve L1 difference: 5.539 x 10^-11

The direct solve was used only as verification for the small graph. It confirms that the power iteration result solves the same fixed-point problem up to numerical tolerance.

## Highest PageRank Values

The highest-ranked pages in the experiment were:

- Page 3: 0.15920
- Page 6: 0.13819
- Page 7: 0.11572
- Page 10: 0.11522
- Page 1: 0.10269

Page 3 has the highest PageRank value because it receives several important incoming links.

Page 10 is dangling, but it still has high PageRank because multiple pages point to it. This shows that PageRank is not the same as simply counting outgoing links or incoming links.

## Damping Factor Comparison

The project also compares different damping factors.

Results:

alpha = 0.70:
- Iterations: 22

alpha = 0.85:
- Iterations: 27

alpha = 0.95:
- Iterations: 32

As alpha increases, convergence becomes slower. Numerically, this is connected to the spectral gap and the behavior of subdominant eigenvalues.

## Numerical Analysis Concepts

This project connects to several topics from Numerical Analysis:

- eigenvalue problems
- dominant eigenvectors
- power iteration
- fixed-point iteration
- sparse matrices
- stochastic matrices
- stopping tolerance
- convergence rate
- floating-point round-off
- probability mass conservation

## Files in This Repository

All project files are stored in the root of the repository.

Expected files:

- README.md
- LICENSE
- NOTICE.md
- MTH_4311_Final_Project_Code.m
- MTH4311_Project2_Report.pdf
- MTH4311_Project2_Slides.pdf

The MATLAB file contains the implementation.

The PDF report contains the full explanation, mathematical formulation, numerical results, and appendix code.

The slides summarize the project for presentation.

## How to View the Project

Open the report PDF to read the full project explanation.

Open the slides PDF to review the presentation version.

Open the MATLAB file to inspect or run the implementation.

The MATLAB script saves figures for:

- directed web graph
- computed PageRank values
- power iteration convergence
- damping factor comparison
- probability mass error

## How to Run

Open the MATLAB file:

MTH_4311_Final_Project_Code.m

Run the script in MATLAB.

The script will:

1. Build the sparse adjacency matrix.
2. Construct the sparse transition matrix.
3. Identify dangling nodes.
4. Run sparse power iteration.
5. Print PageRank results.
6. Generate plots.
7. Compare damping factors.
8. Verify the result using direct solve.
9. Check probability mass conservation.

## Reproducibility Notes

This repository is mainly shared as an academic and portfolio project.

The project does not require external datasets.

The graph is hard-coded inside the MATLAB script so that the numerical experiment can be reproduced directly.

MATLAB is required to run the code.

## Important Notes

This project is not trying to build a real search engine.

The goal is to show how a search-ranking problem becomes a numerical linear algebra problem.

The small graph is used for clarity and verification. Real PageRank systems use much larger sparse graphs, where direct solve is not practical and iterative sparse methods are necessary.

## Limitations

There are several limitations:

- The graph has only 10 pages.
- Runtime scaling is not fully demonstrated.
- The direct solve check is only practical for small graphs.
- The graph is synthetic and does not represent a real web crawl.
- The project focuses on numerical formulation and correctness, not production-scale ranking.

## Main Finding

PageRank is a strong example of applied numerical linear algebra.

A search engine ranking problem becomes an eigenvector and fixed-point problem, and the solution depends on sparse matrix computation, iterative methods, convergence behavior, and numerical stability.

## License

This project is licensed under the MIT License.
