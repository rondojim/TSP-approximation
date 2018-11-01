# Travelling-Salesman-Problem

The Traveling Salesman Problem (TSP) is one of the most fundamental algorithmic problems in Computer Science. Existing exact algorithms (ie., algorithms that produce the optimal solution) have exponential worst-case time-complexity. Due to the difficulty of TSP, approximation algorithms have been proposed that are much more efficient than existing exact algorithms but whose solutions are not necessarily optimal. 

One of the best-known approximation algorithms for TSP is Christofides' algorithm. The purpose of our work is to propose a heuristic for enhancing the accuracy of Christofides' algorithm and to provide an experimental evaluation of its performance. More specifically, we consider a number of cities with integer coordinates in a restricted square of the Euclidean plane (namely in the grid [0,3]x[0,3]). Such a placement of cities may produce, in the general case, many different minimum spanning trees of the complete graph connecting all the cities. We propose a modification of Christofides' algorithm which uses many different minimum spanning trees in order to produce the approximation of the optimal solution. 

We provide an experimental evaluation of our heuristic. In particular, we implement in C++ three algorithms: (a) a brute-force algorithm that returns the optimal TSP solution and has exponential complexity, (b) Christofides' algorithm which has much better complexity but returns an approximate solution, and (c) our extension of Christofides' algorithm that examines many alternative minimum spanning trees. The experimental results suggest that our heuristic returns an approximate solution that is on the average significantly better than the one produced by the original Christofides' algorithm.

Files:

input.txt: contains the 8008 different placements of 10 cities in the grid [0,3]x[0,3].

ratio_comparison.cpp: C++ implementation of a brute-force dynamic programming algorithm for TSP, and of our heuristic for TSP.

results(k=1).txt,...,results(k=10).txt: Each of these files contains the performance of our extension for a different number of minimum spanning trees used. The case k=1 is the original algorithm of Christofides. The case k=10 uses 10 minimum spanning trees and selects the best Hamilton cycle produced by these 10 different trees.

