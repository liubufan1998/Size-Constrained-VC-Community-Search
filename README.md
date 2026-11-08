# Vertex Connectivity-based Size-Constrained Community Search

This repository implements **SVC-TD**, the top-down exact algorithm proposed in our paper:

**Maximizing Robustness to Vertex Failures in Size-Constrained Community Search**

Given an undirected graph, a query vertex, and a size range, the goal is to find a connected query-containing community whose size satisfies the given constraint and whose vertex connectivity is maximized.

The implementation also uses **SVC-Heu** to obtain an initial feasible community and the corresponding lower bound before the exact search.

## Technical report

[Read or download the full technical report (PDF)](docs/SVC-Technical-Report.pdf).

The report includes the complete paper, references, proofs, and additional experiments (17 pages).

## Compile the code

```sh
$ make clean
$ make
```

It generates an executable file named `SVC-TD`, which corresponds to our SVC-TD algorithm.

## Get datasets

The 12 datasets used in our paper can be found from the following sources:

1. UK2002 and Webbase are downloaded from WebGraph [1].

2. All the other graphs are downloaded from SNAP [2].

[1] Paolo Boldi and Sebastiano Vigna. 2004. The WebGraph framework I: compression techniques. In Proceedings of the International Conference on World Wide Web. 595–602.

[2] Jure Leskovec and Andrej Krevl. 2014. SNAP Datasets: Stanford Large Network Dataset Collection.

## Data format

Our program accepts a **binary input graph file** with a `.bin` extension.

If the graph is not in binary format, please use the provided data transformation tool `translate.cpp`, which is contained in the directory `data_transform_tool`, to convert the original graph file into binary format.

We also include an example raw graph file `bio-celegans.mtx` in the same directory to illustrate the input format before conversion.

An example of converting a `.mtx` graph file into binary format is shown below:

```sh
g++ -O3 -std=c++11 translate.cpp -o a.out -w
./a.out bio-celegans.mtx bio-celegans.bin
```

The raw graph file before conversion should follow the format below. The two integers in the first line denote the number of vertices and edges, respectively.

```text
4 5
1 2
1 4
2 3
2 4
3 4
```

Each following line represents one undirected edge.

## Run the code

The usage of our program is as follows:

```sh
$ ./SVC-TD -g {path_to_graph} -l {size_lower_bound} -u {size_upper_bound} -q {query_vertex}
```

An example of solving the SVC problem on the graph `Google.bin` is shown below:

```sh
$ ./SVC-TD -g dataset/Google.bin -l 6 -u 9 -q 10
```

This command searches for an optimal community containing query vertex `10` within the size range `[6, 9]`.

## Result Analysis

If you run the example above, you may see outputs similar to the following:

```text
**** SVC-TD (Release) build at 16:28:06 May 29 2026 ***
The size lower bound LB is: 6 ,the upper bound UB is: 9 ,and the query vertex is 10
File: Google n= 875713 m= 4322051
Graph init ok
The value of vc_min identified by our heuristic method is 1
    The total running time is: 980,467 (microseconds)
    The optimal vertex connectivity is: 3
The members of the found optimal community are shown as follow:
10, 158281, 16434, 7, 17, 18,
Everything is OK!!!
```

From the output, we can obtain the following information:

- **Graph Data Characteristics:**
  - **File Name:** `Google`
  - **Vertices (n):** 875,713
  - **Edges (m):** 4,322,051

- **Heuristic Information:**
  - **Initial vertex connectivity found by SVC-Heu:** 1

- **Search Result:**
  - **Total Runtime:** 980,467 microseconds
  - **Optimal vertex connectivity:** 3
  - **Members of the optimal community:** shown in the output