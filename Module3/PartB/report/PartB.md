# System biology module 3 part B

## File description



## General Purpose Modules

usage:

```python
from draw_MR import *
from degreeAnalysis import *
```

### draw_MR.py

`read_MR_txt`: read from a txt file, for example, fructoseanaerobicfluxgraph.txt
and return the list of edges `el` and nodes `nl` in the graph.

`draw_MR_Graph`: use the list returned by `read_MR_txt` to generate a `networkx`
graph.

### degreeAnalysis.py

`my_dijkstra`: calculate the shortest distance from a source node to the rest of
the nodes using Dijkstra algorithm. `weight` is a optional argument and it can
be `string`, `int` or `float`.

- `string`: use the value in the field of the string as the weight.

- `int` or `float`: use this value as the fixed weight for all edges.

`degree_stats`: plot a histogram of the number of the nodes with different numbers
of degree. Print the average degree, variance of degree and the hubs node (with
greatest degree).

`draw_bar`: draw a bar chart using the given data. It is called  in `degree_stats`
and `compartment_analysis`.

`compartment_analysis`: plot a histogram of the number of reactions and metabolites
with regard to the compartment they belong to.

## B.2, B.3

Please refer to notebook `Exercise_B2B3.py` or report `Exercise_B2B3.html`.

## B.4

To block the product `M_mthgxl_e`, one has to block all the reactions that produces
it. There is only one reaction in the graph that produce `M_mthgxl_e`, which is
`R_EX_mthgxl(e), ID = R_EX_mthgxl_LPAREN_e_RPAREN_`.

I did the **following** work in this script:

1. add two functions `read_block_list` and `block_reactions`.

    - read_block_list: read the reactions to block from a text file and return the
      reactions as a list

    - block_reactions: set the kineticLaw parameters `UPPER_BOUND` and `LOWER_BOUND`
      of all the reactions in block list to zeros so that the reaction is blocked
      from both sides

2. modify your function `max_flux` to `my_max_flux` so that it calls my functions
described above when calculating the maximum flux. Also modify `max_fluxes` to
`my_max_fluxes` in which `max_flux` is called.

The resulted R/M graph can be plotted using script `draw_MR.py`:

```bash
$ python draw_MR.py exclude_toxic.txt
```

Result of script is stored in file `exclude_toxic.txt`. It can be verified that
there are no product `M_mthgxl_b` any more.

## B.5, B.6

First, in script Exercise_B5B6.py, I generate the text file of different graphs.

- Fructose aerobic: EX_fru(e)True.txt
- Glucose anaerobic: EX_glc(e)False.txt
- Glucose aerobic: EX_glc(e)True.txt
- Fructose + glucose anaerobic: EX_glc(e)EX_fru(e)False.txt
- Fructose + glucose aerobic: EX_glc(e)EX_fru(e)True.txt

### Flux result:

- Fructose aerobic: 1600
- Glucose anaerobic: 100
- Glucose aerobic: 1600
- Fructose + glucose anaerobic: 100
- Fructose + glucose aerobic: 1600

### Plotting
Then, in notebook Exercise_B5B6.ipynb, I repeated all the analysis for these new
graphs. Please read the notebook or report Exercise_B5B6.html for more details.
