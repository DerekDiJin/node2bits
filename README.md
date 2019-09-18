# node2bits

**Paper**: Di Jin, Mark Heimann, Ryan A. Rossi, Danai Koutra. node2bits: Compact Time- and Attribute-aware Node Representations for User Stitching. ACM SIGKDD International Conference on Knowledge Discovery and Data Mining (KDD), 2019.

<p align="center">
<img src="https://derekdijin.github.io/assets/projects/node2bits_overview_final.jpg" width="550"  alt="Overview of node2bits">
</p>


**Citation (bibtex)**:
```
@inproceedings{node2bits-ECML19,
   author={Di Jin and Mark Heimann and Ryan A. Rossi and Danai Koutra},
   title={Node2BITS: Compact Time- and Attribute-aware Node Representations for User Stitching},
   booktitle={ECML/PKDD},
   year={2019},
   pages={22},
}

Di Jin, Mark Heimann, Ryan A. Rossi, and Danai Koutra. "Node2BITS: Compact Time- and Attribute-aware Node Representations for User Stitching." ECML/PKDD, pp. 22. 2019.
```


# Code

## Inputs:

node2bits takes two files as input, the graph file and the category file.

### Input graph file
The input graph file can be either static or temporal edge list in the following format separated by tab:
```
<src> <dst> <weight> <timestamp> (optional)
```
node2bits will automatically determine if the input graph is static or temporal. The edge list is assumed to be re-ordered consecutively from 0, i.e., the minimum node ID is 0, and the maximum node ID is <#node - 1>. A toy static graph is under "/graph/" directory.

### Input category file
The category file is a mapping between the node ID and its type (e.g., IP, cookie, web agent) with the following format separated by tab:
```
<category> <id_initial> <id_ending>
```
if the node IDs are grouped by the type, where ```<id_initial>``` and ```<id_ending>``` are the starting and ending node ids in type ```<category>```
For example,
```
0	0	279629
1	279630	283182
```
means node 0, 1, ... 279629 are in type 0, node 279630, 279631, ... 283182 are in type 1.

But if the node IDs are not grouped by the types, this implementation also supports the following format separated by tab:
```
<category> <node_id>
```
which is just the 1-1 mapping. The code accepts either format.

## Inputs

The complete command to run node2bits is as follows.

```
python main.py --input <graph_file_path> --cat <category_file_path> --output <embedding_file_path> --dim <embedding_dimension> 
	--scope <max_temporal_distance> --base <constant of logarithm binning> --walk_num <#walks_per_node> --walk_length <walk_length> 
	--walk_mod <temporal_random_walk_bias> --ignore_time <if_ignore_time>

```

- input, the input graph file stated under the "Graph Input" section above.
- cat, 