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


**Code**: 

*Graph Input*:


*Usage*:

```
python main.py --input <input_graph_file_path> --cat <input_cat_file_path> --output <output_embedding_file_path> \
	--dim <embedding_dimension> --scope <max_temporal_distance> --base <constant of logarithm binning> \
	--walk_num <#walks_per_node> --walk_length <walk_length> --walk_mod <temporal_random_walk_bias> --ignore_time <if_ignore_time>

```