## TCEM: Two-Phase Co-exposure maximization Algorithm

Contact author: Sijing Tu [sijing@kth.se](mailto:sijing@kth.se) 

- Citation information: Tu, S., Aslay, C., Gionis, A. (2020). Co-exposure Maximization in Online Social Networks. In *Advances in Neural Information Processing Systems (NeurIPS)*

- Copyright: Redistribution and use in source and binary forms, with or without modifications, are permitted for academic purposes, provided that the proper acknowledgements are done.

Acknowledgement: Cigdem Aslay.

### Compilation  

**make -f Makefile**


### Configuration of input files and parameters 

**config.txt:**: specify the output folder, set &epsilon = 0.2, and l = 1. It is possible 
to modify these parameters.

### Input file of the graph 
Input directed graph file with one line per arc followed by two probabilities, the node ids should be mapped to 0 to n-1 where n is the total number of nodes. 

format: node_u node_v p^1_uv p^2_uv 

### compare file
Indicate the selected red nodes ids and blue nodes ids. The file should contain two lines; the first 
line begins with `red nodes: `, and the second line begins with `blue nodes: `. The node ids should separate by commas.  

format:

red nodes: x1 x2 ... xk_r 

blue nodes: y1 y2 ... yk_b 


### Running from command line
There are two possible command types, they have the common part e.g. **./main_TCEM -c config.txt -x data/karate.txt**. The flag `-c` is followed by a configure file, in this project, it is `config.txt`; the flag `-x` is followed by the file of input graph. 

For the first one you need to specify the number of red nodes and blue nodes; the number of red nodes and blue nodes are followed by `-r` and 
`-b` respectively.

**./main_TCEM -c config.txt -x data/karate.txt -r 2 -b 4**

For the second one you need to specify the compare file, which is followed by `-p`.

**./main_TCEM -c config.txt -x data/karate.txt -p comparefolder/compare.txt**


**Note:** The implementation also contains the code for baselines used. Comment out line 88-93 of allocator.cc to receive also the results for baselines. 


### Reference
 *Aslay, C., Galbrun, E., Matakos, A., Gionis, A. (2018). Maximizing the diversity of exposure in a social network. IEEE International Conference on Data Mining (ICDM).* 
 
### Note
Update in Oct 2020,
We made a modification on lambda in our code, changing from 
`theta = (2 + 2/3.0 * epsilon) * n * (ell * log(n) + log(2) +  logP(n,k_r,tao)) / (epsilon * epsilon * lb);`
to 
`theta = (8 + 4/3.0 * epsilon) * n * (ell * log(n) + log(2) +  logP(n,k_r,tao)) / (epsilon * epsilon * lb);`.



