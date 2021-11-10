# BWT + MTF + Huffman compressor

File compression/decompression tool based on BWT (Burrows–Wheeler transform), MTF (Move-to-front transform) and Huffman coding algorithms. 

#### Comprassion steps:

1. Burrows–Wheeler transform | O(N) memory complexity
2. Move-to-front transform | O(N) memory complexity
3. Huffman | naive implementation with binary tree

#### Decomprassion steps:

1. Huffman reverse transform | naive implementation with binary tree
2. Move-to-front reverse transform | O(N) memory complexity
3. Burrows–Wheeler reverse transform | O(N) memory complexity

Test files for benchmark are stored in `cmake-build-release/calgarycorpus` folder. **Each initial and decoded files are equal with a bit precision** (to validate run benchmark (`main.cpp`) and then `check_diff.sh` in the test files folder).

## Compression benchmarks

| file name | initial file size (byte) | encoded file size (byte) | average amount of bits per byte (bit) | meta information size (byte) | compression rate |
| :-------: | :----------------------: | :----------------------: | :-----------------------------------: | :--------------------------: | :--------------: |
|    bib    |          111261          |          33205           |                2.38754                |             162              |     0.298442     |
|   book1   |          768771          |          267163          |                2.78016                |             161              |     0.34752      |
|   book2   |          610856          |          186994          |                2.44894                |             176              |     0.306118     |
|    geo    |          102400          |          69563           |                5.43461                |             344              |     0.679326     |
|   news    |          377109          |          133517          |                2.83243                |             178              |     0.354054     |
|   obj1    |          21504           |          11785           |                4.3843                 |             344              |     0.548038     |
|   obj2    |          246814          |          88733           |                2.87611                |             344              |     0.359514     |
|  paper1   |          53161           |          18224           |                2.74246                |             167              |     0.342808     |
|  paper2   |          82199           |          28136           |                2.73833                |             173              |     0.342291     |
|    pic    |          513216          |          101508          |                1.5823                 |             277              |     0.197788     |
|   progc   |          39611           |          13699           |                2.76671                |             172              |     0.345838     |
|   progl   |          71646           |          18745           |                2.09307                |             171              |     0.261634     |
|   progp   |          49379           |          12826           |                2.07797                |             169              |     0.259746     |
|   trans   |          93695           |          22400           |                1.91259                |             178              |     0.239074     |

Disadvantages of the implementation:

* Loading complete file into the RAM
* Inability to work with folders
* Not the best time complexities (e.g. O(N^2logN) for BWT) of the algorithms (but also not the worst)