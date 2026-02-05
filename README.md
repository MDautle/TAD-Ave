# intraTAD: Quantifying intra-TAD interactions from Hi-C Data

If using this code please cite : **PAPER WILL GO HERE**


### Overview
This program (I call it TADbit, but call it what you like) is a Python command-line tool that compares intra-TAD contacts between two conditions (e.g. a control and treatment). Using the control matrix, TADbit identifies TAD boundaries via insulation score analysis, then quanitifies the average interactions within each TAD for both conditions. This program takes Hi-C matrices (in .mcool format) and outputs the following: 

#### The Output
+ TAD boundaries
+ Summary statistics per TAD (mean and median values)
+ Statistical comparison between conditions (Wilcoxon rank-sum)
+ Violin plot visualization

##### File names and Descriptions 
| File | Description| 
|---|---|
|```matrix1_TADboundaries.csv``` | Control TAD boundary calls | 
|```matrix1_TADs_df.tsv``` | TADs identified from boundary file | 
|```matrix1/2_statistics.tsv``` | Per TAD statistics for each condition | 
|```intraTAD_intensity_rawData.tsv``` | Mean intra-TAD contact values | 
|```intraTAD_violins.pdf``` | Violin plot comparing intra-TAD intensity | 
|```statistics.csv``` | Wilcoxon rank sum test results | 

### Requirements 
A conda environment file is available for dependency installation. If it doesn't work, the following packages are required: 
+ pandas
+ numpy
+ cooler
+ cooltools
+ matplotlib
+ seaborn
+ scikit-learn

Start with Python version 3.11 or lower, many of the cooler/cooltools packages do not support more recent versions of Python. 

### Usage 
```
python intraTAD_analysis.py \
  -hic1 path/to/control.mcool \
  -hic2 path/to/treatment.mcool \
  -o results_directory/
  -r 20000
```
#### Arguments 
| Flag | Name | Description | 
|---|---|---|
| ```-hic1```/```--matrix1``` | Control Hi-C cooler file | TAD calling reference | 
| ```-hic2```/```--matrix2``` | Treatment/Condition Hi-C cooler file | Comparison dataset | 
| ```-o```/```--output``` | Output directory | Location to write results to | 
| ```-r```/```--resolution``` | Bin size | Defauly 20,000 (You'll need to make sure the file contains resolutions for 3x, 5x, and 10x resolutions) | 




