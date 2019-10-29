# Integration of Single Cell and Spatial Transcriptomics Data

This repository contains the Python package **stereoscope**, being the implementation of the method presented in the
paper. Additionally scripts used to preprocess, visualize and compare data/results presented in the paper is also found
herewithin.

Below you will find  examples of how to use stereoscope, whilst these are cast as a guide in how to reproduce the
results the procedure can of course be generalized to any data set.

## Installing stereoscope
To perform any of the analyses below you need to install stereoscope, we have also prepared some datasets for you (found
within this github repo) hence the first thing we will do is to clone this repository. Open your terminal and change
to a suitable directory - then enter:

```console
foo@bar:~$ git clone https://github.com/almaan/stereoscope 

```

Now let's start by installing stereoscope, you can easily do this by the following commands
```bash
foo@bar:~$ cd stsc
foo@bar:~$ ./setup install --user

```
This should both give you access to the stereoscope python package (stsc) and provide you with a command line interface
(CLI) meaning you can conduct the analysis from the terminal. To make sure the installation was successfull, we
will run two tests:

```console
foo@bar:~$ python3 -c "import stsc; print(stsc.__version__)"
stereoscope : 0.2.0
foo@bar:~$ stereoscope test
successfully installed stereoscope CLI
```
If you cannot access stereoscope from the command line, and get something like
```console
foo@bar:~$ stereoscope test
bash stereoscope: command not found..

```

make sure that your install location is a part of your PATH variable. Including the ```--user``` flag during installation,
should place your packages in  ~/.local directory. Thus you can likely, 
fix this by entering the follwing into your terminal:


```console
foo@bar:~$ export PATH=$PATH:/user/home/.local/bin
```
(This a line you may want to add to your ~/.bashrc file, or equivalent, for a more seamless usage in the future)

Having installed stereoscope, we can now start with the analysis.


## Reproducing The Mouse Brain Analysis
Let us begin by reproducing the results presented for the mouse brain (hippocampal) region. Here we will go through the
whole workflow from downloading the data, preprocessing, analysis and visualization. If you aren't that keen on doing
the first parts and want to get immediately to the analysis part, you can skip step 1-2 and use the already processed
data found in the data/mousebrain folder, just unzip these files into a folder named data/curated.


### 1. Downloading the data
The data was originally downloaded from [mousebrain.org](http://mousebrain.org/tissues.html), where we downloaded the
loom-file originating from the Hippocampus. You can either download this set via the web browser, or the terminal using

```console
foo@bar:~$ cd data
foo@bar:~$ mkdir raw curated
foo@bar:~$ cd raw
foo@bar:~$ curl -O https://storage.googleapis.com/linnarsson-lab-loom/l1_hippocampus.loom  

```
### 2. Prepare Data

#### 2.1 Subsample Single Cell Data
We will subsample the single cell data, this is not a requirement, but it allows us to run the analysis a bit faster. We will
do this in two steps, first we'll create a modified loom-file where we have added a new column joining the ''Cluster''
and ''Class'' labels together, this adds some interpretability to the cluster labels. Enter the following in the
terminal:

```console
foo@bar:~$ ../../preprocess/hippocampus/create-mod-loom.py l1_hippocampus.loom .
successfully created modifed loom-file >> mod_l1_hippocampus.loom
```
Next we subsample out dataset, using a lower bound and upper bound (See Methods) of 25 repsectively 250 cells per celltype. Do this by
entering the following into the terminal:

```console
foo@bar:~$ ../../preprocess/hippocampus/subsample-data.py -lf mod_l1_hippocampus.loom -o ../curated -lb 25 -ub 250 -cn
"Celltype_1"

Unique Identifier for set >> 20191029082957812352            
Astrocytes_13 | was discarded due to insufficient number of cells
Astrocytes_14 | Used 250 cells                               
Astrocytes_38 | was discarded due to insufficient number of cells
Astrocytes_39 | was discarded due to insufficient number of cells
Astrocytes_40 | Used 250 cells                               
Astrocytes_41 | Used 31 cells
Astrocytes_42 | Used 250 cells
Astrocytes_44 | Used 41 cells
....

```
This will create three files in the data/curated folder - a count matrix of the cells included in the set, a meta-file
containing their respective labels and a ''stats-file'' which displays the composition of the set. All belonging to the
same set are marked with a unique identifier, which is time and date-based, your identifier will therefore be different
from what you see here. As you also can see not all ''cell types'' were included, due to having very few members. 

#### 2.2 Unzip ST-data
We have included the two sections of the mouse brain used in the paper as .tsv files in the repo, there is no need to
preprocess these - but we have zipped them to save some space. Unzip these files and place in the data/curated folder as
well - either interactively or for example entering the following into the terminal:

```console
foo@bar:~$ basename $( pwd )
curated
foo@bar:~$ unzip ../mouse_brain/mouse-st-data.zip -d .
Archive:  ../mouse_brain/mouse-st-data.zip
  inflating: ./st-hippo1.tsv         
  inflating: ./st-hippo2.tsv         
```

If all steps have been successfull, this (below) is the content which you should have within the data/curated folder

```console
foo@bar:~$ ls -1
20191029082957812352.cnt_data.tsv
20191029082957812352.mta_data.tsv
20191029082957812352.stats.tsv
st-hippo1.tsv
st-hippo2.tsv

```
####2.3 Format data
The ouput from the subsampling and the ST data that we provide you with will already be given in the correct format,
hence you will not have to do any additional work here; still we here list what type of input files and what format
that is required in order to run stereoscope.

* **Single Cell Count Data File** - A .tsv file with cells as rows and genes as columns, each cell (row) should have a unique label
* **Single Cell Annotation Data** - A .tsv file with the same rownames as the count data file, either with one single column
  listing the annotations, or multiple columns where the column containing the labels should be named 'bio_celltype'
* **Spatial Transciptomics (ST) count data** - A .tsv file containing with spots as rows and genes as columns

Some additional thing to keep in mind are:
* Make sure that the ST and single cell data uses the same gene identifiers, i.e. that not one is using ENSEBL ids
  whilst the other one has HGNC gene symbols.
* Do not normalize your data - all counts should be given in raw counts




## 3. Analysis

### 3.1 Run the analysis
Now when the data is prepared, we can run the actual analysis. We will run the complete analysis,  estimating rates and
logits from the single cell data and then using these to infer the proportion values in our spatial data. We will use
the following specifics in our analysis

| parameter | values |
| --- | --- |
| number of genes | 5000 |
| sc epochs | 75000 |
| sc batch size | 100 |
| st epochs | 75000
| st batch size | 100|
| gpu | True |

To run the analysis enter the following into your terminal

```console
foo@bar:~$ cd ../../res
foo@bar:~$ stereoscope run --sc_cnt ../data/curated/*cnt*.tsv --sc_labels ../data/curated/*mta*.tsv -sce 75000  -o hippo_1 -n 5000 --st_cnt ../data/curated/st-hippo*tsv -ste 75000 --gpu -stb 100 -scb 100
[2019-10-29 09:07:50,891 - stsc - INFO ] >> Using device cuda
[2019-10-29 09:07:50,891 - stsc - INFO ] >> fitting sc data | count file : ../data/curated/20191029082957812352.cnt_data.tsv | labels file : ../data/curated/20191029082957812352.mta_data.tsv                                                                                 
[2019-10-29 09:09:51,527 - stsc - INFO ] >> SC data GENES : 5000  SC data CELLS : 8449  SC data TYPES : 56

Epoch : 211  /75000 | Loss : 3.004263E+07 | [                     ]
```
This will create a subfolder named "hippo_1" in the res folder, where all results and data related to this analysis will
be found.

For more information regarding which arguments and configurations you can make to your analysis use
```console
foo@bar:~$ stereoscope run -h
```

### 3.2 Insepect progress
Even with a GPU the analysis will take some time to complete, and whilst you have the progress bar show the current
status of the process - it's also of interest to put this into a broader perspective and look at the progress
over time. We can do so by the following command :

```
foo@bar:~$ stereoscope progress -lf hippo_1/sc_loss*txt & disown

```
This will open up a matplotlib interactive window, where you can zoom and move around, being updated every 10 seconds where you can see how the loss has
changed over time. By appending ```& disown``` to the command, allows us to still use the terminal and won't kill the job
if we were to exit the shell, if your shell does not support this command you can use ```nohup``` instead (to be put before
the command rather than after).

![alt text](https://github.com/almaan/stscpaper/blob/master/imgs/lossprog.png?raw=true "Progress")

### 4. Inspect Results

**4.1 Orienting the of results**

When the complete analysis is done your res folder should contain the following content
```console
foo@bar:~$ ls -1
logits.2019-10-29090750.880065.tsv
R.2019-10-29090750.880065.tsv
sc_loss.2019-10-29090750.880065.txt
sc_model.2019-10-29090750.880065.pt
st-hippo1
st-hippo2
st_loss.2019-10-29090750.880065.txt
st_model.2019-10-29090750.880065.pt
stsc.2019-10-29090750.880065.log
```
Where we also have
```console
foo@bar:~$ ls st-hippo*/ 
st-hippo1/:
W.2019-10-29090750.880065.tsv

st-hippo2/:
W.2019-10-29090750.880065.tsv

```

The content of these ''W-files'' are the results that we are interested in, every section has its own ouput folder,
whilst the proportion estimate names are shared. These files are given in a matrix format [n_spots x n_types] where each
element represent the proportion of a cell type within a specific spot. To further illustrate

```console
foo@bar:~$ head st-hippo1/W*tsv -n 5 | cut -f 1-5 | column -t
               Astrocytes_14  Astrocytes_40  Astrocytes_41  Astrocytes_42  
4.83x31.08     0.017654551    0.019365933    0.021137744    0.019530216
19.98x24.93    0.020201972    0.017161706    0.016906356    0.017648073
15.87x9.01     0.017025681    0.017077405    0.016912553    0.016671246
5.83x27.97     0.019220736    0.018679986    0.017026993    0.017045338
```
Whilst it would be legitimate to put a threshold on some of these values; for example setting all proportions lower than
a certain value to zero and renormalize. We choose to keep all values in the analysis - choose whatever strategy you
feel is most appealing.

**4.2 Visualizing results**

We include a tool for quick visualization of the results in the stereoscope package - to run this simply do:

```console
foo@bar:~$ stereoscope look -pp st-hippo*/W*tsv -o viz -sc i -sb s -nc 7 -c "umap" -g -ms 40
```
This will generate two types of images, saved to the folder ''viz''

1. _Separate visualizations_ : spots are plotted according to their array coordinates, and the intensity of their
blue facecolor corresponds to the proportion value of _each_ celltype. All of these plots are scaled internally (this highlights the
spatial patterns, but does not allow for comparison of abundance), but chancing the argument  ```-sc i``` to ```-sc s```
will scale all values within each section.

2.  _Joint visualizations_ : These are the type of ''compressed'' images described in the Method and found within the
    supplementary. Here regions of similar colors have a similar cell type composition. The method used for the
    dimensionality reduction is umap, but you can also choose between pca (slightly faster) and tsne (much slower).

Below you see some examples of the generated images:

<p align="center">
![alt text](imgs/propviz.png "single proportion visualization")
![alt text](imgs/compressed.png "joint visualization")
<p> 

You can customize the output of ```stereoscope look```in multiple ways by providing certain arguments, for more
information run 
```console
foo@bar:~$ stereoscope look -h

```

For a more sopisticated visualization you could also overlay the spots on the tissue, to better see how the spatial patterns
matches with the morphology, like we did in our figures. This is not a part of the stereoscope package, but we do
provide scripts for this.

The material that you need for such visualization is:
* HE-image - image taken of the tissue, can be scaled but not cropped
* Transformation Matrix - matrix that maps spot coordinates to pixel coordinates
* Proportion Estimates - output from stereoscope analysis
* Mask (optional)- mask to indicate which parts of the are to be included (transparent) resepectively excluded (black)

We actually use resized images (30% downscaled), since the original images are unnecessarily large for our
purpose of visualization. Still being in the res folder, all you have to do is run:

```console
foo@bar:~$ ./map2he.py -i ../data/mouse/rsc/st-hippo1.jpg -t ..data/mouse/rsc/st-hippo1-tmat.txt -p st-hippo1/W*tsv -sf 0.3 -si -o he_overlay

using image file ../data/mouse/rsc/st-hippo1.jpg
using proportions file st-hippo1/W.2019-10-29090750.880065.tsv
Rendering type Astrocytes_14
Rendering type Astrocytes_40
....
```

Resulting in images like these:

![alt text](imgs/overlay-example.png "HE-overlay")



