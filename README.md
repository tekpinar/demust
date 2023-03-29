# demust

*A Python toolkit to prepare, modify, visualize and analyze deep mutational scanning (DMS) data of proteins.*

demust contains several scripts: plots, maps, convert, compare and removegaps, extract, compare, inputgenerator etc. 

## Installation
Download the repository:
```bash
git clone https://github.com/tekpinar/demust.git
```


Go to the demust folder. Then, run the following command:

```bash
pip install -e .
```


## Quick Start

You can invoke help about each demust script as follows:

demust maps -h
demust plots -h
demust convert -h
demust compare -h
demust removegaps -h

demust maps tool can help visualize GEMME, FoldX, Rhapsody and various other Deep Mutational Scanning (DMS)
data.

demust plots tool extracts data (such as averages etc) from the text data output of GEMME or FoldX.
It is supposed to be a 2D plots tool.

demust convert will convert text output from (to) GEMME to (from) FoldX and Rhapsody formats. 

demust compare will give you Spearman correlation between two DMS data sets.

demust inputgenerator will produce input files for polyphen2 from fasta file. 

## Cite
We recommend you to put a link to the github repository for now. 

## Licensing

*demust* is developed and released under [GNU Lesser GPL Licence](https://www.gnu.org/licenses/lgpl-3.0.en.html). 
Please read to the **COPYING** and **COPYING.LESSER** files to know more. 
