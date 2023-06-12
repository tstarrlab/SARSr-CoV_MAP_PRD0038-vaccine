# Count variants in each sample
This Python Jupyter notebook counts occurrences of each barcode in each sample from Illumina barcode sequencing, and adds these counts to the codon variant table.

## Set up analysis
### Import Python modules.
Use [plotnine](https://plotnine.readthedocs.io/en/stable/) for ggplot2-like plotting.

The analysis relies heavily on the Bloom lab's [dms_variants](https://jbloomlab.github.io/dms_variants) package:


```python
import glob
import itertools
import multiprocessing
import multiprocessing.pool
import os
import warnings

import alignparse
import alignparse.targets

import dms_variants.codonvarianttable
from dms_variants.constants import CBPALETTE
import dms_variants.illuminabarcodeparser
import dms_variants.utils
import dms_variants.plotnine_themes

from IPython.display import display, HTML

import pandas as pd

from plotnine import *

import yaml
```

Set [plotnine](https://plotnine.readthedocs.io/en/stable/) theme to the gray-grid one defined in `dms_variants`:


```python
theme_set(dms_variants.plotnine_themes.theme_graygrid())
```

Versions of key software:


```python
print(f"Using alignparse version {alignparse.__version__}")
print(f"Using dms_variants version {dms_variants.__version__}")
```

    Using alignparse version 0.6.1
    Using dms_variants version 1.4.3


Ignore warnings that clutter output:


```python
warnings.simplefilter('ignore')
```

### Parameters for notebook
Read the configuration file:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Make output directory if needed:


```python
os.makedirs(config['counts_dir'], exist_ok=True)
os.makedirs(config['figs_dir'], exist_ok=True)
```

## Input variant tables
Initialize the table of barcode-variant pairs from the respective bc-variant lookup tables assembled for each variant previously.


```python
variants = pd.read_csv(config['codon_variant_table'], na_filter=None)

display(HTML(variants.tail().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>target</th>
      <th>library</th>
      <th>sublibrary</th>
      <th>barcode</th>
      <th>variant_call_support</th>
      <th>codon_substitutions</th>
      <th>aa_substitutions</th>
      <th>n_codon_substitutions</th>
      <th>n_aa_substitutions</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>ZC45</td>
      <td>pool1</td>
      <td>lib47_SARSr-wts</td>
      <td>TTTTTTCTAATGGAAT</td>
      <td>1</td>
      <td></td>
      <td></td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <td>GD-Pangolin</td>
      <td>pool1</td>
      <td>lib47_SARSr-wts</td>
      <td>TTTTTTCTAGCTGGAG</td>
      <td>2</td>
      <td></td>
      <td></td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <td>AncSARS-CoV-1_alt</td>
      <td>pool1</td>
      <td>lib47_SARSr-wts</td>
      <td>TTTTTTTAGAAGACAT</td>
      <td>10</td>
      <td></td>
      <td></td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <td>RhGB01</td>
      <td>pool1</td>
      <td>lib47_SARSr-wts</td>
      <td>TTTTTTTGCGTGACAT</td>
      <td>4</td>
      <td></td>
      <td></td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <td>ZXC21</td>
      <td>pool1</td>
      <td>lib47_SARSr-wts</td>
      <td>TTTTTTTGTACATAGC</td>
      <td>2</td>
      <td></td>
      <td></td>
      <td>0</td>
      <td>0</td>
    </tr>
  </tbody>
</table>


Are there any barcodes in the same library that are shared across targets?
If so, we need to get rid of those as they will be confounded in barcode parsing:


```python
dup_barcodes = (
    variants
    .groupby(['library', 'barcode'])
    .size()
    .rename('duplicate_count')
    .reset_index()
    .query('duplicate_count > 1')
    )

print('Here are duplicated barcodes:')
display(HTML(dup_barcodes.head().to_html(index=False)))

print(f"\nRemoving the {len(dup_barcodes)} duplicated barcodes."
      f"Started with {len(variants)} barcodes:")
variants = (
    variants
    .merge(dup_barcodes, on=['library', 'barcode'], how='outer')
    .query('duplicate_count.isnull()', engine='python')
    )
print(f"After removing duplicates, there are {len(variants)} barcodes.")
```

    Here are duplicated barcodes:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>barcode</th>
      <th>duplicate_count</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>pool1</td>
      <td>AAAAGAAACATTGAAA</td>
      <td>2</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AAAAGAATCTAGAAAA</td>
      <td>2</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AACAAAACTTAAAACG</td>
      <td>2</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AACTACAATAAATGGA</td>
      <td>2</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AACTATGAACCAATAA</td>
      <td>2</td>
    </tr>
  </tbody>
</table>


    
    Removing the 29 duplicated barcodes.Started with 268757 barcodes:
    After removing duplicates, there are 268699 barcodes.


Pull out a target sequence for matching to the barcode and flanking sequence regions. Note, in this pipeline this is ok because our different backgrounds don't have differing flanks or other features within the actual N16 region covered in Illumina sequencing. If ever placing in-line barcodes here in the future, we would need to modify this.


```python
# get wildtype gene sequence for primary target
targets = alignparse.targets.Targets(seqsfile=config['amplicons_PRD-0038'],
                                     feature_parse_specs=config['feature_parse_specs_PRD-0038'])
```

## Setup to parse barcodes
Read data frame with list of all barcode runs.


```python
# barcode runs with R1 files by semicolon string split
barcode_runs = (pd.read_csv(config['barcode_runs'])
                .assign(R1=lambda x: x['R1'].str.split('; '))
                )
    
display(HTML(barcode_runs.to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>sample</th>
      <th>sample_type</th>
      <th>sort_bin</th>
      <th>concentration</th>
      <th>date</th>
      <th>number_cells</th>
      <th>R1</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>pool1</td>
      <td>mouse1-3_01_bin1</td>
      <td>mouse1-3</td>
      <td>1</td>
      <td>1</td>
      <td>230221</td>
      <td>233777</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s1_b1_S1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-3_01_bin2</td>
      <td>mouse1-3</td>
      <td>2</td>
      <td>1</td>
      <td>230221</td>
      <td>295030</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s1_b2_S2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-3_01_bin3</td>
      <td>mouse1-3</td>
      <td>3</td>
      <td>1</td>
      <td>230221</td>
      <td>1639610</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s1_b3_S3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-3_01_bin4</td>
      <td>mouse1-3</td>
      <td>4</td>
      <td>1</td>
      <td>230221</td>
      <td>4419778</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s1_b4_S4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-3_02_bin1</td>
      <td>mouse1-3</td>
      <td>1</td>
      <td>2</td>
      <td>230221</td>
      <td>486969</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s2_b1_S5_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-3_02_bin2</td>
      <td>mouse1-3</td>
      <td>2</td>
      <td>2</td>
      <td>230221</td>
      <td>898382</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s2_b2_S6_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-3_02_bin3</td>
      <td>mouse1-3</td>
      <td>3</td>
      <td>2</td>
      <td>230221</td>
      <td>3537552</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s2_b3_S7_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-3_02_bin4</td>
      <td>mouse1-3</td>
      <td>4</td>
      <td>2</td>
      <td>230221</td>
      <td>1929006</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s2_b4_S8_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-3_03_bin1</td>
      <td>mouse1-3</td>
      <td>1</td>
      <td>3</td>
      <td>230221</td>
      <td>984160</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s3_b1_S9_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-3_03_bin2</td>
      <td>mouse1-3</td>
      <td>2</td>
      <td>3</td>
      <td>230221</td>
      <td>3096110</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s3_b2_S10_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-3_03_bin3</td>
      <td>mouse1-3</td>
      <td>3</td>
      <td>3</td>
      <td>230221</td>
      <td>1447976</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s3_b3_S11_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-3_03_bin4</td>
      <td>mouse1-3</td>
      <td>4</td>
      <td>3</td>
      <td>230221</td>
      <td>981244</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s3_b4_S12_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-3_04_bin1</td>
      <td>mouse1-3</td>
      <td>1</td>
      <td>4</td>
      <td>230221</td>
      <td>2748507</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s4_b1_S13_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-3_04_bin2</td>
      <td>mouse1-3</td>
      <td>2</td>
      <td>4</td>
      <td>230221</td>
      <td>2060844</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s4_b2_S14_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-3_04_bin3</td>
      <td>mouse1-3</td>
      <td>3</td>
      <td>4</td>
      <td>230221</td>
      <td>1369837</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s4_b3_S15_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-3_04_bin4</td>
      <td>mouse1-3</td>
      <td>4</td>
      <td>4</td>
      <td>230221</td>
      <td>1250</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s4_b4_S16_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-3_05_bin1</td>
      <td>mouse1-3</td>
      <td>1</td>
      <td>5</td>
      <td>230221</td>
      <td>5729250</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s9_b1_S33_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-3_05_bin2</td>
      <td>mouse1-3</td>
      <td>2</td>
      <td>5</td>
      <td>230221</td>
      <td>324286</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s9_b2_S34_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-3_05_bin3</td>
      <td>mouse1-3</td>
      <td>3</td>
      <td>5</td>
      <td>230221</td>
      <td>70</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s9_b3_S35_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-3_05_bin4</td>
      <td>mouse1-3</td>
      <td>4</td>
      <td>5</td>
      <td>230221</td>
      <td>8</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s9_b4_S36_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-4_01_bin1</td>
      <td>mouse1-4</td>
      <td>1</td>
      <td>1</td>
      <td>230221</td>
      <td>199399</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s5_b1_S17_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-4_01_bin2</td>
      <td>mouse1-4</td>
      <td>2</td>
      <td>1</td>
      <td>230221</td>
      <td>366573</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s5_b2_S18_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-4_01_bin3</td>
      <td>mouse1-4</td>
      <td>3</td>
      <td>1</td>
      <td>230221</td>
      <td>3244015</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s5_b3_S19_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-4_01_bin4</td>
      <td>mouse1-4</td>
      <td>4</td>
      <td>1</td>
      <td>230221</td>
      <td>2520585</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s5_b4_S20_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-4_02_bin1</td>
      <td>mouse1-4</td>
      <td>1</td>
      <td>2</td>
      <td>230221</td>
      <td>458607</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s6_b1_S21_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-4_02_bin2</td>
      <td>mouse1-4</td>
      <td>2</td>
      <td>2</td>
      <td>230221</td>
      <td>2533999</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s6_b2_S22_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-4_02_bin3</td>
      <td>mouse1-4</td>
      <td>3</td>
      <td>2</td>
      <td>230221</td>
      <td>2083993</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s6_b3_S23_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-4_02_bin4</td>
      <td>mouse1-4</td>
      <td>4</td>
      <td>2</td>
      <td>230221</td>
      <td>1213113</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s6_b4_S24_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-4_03_bin1</td>
      <td>mouse1-4</td>
      <td>1</td>
      <td>3</td>
      <td>230221</td>
      <td>2188231</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s7_b1_S25_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-4_03_bin2</td>
      <td>mouse1-4</td>
      <td>2</td>
      <td>3</td>
      <td>230221</td>
      <td>2928720</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s7_b2_S26_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-4_03_bin3</td>
      <td>mouse1-4</td>
      <td>3</td>
      <td>3</td>
      <td>230221</td>
      <td>1716947</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s7_b3_S27_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-4_03_bin4</td>
      <td>mouse1-4</td>
      <td>4</td>
      <td>3</td>
      <td>230221</td>
      <td>36698</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s7_b4_S28_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-4_04_bin1</td>
      <td>mouse1-4</td>
      <td>1</td>
      <td>4</td>
      <td>230221</td>
      <td>3979399</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s8_b1_S29_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-4_04_bin2</td>
      <td>mouse1-4</td>
      <td>2</td>
      <td>4</td>
      <td>230221</td>
      <td>1879869</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s8_b2_S30_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-4_04_bin3</td>
      <td>mouse1-4</td>
      <td>3</td>
      <td>4</td>
      <td>230221</td>
      <td>234863</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s8_b3_S31_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-4_04_bin4</td>
      <td>mouse1-4</td>
      <td>4</td>
      <td>4</td>
      <td>230221</td>
      <td>9</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s8_b4_S32_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-4_05_bin1</td>
      <td>mouse1-4</td>
      <td>1</td>
      <td>5</td>
      <td>230221</td>
      <td>5729250</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s9_b1_S33_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-4_05_bin2</td>
      <td>mouse1-4</td>
      <td>2</td>
      <td>5</td>
      <td>230221</td>
      <td>324286</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s9_b2_S34_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-4_05_bin3</td>
      <td>mouse1-4</td>
      <td>3</td>
      <td>5</td>
      <td>230221</td>
      <td>70</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s9_b3_S35_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-4_05_bin4</td>
      <td>mouse1-4</td>
      <td>4</td>
      <td>5</td>
      <td>230221</td>
      <td>8</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s9_b4_S36_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-5_01_bin1</td>
      <td>mouse1-5</td>
      <td>1</td>
      <td>1</td>
      <td>230223</td>
      <td>146516</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s1_b1_S1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-5_01_bin2</td>
      <td>mouse1-5</td>
      <td>2</td>
      <td>1</td>
      <td>230223</td>
      <td>204956</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s1_b2_S2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-5_01_bin3</td>
      <td>mouse1-5</td>
      <td>3</td>
      <td>1</td>
      <td>230223</td>
      <td>1251240</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s1_b3_S3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-5_01_bin4</td>
      <td>mouse1-5</td>
      <td>4</td>
      <td>1</td>
      <td>230223</td>
      <td>3621426</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s1_b4_S4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-5_02_bin1</td>
      <td>mouse1-5</td>
      <td>1</td>
      <td>2</td>
      <td>230223</td>
      <td>259602</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s2_b1_S5_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-5_02_bin2</td>
      <td>mouse1-5</td>
      <td>2</td>
      <td>2</td>
      <td>230223</td>
      <td>825408</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s2_b2_S6_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-5_02_bin3</td>
      <td>mouse1-5</td>
      <td>3</td>
      <td>2</td>
      <td>230223</td>
      <td>3436269</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s2_b3_S7_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-5_02_bin4</td>
      <td>mouse1-5</td>
      <td>4</td>
      <td>2</td>
      <td>230223</td>
      <td>1649574</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s2_b4_S8_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-5_03_bin1</td>
      <td>mouse1-5</td>
      <td>1</td>
      <td>3</td>
      <td>230223</td>
      <td>1143115</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s3_b1_S9_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-5_03_bin2</td>
      <td>mouse1-5</td>
      <td>2</td>
      <td>3</td>
      <td>230223</td>
      <td>2911949</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s3_b2_S10_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-5_03_bin3</td>
      <td>mouse1-5</td>
      <td>3</td>
      <td>3</td>
      <td>230223</td>
      <td>1797158</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s3_b3_S11_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-5_03_bin4</td>
      <td>mouse1-5</td>
      <td>4</td>
      <td>3</td>
      <td>230223</td>
      <td>482875</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s3_b4_S12_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-5_04_bin1</td>
      <td>mouse1-5</td>
      <td>1</td>
      <td>4</td>
      <td>230223</td>
      <td>4903285</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s4_b1_S13_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-5_04_bin2</td>
      <td>mouse1-5</td>
      <td>2</td>
      <td>4</td>
      <td>230223</td>
      <td>902921</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s4_b2_S14_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-5_04_bin3</td>
      <td>mouse1-5</td>
      <td>3</td>
      <td>4</td>
      <td>230223</td>
      <td>339685</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s4_b3_S15_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-5_04_bin4</td>
      <td>mouse1-5</td>
      <td>4</td>
      <td>4</td>
      <td>230223</td>
      <td>37</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s4_b4_S16_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-5_05_bin1</td>
      <td>mouse1-5</td>
      <td>1</td>
      <td>5</td>
      <td>230223</td>
      <td>5729250</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s9_b1_S33_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-5_05_bin2</td>
      <td>mouse1-5</td>
      <td>2</td>
      <td>5</td>
      <td>230223</td>
      <td>324286</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s9_b2_S34_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-5_05_bin3</td>
      <td>mouse1-5</td>
      <td>3</td>
      <td>5</td>
      <td>230223</td>
      <td>70</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s9_b3_S35_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-5_05_bin4</td>
      <td>mouse1-5</td>
      <td>4</td>
      <td>5</td>
      <td>230223</td>
      <td>8</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s9_b4_S36_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-6_01_bin1</td>
      <td>mouse1-6</td>
      <td>1</td>
      <td>1</td>
      <td>230223</td>
      <td>108580</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s5_b1_S17_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-6_01_bin2</td>
      <td>mouse1-6</td>
      <td>2</td>
      <td>1</td>
      <td>230223</td>
      <td>498663</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s5_b2_S18_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-6_01_bin3</td>
      <td>mouse1-6</td>
      <td>3</td>
      <td>1</td>
      <td>230223</td>
      <td>2633385</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s5_b3_S19_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-6_01_bin4</td>
      <td>mouse1-6</td>
      <td>4</td>
      <td>1</td>
      <td>230223</td>
      <td>2666600</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s5_b4_S20_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-6_02_bin1</td>
      <td>mouse1-6</td>
      <td>1</td>
      <td>2</td>
      <td>230223</td>
      <td>582196</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s6_b1_S21_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-6_02_bin2</td>
      <td>mouse1-6</td>
      <td>2</td>
      <td>2</td>
      <td>230223</td>
      <td>2207214</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s6_b2_S22_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-6_02_bin3</td>
      <td>mouse1-6</td>
      <td>3</td>
      <td>2</td>
      <td>230223</td>
      <td>1922889</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s6_b3_S23_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-6_02_bin4</td>
      <td>mouse1-6</td>
      <td>4</td>
      <td>2</td>
      <td>230223</td>
      <td>1679815</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s6_b4_S24_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-6_03_bin1</td>
      <td>mouse1-6</td>
      <td>1</td>
      <td>3</td>
      <td>230223</td>
      <td>2007327</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s7_b1_S25_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-6_03_bin2</td>
      <td>mouse1-6</td>
      <td>2</td>
      <td>3</td>
      <td>230223</td>
      <td>1991086</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s7_b2_S26_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-6_03_bin3</td>
      <td>mouse1-6</td>
      <td>3</td>
      <td>3</td>
      <td>230223</td>
      <td>1210485</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s7_b3_S27_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-6_03_bin4</td>
      <td>mouse1-6</td>
      <td>4</td>
      <td>3</td>
      <td>230223</td>
      <td>884601</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s7_b4_S28_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-6_04_bin1</td>
      <td>mouse1-6</td>
      <td>1</td>
      <td>4</td>
      <td>230223</td>
      <td>3529434</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s8_b1_S29_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-6_04_bin2</td>
      <td>mouse1-6</td>
      <td>2</td>
      <td>4</td>
      <td>230223</td>
      <td>1267063</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s8_b2_S30_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-6_04_bin3</td>
      <td>mouse1-6</td>
      <td>3</td>
      <td>4</td>
      <td>230223</td>
      <td>1489152</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s8_b3_S31_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-6_04_bin4</td>
      <td>mouse1-6</td>
      <td>4</td>
      <td>4</td>
      <td>230223</td>
      <td>344</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s8_b4_S32_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-6_05_bin1</td>
      <td>mouse1-6</td>
      <td>1</td>
      <td>5</td>
      <td>230223</td>
      <td>5729250</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s9_b1_S33_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-6_05_bin2</td>
      <td>mouse1-6</td>
      <td>2</td>
      <td>5</td>
      <td>230223</td>
      <td>324286</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s9_b2_S34_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-6_05_bin3</td>
      <td>mouse1-6</td>
      <td>3</td>
      <td>5</td>
      <td>230223</td>
      <td>70</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s9_b3_S35_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse1-6_05_bin4</td>
      <td>mouse1-6</td>
      <td>4</td>
      <td>5</td>
      <td>230223</td>
      <td>8</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s9_b4_S36_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-3_01_bin1</td>
      <td>mouse2-3</td>
      <td>1</td>
      <td>1</td>
      <td>230223</td>
      <td>96146</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s9_b1_S33_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-3_01_bin2</td>
      <td>mouse2-3</td>
      <td>2</td>
      <td>1</td>
      <td>230223</td>
      <td>151199</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s9_b2_S34_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-3_01_bin3</td>
      <td>mouse2-3</td>
      <td>3</td>
      <td>1</td>
      <td>230223</td>
      <td>924775</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s9_b3_S35_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-3_01_bin4</td>
      <td>mouse2-3</td>
      <td>4</td>
      <td>1</td>
      <td>230223</td>
      <td>4643174</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s9_b4_S36_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-3_02_bin1</td>
      <td>mouse2-3</td>
      <td>1</td>
      <td>2</td>
      <td>230223</td>
      <td>185999</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s10_b1_S37_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-3_02_bin2</td>
      <td>mouse2-3</td>
      <td>2</td>
      <td>2</td>
      <td>230223</td>
      <td>349014</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s10_b2_S38_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-3_02_bin3</td>
      <td>mouse2-3</td>
      <td>3</td>
      <td>2</td>
      <td>230223</td>
      <td>3142040</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s10_b3_S39_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-3_02_bin4</td>
      <td>mouse2-3</td>
      <td>4</td>
      <td>2</td>
      <td>230223</td>
      <td>2692306</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s10_b4_S40_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-3_03_bin1</td>
      <td>mouse2-3</td>
      <td>1</td>
      <td>3</td>
      <td>230223</td>
      <td>318109</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s11_b1_S41_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-3_03_bin2</td>
      <td>mouse2-3</td>
      <td>2</td>
      <td>3</td>
      <td>230223</td>
      <td>2013025</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s11_b2_S42_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-3_03_bin3</td>
      <td>mouse2-3</td>
      <td>3</td>
      <td>3</td>
      <td>230223</td>
      <td>3088195</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s11_b3_S43_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-3_03_bin4</td>
      <td>mouse2-3</td>
      <td>4</td>
      <td>3</td>
      <td>230223</td>
      <td>968502</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s11_b4_S44_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-3_04_bin1</td>
      <td>mouse2-3</td>
      <td>1</td>
      <td>4</td>
      <td>230223</td>
      <td>1347352</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s12_b1_S45_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-3_04_bin2</td>
      <td>mouse2-3</td>
      <td>2</td>
      <td>4</td>
      <td>230223</td>
      <td>3591320</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s12_b2_S46_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-3_04_bin3</td>
      <td>mouse2-3</td>
      <td>3</td>
      <td>4</td>
      <td>230223</td>
      <td>1287363</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s12_b3_S47_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-3_04_bin4</td>
      <td>mouse2-3</td>
      <td>4</td>
      <td>4</td>
      <td>230223</td>
      <td>3292</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s12_b4_S48_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-3_05_bin1</td>
      <td>mouse2-3</td>
      <td>1</td>
      <td>5</td>
      <td>230223</td>
      <td>5729250</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s9_b1_S33_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-3_05_bin2</td>
      <td>mouse2-3</td>
      <td>2</td>
      <td>5</td>
      <td>230223</td>
      <td>324286</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s9_b2_S34_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-3_05_bin3</td>
      <td>mouse2-3</td>
      <td>3</td>
      <td>5</td>
      <td>230223</td>
      <td>70</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s9_b3_S35_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-3_05_bin4</td>
      <td>mouse2-3</td>
      <td>4</td>
      <td>5</td>
      <td>230223</td>
      <td>8</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s9_b4_S36_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-6_01_bin1</td>
      <td>mouse2-6</td>
      <td>1</td>
      <td>1</td>
      <td>230223</td>
      <td>90014</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s13_b1_S49_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-6_01_bin2</td>
      <td>mouse2-6</td>
      <td>2</td>
      <td>1</td>
      <td>230223</td>
      <td>185048</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s13_b2_S50_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-6_01_bin3</td>
      <td>mouse2-6</td>
      <td>3</td>
      <td>1</td>
      <td>230223</td>
      <td>1107553</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s13_b3_S51_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-6_01_bin4</td>
      <td>mouse2-6</td>
      <td>4</td>
      <td>1</td>
      <td>230223</td>
      <td>4833104</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s13_b4_S52_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-6_02_bin1</td>
      <td>mouse2-6</td>
      <td>1</td>
      <td>2</td>
      <td>230223</td>
      <td>149103</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s14_b1_S53_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-6_02_bin2</td>
      <td>mouse2-6</td>
      <td>2</td>
      <td>2</td>
      <td>230223</td>
      <td>323787</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s14_b2_S54_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-6_02_bin3</td>
      <td>mouse2-6</td>
      <td>3</td>
      <td>2</td>
      <td>230223</td>
      <td>3133720</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s14_b3_S55_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-6_02_bin4</td>
      <td>mouse2-6</td>
      <td>4</td>
      <td>2</td>
      <td>230223</td>
      <td>2784028</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s14_b4_S56_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-6_03_bin1</td>
      <td>mouse2-6</td>
      <td>1</td>
      <td>3</td>
      <td>230223</td>
      <td>324929</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s15_b1_S57_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-6_03_bin2</td>
      <td>mouse2-6</td>
      <td>2</td>
      <td>3</td>
      <td>230223</td>
      <td>2169845</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s15_b2_S58_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-6_03_bin3</td>
      <td>mouse2-6</td>
      <td>3</td>
      <td>3</td>
      <td>230223</td>
      <td>3636484</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s15_b3_S59_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-6_03_bin4</td>
      <td>mouse2-6</td>
      <td>4</td>
      <td>3</td>
      <td>230223</td>
      <td>536780</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s15_b4_S60_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-6_04_bin1</td>
      <td>mouse2-6</td>
      <td>1</td>
      <td>4</td>
      <td>230223</td>
      <td>1373314</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s16_b1_S61_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-6_04_bin2</td>
      <td>mouse2-6</td>
      <td>2</td>
      <td>4</td>
      <td>230223</td>
      <td>3793805</td>
      <td>[/shared/ngs/illumina/tstarr/230530_VH00699_318_AACVCNWM5/Unaligned/Project_tstarr/230223_s16_b2_S37_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-6_04_bin3</td>
      <td>mouse2-6</td>
      <td>3</td>
      <td>4</td>
      <td>230223</td>
      <td>968064</td>
      <td>[/shared/ngs/illumina/tstarr/230530_VH00699_318_AACVCNWM5/Unaligned/Project_tstarr/230223_s16_b3_S38_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-6_04_bin4</td>
      <td>mouse2-6</td>
      <td>4</td>
      <td>4</td>
      <td>230223</td>
      <td>107</td>
      <td>[/shared/ngs/illumina/tstarr/230306_VH01189_73_AACKH7JM5/Unaligned/Project_tstarr/230223_s16_b4_S64_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-6_05_bin1</td>
      <td>mouse2-6</td>
      <td>1</td>
      <td>5</td>
      <td>230223</td>
      <td>5729250</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s9_b1_S33_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-6_05_bin2</td>
      <td>mouse2-6</td>
      <td>2</td>
      <td>5</td>
      <td>230223</td>
      <td>324286</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s9_b2_S34_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-6_05_bin3</td>
      <td>mouse2-6</td>
      <td>3</td>
      <td>5</td>
      <td>230223</td>
      <td>70</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s9_b3_S35_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>mouse2-6_05_bin4</td>
      <td>mouse2-6</td>
      <td>4</td>
      <td>5</td>
      <td>230223</td>
      <td>8</td>
      <td>[/shared/ngs/illumina/tstarr/230303_VH00699_273_AACKH3GM5/Unaligned/Project_tstarr/230221_s9_b4_S36_R1_001.fastq.gz]</td>
    </tr>
  </tbody>
</table>


Make sure library / sample combinations are unique:


```python
assert len(barcode_runs) == len(barcode_runs.groupby(['library', 'sample']))
```

Make sure the the libraries for which we have barcode runs are all in our variant table:


```python
unknown_libs = set(barcode_runs['library']) - set(variants['library'])
if unknown_libs:
    raise ValueError(f"Libraries with barcode runs not in variant table: {unknown_libs}")
```

Now we initialize an [IlluminaBarcodeParser](https://jbloomlab.github.io/dms_variants/dms_variants.illuminabarcodeparser.html#dms_variants.illuminabarcodeparser.IlluminaBarcodeParser) for each library.

First, get the length of the barcode from the alignment target after making sure the same length for all targets:


```python
bclen = len(targets.targets[0].get_feature('barcode').seq)

assert (bclen == len(target.get_feature('barcode').seq) for target in targets.targets)

print(f"Barcodes of length {bclen}")
```

    Barcodes of length 16


The other barcode parsing params come from the config file:


```python
parser_params = config['illumina_barcode_parser_params']

display(HTML(
    pd.Series(parser_params, name='value')
    .rename_axis(index='parameter')
    .reset_index()
    .to_html(index=False)
    ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>parameter</th>
      <th>value</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>upstream</td>
      <td>GGCCGC</td>
    </tr>
    <tr>
      <td>downstream</td>
      <td></td>
    </tr>
    <tr>
      <td>minq</td>
      <td>20</td>
    </tr>
    <tr>
      <td>upstream_mismatch</td>
      <td>1</td>
    </tr>
    <tr>
      <td>downstream_mismatch</td>
      <td>0</td>
    </tr>
  </tbody>
</table>


The parser needs to know the set of valid barcodes, which are stored in the variant table and are different for each library.
So we create a different parser for each library using these valid barcode sets:


```python
# create dict keyed by library, value is parser for library
parsers = {lib: dms_variants.illuminabarcodeparser.IlluminaBarcodeParser(
                    bclen=bclen,
                    valid_barcodes=variants.loc[variants['library']==lib]['barcode'],
                    **parser_params)
           for lib in set(variants['library'])}

print('Number of valid barcodes searched for by each parser:')
display(HTML(
    pd.DataFrame([(lib, len(p.valid_barcodes)) for lib, p in parsers.items()],
                 columns=['library', 'number of valid barcodes'])
    .to_html(index=False)
    ))
```

    Number of valid barcodes searched for by each parser:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>number of valid barcodes</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>pool1</td>
      <td>268699</td>
    </tr>
  </tbody>
</table>


## Parse barcodes
We now parse the barcodes.
Since this will take a while, we utilize multiple CPUs via the Python [multiprocessing](https://docs.python.org/3.6/library/multiprocessing.html) module.
First, determine how many CPUs to use.
We use the minimum of the user-specified number hardcoded below and the number actually available.
(If you are running *interactively* on the Hutch cluster, you may need to reduce the number below in order to avoid an error as there is an enforced CPU limit on the home `rhino` nodes):


```python
ncpus = min(config['max_cpus'], multiprocessing.cpu_count())
print(f"Using {ncpus} CPUs")
```

    Using 16 CPUs


Parse the barcodes in parallel via a [multiprocessing.Pool](https://docs.python.org/3.6/library/multiprocessing.html#multiprocessing.pool.Pool) using all the available CPUs to get a list of the data frames with barcode counts / fates for each sample:


```python
def process_func(parser, r1files, library, sample):
    """Convenience function to be starmapped to multiprocessing pool."""
    return parser.parse(r1files, add_cols={'library': library, 'sample': sample})

# parallel computation of list of data frames
with multiprocessing.pool.Pool(processes=ncpus) as pool:
    bclist = pool.starmap(
                process_func,
                [(parsers[run.library], run.R1, run.library, run.sample)
                  for run in barcode_runs.itertuples()],
                )
```

Now concatenate the list into data frames of barcode counts and barcode fates:


```python
counts = pd.concat([samplecounts for samplecounts, _ in bclist],
                   sort=False,
                   ignore_index=True)

print('First few lines of counts data frame:')
display(HTML(counts.head().to_html(index=False)))

fates = pd.concat([samplefates for _, samplefates in bclist],
                  sort=False,
                  ignore_index=True)

print('First few lines of fates data frame:')
display(HTML(fates.head().to_html(index=False)))
```

    First few lines of counts data frame:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>barcode</th>
      <th>count</th>
      <th>library</th>
      <th>sample</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>AGTTGGCCTCACGCCA</td>
      <td>1141</td>
      <td>pool1</td>
      <td>mouse1-3_01_bin1</td>
    </tr>
    <tr>
      <td>AGAAAGCCCTGCGGTA</td>
      <td>1079</td>
      <td>pool1</td>
      <td>mouse1-3_01_bin1</td>
    </tr>
    <tr>
      <td>AGGAACCCACAGAACC</td>
      <td>937</td>
      <td>pool1</td>
      <td>mouse1-3_01_bin1</td>
    </tr>
    <tr>
      <td>ATCTAGCTCCTCAAAA</td>
      <td>933</td>
      <td>pool1</td>
      <td>mouse1-3_01_bin1</td>
    </tr>
    <tr>
      <td>ACAAATAACGCAGATT</td>
      <td>910</td>
      <td>pool1</td>
      <td>mouse1-3_01_bin1</td>
    </tr>
  </tbody>
</table>


    First few lines of fates data frame:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>fate</th>
      <th>count</th>
      <th>library</th>
      <th>sample</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>valid barcode</td>
      <td>1296455</td>
      <td>pool1</td>
      <td>mouse1-3_01_bin1</td>
    </tr>
    <tr>
      <td>invalid barcode</td>
      <td>372695</td>
      <td>pool1</td>
      <td>mouse1-3_01_bin1</td>
    </tr>
    <tr>
      <td>low quality barcode</td>
      <td>289061</td>
      <td>pool1</td>
      <td>mouse1-3_01_bin1</td>
    </tr>
    <tr>
      <td>unparseable barcode</td>
      <td>19763</td>
      <td>pool1</td>
      <td>mouse1-3_01_bin1</td>
    </tr>
    <tr>
      <td>failed chastity filter</td>
      <td>0</td>
      <td>pool1</td>
      <td>mouse1-3_01_bin1</td>
    </tr>
  </tbody>
</table>


## Examine fates of parsed barcodes
First, we'll analyze the "fates" of the parsed barcodes.
These fates represent what happened to each Illumina read we parsed:
 - Did the barcode read fail the Illumina chastity filter?
 - Was the barcode *unparseable* (i.e., the read didn't appear to be a valid barcode based on flanking regions)?
 - Was the barcode sequence too *low quality* based on the Illumina quality scores?
 - Was the barcode parseable but *invalid* (i.e., not in our list of variant-associated barcodes in the codon variant table)?
 - Was the barcode *valid*, and so will be added to variant counts.
 
First, we just write a CSV file with all the barcode fates:


```python
fatesfile = os.path.join(config['counts_dir'], 'barcode_fates.csv')
print(f"Writing barcode fates to {fatesfile}")
fates.to_csv(fatesfile, index=False)
```

    Writing barcode fates to results/counts/barcode_fates.csv


Next, we tabulate the barcode fates in wide format:


```python
display(HTML(fates
             .pivot_table(columns='fate',
                          values='count',
                          index=['library', 'sample'])
             .to_html()
             ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>fate</th>
      <th>failed chastity filter</th>
      <th>invalid barcode</th>
      <th>low quality barcode</th>
      <th>unparseable barcode</th>
      <th>valid barcode</th>
    </tr>
    <tr>
      <th>library</th>
      <th>sample</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="120" valign="top">pool1</th>
      <th>mouse1-3_01_bin1</th>
      <td>0</td>
      <td>372695</td>
      <td>289061</td>
      <td>19763</td>
      <td>1296455</td>
    </tr>
    <tr>
      <th>mouse1-3_01_bin2</th>
      <td>0</td>
      <td>578923</td>
      <td>485652</td>
      <td>35411</td>
      <td>2228211</td>
    </tr>
    <tr>
      <th>mouse1-3_01_bin3</th>
      <td>0</td>
      <td>1522587</td>
      <td>1314449</td>
      <td>88769</td>
      <td>5840982</td>
    </tr>
    <tr>
      <th>mouse1-3_01_bin4</th>
      <td>0</td>
      <td>6151156</td>
      <td>4996841</td>
      <td>332667</td>
      <td>21007603</td>
    </tr>
    <tr>
      <th>mouse1-3_02_bin1</th>
      <td>0</td>
      <td>500124</td>
      <td>388521</td>
      <td>27035</td>
      <td>1740436</td>
    </tr>
    <tr>
      <th>mouse1-3_02_bin2</th>
      <td>0</td>
      <td>870548</td>
      <td>706705</td>
      <td>50451</td>
      <td>3181864</td>
    </tr>
    <tr>
      <th>mouse1-3_02_bin3</th>
      <td>0</td>
      <td>4599417</td>
      <td>3462274</td>
      <td>235348</td>
      <td>15394693</td>
    </tr>
    <tr>
      <th>mouse1-3_02_bin4</th>
      <td>0</td>
      <td>1918653</td>
      <td>1610560</td>
      <td>119071</td>
      <td>7459371</td>
    </tr>
    <tr>
      <th>mouse1-3_03_bin1</th>
      <td>0</td>
      <td>1440266</td>
      <td>1154175</td>
      <td>83707</td>
      <td>5190333</td>
    </tr>
    <tr>
      <th>mouse1-3_03_bin2</th>
      <td>0</td>
      <td>4312218</td>
      <td>3349653</td>
      <td>215732</td>
      <td>14432986</td>
    </tr>
    <tr>
      <th>mouse1-3_03_bin3</th>
      <td>0</td>
      <td>1251521</td>
      <td>1128165</td>
      <td>79232</td>
      <td>5150222</td>
    </tr>
    <tr>
      <th>mouse1-3_03_bin4</th>
      <td>0</td>
      <td>998062</td>
      <td>765424</td>
      <td>57165</td>
      <td>3466180</td>
    </tr>
    <tr>
      <th>mouse1-3_04_bin1</th>
      <td>0</td>
      <td>3029474</td>
      <td>2503205</td>
      <td>165267</td>
      <td>10539544</td>
    </tr>
    <tr>
      <th>mouse1-3_04_bin2</th>
      <td>0</td>
      <td>3331679</td>
      <td>2740204</td>
      <td>184205</td>
      <td>12104693</td>
    </tr>
    <tr>
      <th>mouse1-3_04_bin3</th>
      <td>0</td>
      <td>1944530</td>
      <td>1535238</td>
      <td>111485</td>
      <td>6683197</td>
    </tr>
    <tr>
      <th>mouse1-3_04_bin4</th>
      <td>0</td>
      <td>2481</td>
      <td>2750</td>
      <td>132</td>
      <td>8999</td>
    </tr>
    <tr>
      <th>mouse1-3_05_bin1</th>
      <td>0</td>
      <td>5277023</td>
      <td>4265057</td>
      <td>291300</td>
      <td>18447823</td>
    </tr>
    <tr>
      <th>mouse1-3_05_bin2</th>
      <td>0</td>
      <td>390900</td>
      <td>321063</td>
      <td>21809</td>
      <td>1416811</td>
    </tr>
    <tr>
      <th>mouse1-3_05_bin3</th>
      <td>0</td>
      <td>57</td>
      <td>1270</td>
      <td>7</td>
      <td>296</td>
    </tr>
    <tr>
      <th>mouse1-3_05_bin4</th>
      <td>0</td>
      <td>40</td>
      <td>605</td>
      <td>4</td>
      <td>130</td>
    </tr>
    <tr>
      <th>mouse1-4_01_bin1</th>
      <td>0</td>
      <td>281544</td>
      <td>229436</td>
      <td>16293</td>
      <td>1004744</td>
    </tr>
    <tr>
      <th>mouse1-4_01_bin2</th>
      <td>0</td>
      <td>551464</td>
      <td>459730</td>
      <td>33728</td>
      <td>2004621</td>
    </tr>
    <tr>
      <th>mouse1-4_01_bin3</th>
      <td>0</td>
      <td>5271222</td>
      <td>4224178</td>
      <td>282462</td>
      <td>17795407</td>
    </tr>
    <tr>
      <th>mouse1-4_01_bin4</th>
      <td>0</td>
      <td>3234402</td>
      <td>2680693</td>
      <td>181344</td>
      <td>11824261</td>
    </tr>
    <tr>
      <th>mouse1-4_02_bin1</th>
      <td>0</td>
      <td>506038</td>
      <td>424731</td>
      <td>29503</td>
      <td>1825975</td>
    </tr>
    <tr>
      <th>mouse1-4_02_bin2</th>
      <td>0</td>
      <td>3578532</td>
      <td>2878848</td>
      <td>196183</td>
      <td>12398472</td>
    </tr>
    <tr>
      <th>mouse1-4_02_bin3</th>
      <td>0</td>
      <td>2842507</td>
      <td>2341210</td>
      <td>156822</td>
      <td>10229447</td>
    </tr>
    <tr>
      <th>mouse1-4_02_bin4</th>
      <td>0</td>
      <td>1357916</td>
      <td>1085971</td>
      <td>77360</td>
      <td>4733565</td>
    </tr>
    <tr>
      <th>mouse1-4_03_bin1</th>
      <td>0</td>
      <td>1884122</td>
      <td>1541324</td>
      <td>104813</td>
      <td>6611322</td>
    </tr>
    <tr>
      <th>mouse1-4_03_bin2</th>
      <td>0</td>
      <td>2486740</td>
      <td>2020269</td>
      <td>133844</td>
      <td>8820753</td>
    </tr>
    <tr>
      <th>mouse1-4_03_bin3</th>
      <td>0</td>
      <td>1535677</td>
      <td>1246046</td>
      <td>87018</td>
      <td>5367016</td>
    </tr>
    <tr>
      <th>mouse1-4_03_bin4</th>
      <td>0</td>
      <td>149974</td>
      <td>119233</td>
      <td>8910</td>
      <td>526215</td>
    </tr>
    <tr>
      <th>mouse1-4_04_bin1</th>
      <td>0</td>
      <td>5257745</td>
      <td>4189379</td>
      <td>284205</td>
      <td>18410683</td>
    </tr>
    <tr>
      <th>mouse1-4_04_bin2</th>
      <td>0</td>
      <td>1617748</td>
      <td>1324735</td>
      <td>93260</td>
      <td>5759953</td>
    </tr>
    <tr>
      <th>mouse1-4_04_bin3</th>
      <td>0</td>
      <td>315018</td>
      <td>249301</td>
      <td>18372</td>
      <td>1104349</td>
    </tr>
    <tr>
      <th>mouse1-4_04_bin4</th>
      <td>0</td>
      <td>6115</td>
      <td>1338</td>
      <td>85</td>
      <td>132</td>
    </tr>
    <tr>
      <th>mouse1-4_05_bin1</th>
      <td>0</td>
      <td>5277023</td>
      <td>4265057</td>
      <td>291300</td>
      <td>18447823</td>
    </tr>
    <tr>
      <th>mouse1-4_05_bin2</th>
      <td>0</td>
      <td>390900</td>
      <td>321063</td>
      <td>21809</td>
      <td>1416811</td>
    </tr>
    <tr>
      <th>mouse1-4_05_bin3</th>
      <td>0</td>
      <td>57</td>
      <td>1270</td>
      <td>7</td>
      <td>296</td>
    </tr>
    <tr>
      <th>mouse1-4_05_bin4</th>
      <td>0</td>
      <td>40</td>
      <td>605</td>
      <td>4</td>
      <td>130</td>
    </tr>
    <tr>
      <th>mouse1-5_01_bin1</th>
      <td>0</td>
      <td>126079</td>
      <td>98398</td>
      <td>6344</td>
      <td>412836</td>
    </tr>
    <tr>
      <th>mouse1-5_01_bin2</th>
      <td>0</td>
      <td>211115</td>
      <td>189199</td>
      <td>12698</td>
      <td>803636</td>
    </tr>
    <tr>
      <th>mouse1-5_01_bin3</th>
      <td>0</td>
      <td>1189743</td>
      <td>1114178</td>
      <td>72363</td>
      <td>4627787</td>
    </tr>
    <tr>
      <th>mouse1-5_01_bin4</th>
      <td>0</td>
      <td>2980366</td>
      <td>2669488</td>
      <td>167080</td>
      <td>10570228</td>
    </tr>
    <tr>
      <th>mouse1-5_02_bin1</th>
      <td>0</td>
      <td>261285</td>
      <td>208859</td>
      <td>14208</td>
      <td>868572</td>
    </tr>
    <tr>
      <th>mouse1-5_02_bin2</th>
      <td>0</td>
      <td>628075</td>
      <td>521529</td>
      <td>35749</td>
      <td>2199078</td>
    </tr>
    <tr>
      <th>mouse1-5_02_bin3</th>
      <td>0</td>
      <td>3009569</td>
      <td>2453942</td>
      <td>159329</td>
      <td>10450048</td>
    </tr>
    <tr>
      <th>mouse1-5_02_bin4</th>
      <td>0</td>
      <td>1546260</td>
      <td>1370565</td>
      <td>98687</td>
      <td>5973053</td>
    </tr>
    <tr>
      <th>mouse1-5_03_bin1</th>
      <td>0</td>
      <td>1252277</td>
      <td>1032241</td>
      <td>68653</td>
      <td>4316106</td>
    </tr>
    <tr>
      <th>mouse1-5_03_bin2</th>
      <td>0</td>
      <td>2809578</td>
      <td>2340582</td>
      <td>148100</td>
      <td>9506982</td>
    </tr>
    <tr>
      <th>mouse1-5_03_bin3</th>
      <td>0</td>
      <td>1564041</td>
      <td>1477855</td>
      <td>96698</td>
      <td>6224148</td>
    </tr>
    <tr>
      <th>mouse1-5_03_bin4</th>
      <td>0</td>
      <td>422984</td>
      <td>352638</td>
      <td>25363</td>
      <td>1504642</td>
    </tr>
    <tr>
      <th>mouse1-5_04_bin1</th>
      <td>0</td>
      <td>3838323</td>
      <td>3285027</td>
      <td>212690</td>
      <td>13152331</td>
    </tr>
    <tr>
      <th>mouse1-5_04_bin2</th>
      <td>0</td>
      <td>954920</td>
      <td>830802</td>
      <td>53507</td>
      <td>3463833</td>
    </tr>
    <tr>
      <th>mouse1-5_04_bin3</th>
      <td>0</td>
      <td>297659</td>
      <td>256563</td>
      <td>16741</td>
      <td>1050642</td>
    </tr>
    <tr>
      <th>mouse1-5_04_bin4</th>
      <td>0</td>
      <td>38</td>
      <td>373</td>
      <td>9</td>
      <td>95</td>
    </tr>
    <tr>
      <th>mouse1-5_05_bin1</th>
      <td>0</td>
      <td>5277023</td>
      <td>4265057</td>
      <td>291300</td>
      <td>18447823</td>
    </tr>
    <tr>
      <th>mouse1-5_05_bin2</th>
      <td>0</td>
      <td>390900</td>
      <td>321063</td>
      <td>21809</td>
      <td>1416811</td>
    </tr>
    <tr>
      <th>mouse1-5_05_bin3</th>
      <td>0</td>
      <td>57</td>
      <td>1270</td>
      <td>7</td>
      <td>296</td>
    </tr>
    <tr>
      <th>mouse1-5_05_bin4</th>
      <td>0</td>
      <td>40</td>
      <td>605</td>
      <td>4</td>
      <td>130</td>
    </tr>
    <tr>
      <th>mouse1-6_01_bin1</th>
      <td>0</td>
      <td>121454</td>
      <td>97885</td>
      <td>6517</td>
      <td>393164</td>
    </tr>
    <tr>
      <th>mouse1-6_01_bin2</th>
      <td>0</td>
      <td>511081</td>
      <td>460398</td>
      <td>29697</td>
      <td>1901489</td>
    </tr>
    <tr>
      <th>mouse1-6_01_bin3</th>
      <td>0</td>
      <td>2838154</td>
      <td>2517411</td>
      <td>160833</td>
      <td>10064089</td>
    </tr>
    <tr>
      <th>mouse1-6_01_bin4</th>
      <td>0</td>
      <td>2908027</td>
      <td>2530237</td>
      <td>164060</td>
      <td>10363743</td>
    </tr>
    <tr>
      <th>mouse1-6_02_bin1</th>
      <td>0</td>
      <td>535274</td>
      <td>454837</td>
      <td>29471</td>
      <td>1810543</td>
    </tr>
    <tr>
      <th>mouse1-6_02_bin2</th>
      <td>0</td>
      <td>2281443</td>
      <td>1930222</td>
      <td>122342</td>
      <td>7767851</td>
    </tr>
    <tr>
      <th>mouse1-6_02_bin3</th>
      <td>0</td>
      <td>2170464</td>
      <td>1900129</td>
      <td>121135</td>
      <td>7816672</td>
    </tr>
    <tr>
      <th>mouse1-6_02_bin4</th>
      <td>0</td>
      <td>1710885</td>
      <td>1566362</td>
      <td>108319</td>
      <td>6475162</td>
    </tr>
    <tr>
      <th>mouse1-6_03_bin1</th>
      <td>0</td>
      <td>1793752</td>
      <td>1515268</td>
      <td>98080</td>
      <td>6045852</td>
    </tr>
    <tr>
      <th>mouse1-6_03_bin2</th>
      <td>0</td>
      <td>1789759</td>
      <td>1507580</td>
      <td>93981</td>
      <td>6147138</td>
    </tr>
    <tr>
      <th>mouse1-6_03_bin3</th>
      <td>0</td>
      <td>770574</td>
      <td>833268</td>
      <td>53073</td>
      <td>3531491</td>
    </tr>
    <tr>
      <th>mouse1-6_03_bin4</th>
      <td>0</td>
      <td>763597</td>
      <td>659677</td>
      <td>46805</td>
      <td>2732094</td>
    </tr>
    <tr>
      <th>mouse1-6_04_bin1</th>
      <td>0</td>
      <td>2789291</td>
      <td>2271529</td>
      <td>146053</td>
      <td>9271845</td>
    </tr>
    <tr>
      <th>mouse1-6_04_bin2</th>
      <td>0</td>
      <td>811787</td>
      <td>829877</td>
      <td>51970</td>
      <td>3486115</td>
    </tr>
    <tr>
      <th>mouse1-6_04_bin3</th>
      <td>0</td>
      <td>1230693</td>
      <td>1060877</td>
      <td>70481</td>
      <td>4376027</td>
    </tr>
    <tr>
      <th>mouse1-6_04_bin4</th>
      <td>0</td>
      <td>537</td>
      <td>818</td>
      <td>22</td>
      <td>1178</td>
    </tr>
    <tr>
      <th>mouse1-6_05_bin1</th>
      <td>0</td>
      <td>5277023</td>
      <td>4265057</td>
      <td>291300</td>
      <td>18447823</td>
    </tr>
    <tr>
      <th>mouse1-6_05_bin2</th>
      <td>0</td>
      <td>390900</td>
      <td>321063</td>
      <td>21809</td>
      <td>1416811</td>
    </tr>
    <tr>
      <th>mouse1-6_05_bin3</th>
      <td>0</td>
      <td>57</td>
      <td>1270</td>
      <td>7</td>
      <td>296</td>
    </tr>
    <tr>
      <th>mouse1-6_05_bin4</th>
      <td>0</td>
      <td>40</td>
      <td>605</td>
      <td>4</td>
      <td>130</td>
    </tr>
    <tr>
      <th>mouse2-3_01_bin1</th>
      <td>0</td>
      <td>82172</td>
      <td>76095</td>
      <td>5203</td>
      <td>309124</td>
    </tr>
    <tr>
      <th>mouse2-3_01_bin2</th>
      <td>0</td>
      <td>129531</td>
      <td>124374</td>
      <td>8064</td>
      <td>534971</td>
    </tr>
    <tr>
      <th>mouse2-3_01_bin3</th>
      <td>0</td>
      <td>663314</td>
      <td>673803</td>
      <td>44594</td>
      <td>2910087</td>
    </tr>
    <tr>
      <th>mouse2-3_01_bin4</th>
      <td>0</td>
      <td>4523335</td>
      <td>3893733</td>
      <td>253317</td>
      <td>15836248</td>
    </tr>
    <tr>
      <th>mouse2-3_02_bin1</th>
      <td>0</td>
      <td>145675</td>
      <td>135570</td>
      <td>9197</td>
      <td>568076</td>
    </tr>
    <tr>
      <th>mouse2-3_02_bin2</th>
      <td>0</td>
      <td>236652</td>
      <td>236221</td>
      <td>15827</td>
      <td>1013704</td>
    </tr>
    <tr>
      <th>mouse2-3_02_bin3</th>
      <td>0</td>
      <td>2053993</td>
      <td>1972221</td>
      <td>128279</td>
      <td>8340427</td>
    </tr>
    <tr>
      <th>mouse2-3_02_bin4</th>
      <td>0</td>
      <td>2339535</td>
      <td>1909730</td>
      <td>119631</td>
      <td>7744417</td>
    </tr>
    <tr>
      <th>mouse2-3_03_bin1</th>
      <td>0</td>
      <td>201366</td>
      <td>195644</td>
      <td>13284</td>
      <td>811430</td>
    </tr>
    <tr>
      <th>mouse2-3_03_bin2</th>
      <td>0</td>
      <td>1197456</td>
      <td>1140562</td>
      <td>75182</td>
      <td>4893497</td>
    </tr>
    <tr>
      <th>mouse2-3_03_bin3</th>
      <td>0</td>
      <td>2113259</td>
      <td>1910802</td>
      <td>124349</td>
      <td>8046852</td>
    </tr>
    <tr>
      <th>mouse2-3_03_bin4</th>
      <td>0</td>
      <td>1074173</td>
      <td>737076</td>
      <td>46978</td>
      <td>2757209</td>
    </tr>
    <tr>
      <th>mouse2-3_04_bin1</th>
      <td>0</td>
      <td>720263</td>
      <td>679086</td>
      <td>45119</td>
      <td>2907524</td>
    </tr>
    <tr>
      <th>mouse2-3_04_bin2</th>
      <td>0</td>
      <td>2297027</td>
      <td>2135493</td>
      <td>144214</td>
      <td>8950129</td>
    </tr>
    <tr>
      <th>mouse2-3_04_bin3</th>
      <td>0</td>
      <td>1606116</td>
      <td>1093149</td>
      <td>68446</td>
      <td>4220438</td>
    </tr>
    <tr>
      <th>mouse2-3_04_bin4</th>
      <td>0</td>
      <td>4434</td>
      <td>3440</td>
      <td>241</td>
      <td>12748</td>
    </tr>
    <tr>
      <th>mouse2-3_05_bin1</th>
      <td>0</td>
      <td>5277023</td>
      <td>4265057</td>
      <td>291300</td>
      <td>18447823</td>
    </tr>
    <tr>
      <th>mouse2-3_05_bin2</th>
      <td>0</td>
      <td>390900</td>
      <td>321063</td>
      <td>21809</td>
      <td>1416811</td>
    </tr>
    <tr>
      <th>mouse2-3_05_bin3</th>
      <td>0</td>
      <td>57</td>
      <td>1270</td>
      <td>7</td>
      <td>296</td>
    </tr>
    <tr>
      <th>mouse2-3_05_bin4</th>
      <td>0</td>
      <td>40</td>
      <td>605</td>
      <td>4</td>
      <td>130</td>
    </tr>
    <tr>
      <th>mouse2-6_01_bin1</th>
      <td>0</td>
      <td>77568</td>
      <td>70119</td>
      <td>5167</td>
      <td>291926</td>
    </tr>
    <tr>
      <th>mouse2-6_01_bin2</th>
      <td>0</td>
      <td>130734</td>
      <td>128821</td>
      <td>8537</td>
      <td>545974</td>
    </tr>
    <tr>
      <th>mouse2-6_01_bin3</th>
      <td>0</td>
      <td>614256</td>
      <td>623193</td>
      <td>41954</td>
      <td>2685275</td>
    </tr>
    <tr>
      <th>mouse2-6_01_bin4</th>
      <td>0</td>
      <td>4763536</td>
      <td>3074569</td>
      <td>196213</td>
      <td>11463190</td>
    </tr>
    <tr>
      <th>mouse2-6_02_bin1</th>
      <td>0</td>
      <td>105026</td>
      <td>97194</td>
      <td>6682</td>
      <td>403449</td>
    </tr>
    <tr>
      <th>mouse2-6_02_bin2</th>
      <td>0</td>
      <td>203489</td>
      <td>208281</td>
      <td>14277</td>
      <td>884915</td>
    </tr>
    <tr>
      <th>mouse2-6_02_bin3</th>
      <td>0</td>
      <td>1768228</td>
      <td>1727081</td>
      <td>113027</td>
      <td>7145410</td>
    </tr>
    <tr>
      <th>mouse2-6_02_bin4</th>
      <td>0</td>
      <td>3575460</td>
      <td>2976947</td>
      <td>188759</td>
      <td>11751023</td>
    </tr>
    <tr>
      <th>mouse2-6_03_bin1</th>
      <td>0</td>
      <td>161204</td>
      <td>153638</td>
      <td>10290</td>
      <td>654187</td>
    </tr>
    <tr>
      <th>mouse2-6_03_bin2</th>
      <td>0</td>
      <td>1528075</td>
      <td>1495031</td>
      <td>99516</td>
      <td>6332850</td>
    </tr>
    <tr>
      <th>mouse2-6_03_bin3</th>
      <td>0</td>
      <td>4169800</td>
      <td>3654540</td>
      <td>238093</td>
      <td>15071135</td>
    </tr>
    <tr>
      <th>mouse2-6_03_bin4</th>
      <td>0</td>
      <td>764619</td>
      <td>507002</td>
      <td>31678</td>
      <td>1942178</td>
    </tr>
    <tr>
      <th>mouse2-6_04_bin1</th>
      <td>0</td>
      <td>1214369</td>
      <td>1173524</td>
      <td>78322</td>
      <td>4900839</td>
    </tr>
    <tr>
      <th>mouse2-6_04_bin2</th>
      <td>0</td>
      <td>1880905</td>
      <td>1619782</td>
      <td>111100</td>
      <td>7068258</td>
    </tr>
    <tr>
      <th>mouse2-6_04_bin3</th>
      <td>0</td>
      <td>680853</td>
      <td>428773</td>
      <td>28580</td>
      <td>1714308</td>
    </tr>
    <tr>
      <th>mouse2-6_04_bin4</th>
      <td>0</td>
      <td>24</td>
      <td>30</td>
      <td>0</td>
      <td>86</td>
    </tr>
    <tr>
      <th>mouse2-6_05_bin1</th>
      <td>0</td>
      <td>5277023</td>
      <td>4265057</td>
      <td>291300</td>
      <td>18447823</td>
    </tr>
    <tr>
      <th>mouse2-6_05_bin2</th>
      <td>0</td>
      <td>390900</td>
      <td>321063</td>
      <td>21809</td>
      <td>1416811</td>
    </tr>
    <tr>
      <th>mouse2-6_05_bin3</th>
      <td>0</td>
      <td>57</td>
      <td>1270</td>
      <td>7</td>
      <td>296</td>
    </tr>
    <tr>
      <th>mouse2-6_05_bin4</th>
      <td>0</td>
      <td>40</td>
      <td>605</td>
      <td>4</td>
      <td>130</td>
    </tr>
  </tbody>
</table>


Now we plot the barcode-read fates for each library / sample, showing the bars for valid barcodes in orange and the others in gray.
We see that the largest fraction of barcode reads correspond to valid barcodes, and most of the others are invalid barcodes (probably because the map to variants that aren't present in our variant table since we didn't associate all variants with barcodes). The exception to this is lib2 Titeseq_03_bin3; the PCR for this sample in the original sequencing run failed, so we followed it up with a single MiSeq lane. We did not filter out the PhiX reads from this data before parsing, so these PhiX reads will deflate the fraction of valid barcode reads as expected, but does not indicate any problems.


```python
barcode_fate_plot = (
    ggplot(
        fates
        .assign(sample=lambda x: pd.Categorical(x['sample'],
                                                x['sample'].unique(),
                                                ordered=True),
                fate=lambda x: pd.Categorical(x['fate'],
                                              x['fate'].unique(),
                                              ordered=True),
                is_valid=lambda x: x['fate'] == 'valid barcode'
                ), 
        aes('fate', 'count', fill='is_valid')) +
    geom_bar(stat='identity') +
    facet_grid('sample ~ library') +
    facet_grid('sample ~ library') +
    scale_fill_manual(CBPALETTE, guide=False) +
    theme(figure_size=(1.4 * (1 + fates['library'].nunique()),
                       1.7 * (1.2 + fates['sample'].nunique())),
          axis_text_x=element_text(angle=90),
          panel_grid_major_x=element_blank()
          ) +
    scale_y_continuous(labels=dms_variants.utils.latex_sci_not,
                       name='number of reads')
    )

_ = barcode_fate_plot.draw()
```


    
![png](count_variants_files/count_variants_42_0.png)
    


## Output csv of barcode counts in variant-barcode lookup table


```python
print(f"Writing variant counts to {config['variant_counts_file']}")
counts.to_csv(config['variant_counts_file'], index=False)
```

    Writing variant counts to results/counts/variant_counts.csv


The [CodonVariantTable](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html#dms_variants.codonvarianttable.CodonVariantTable) has lots of nice functions that can be used to analyze the counts it contains.
However, we do that in the next notebook so we don't have to re-run this entire (rather computationally intensive) notebook every time we want to analyze a new aspect of the counts.


```python

```
