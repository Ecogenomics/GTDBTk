
# Prodigal

## Translation table selection

### About
Due to alternative codon usage, the genes called may vary.
That variation can vary greatly depending on the translation table selected and have a significant impact on
 subsequent processes. It is therefore essential to select the correct translation table.

GTDB-Tk concerns itself with the following translation tables:

* [[4] The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG4)
* [[11] The Bacterial, Archaeal and Plant Plastid Code](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG11)
 
All gene calling is done through the external tool [Prodigal](https://github.com/hyattpd/Prodigal). As of writing this,
Prodigal does not yet support automatic translation table selection. Therefore, we have developed a binary classification
model to predict the correct translation table. That model was tested on NCBI release 86.


### Why was translation table 25 ignored?
In short, it's indistinguishable from translation table 4 (using Prodigal gene scores).

The gene scores output by Prodigal were identical for both translation tables 4 and 25. 
Translation tables 4 and 25 are very similar in that UGA no longer encodes a STOP codon, and instead encodes
Trp (table 4) and Gly (table 11). As a result, the called genes must be directly compared, bringing its own
unique challenges and computational burden.
 
The decision was made to treat translation table 25 as equivalent to translation table 4, due to the under-representation of 
organisms using this translation table (0.013% of NCBI release 86) and the single codon difference.

In a future release of Prodigal, a score for the stop codon will be available (`mscore`) which will hopefully lead 
to an obvious distinction between tables 4 and 25.

### Prodigal gene statistics

Please refer to the [Prodigal documentation](https://github.com/hyattpd/prodigal/wiki/understanding-the-prodigal-output#gene-coordinates) 
for an explanation of the features used. 

Given that only translation tables 4 and 11 were considered, the model was built as a comparison of these distributions.
For each translation table, the median, standard deviation, and maximum were 
calculated for each of the following features: `conf`, `cscore`, `gc_cont`, `rscore`, `tscore`, and `uscore`. 
Additionally, the number of genes called and the coding density were included.


### Model performance
The coding density is a very good indicator of which translation table is true, typically, a higher coding density
indicates the correct translation table.  Selecting the translation table which produced the highest coding density 
resulted in a classification accuracy of: 99.9591% (n=117482).

In addition to using the coding density as a predictor, the aforementioned statistics were included in a logistic
classifier. Using the logistic classifier model to select the highest probability translation table resulted in an accuracy of: 
99.9881% (n=117482, cv=10).

### Caveats
If the number of nucleotides are less than 100,000 Prodigal is run using the `meta` parameter which
produces an identical result for translation tables 4, and 11. In this case, genes will be called using
translation table 11.