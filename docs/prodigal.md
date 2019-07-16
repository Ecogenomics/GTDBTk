
# Prodigal

## Translation table selection

### About
Due to alternative codon usage, the genes called may vary. That variation can 
vary greatly depending on the translation table selected and have a significant 
impact on subsequent processes. It is therefore essential to select the correct 
translation table.

GTDB-Tk concerns itself with the following translation tables:

* [[4] The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG4)
* [[11] The Bacterial, Archaeal and Plant Plastid Code](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG11)
 
All gene calling is done through the external tool [Prodigal](https://github.com/hyattpd/Prodigal). 
As of writing this, Prodigal does not yet support automatic translation table 
selection. Therefore, we have developed a binary classification model to predict 
the correct translation table. That model was trained on NCBI release 86.


### Why was translation table 25 ignored?
In short, it's indistinguishable from translation table 4 (using Prodigal gene 
scores).

The gene scores output by Prodigal were identical for both translation tables 4 
and 25. Translation tables 4 and 25 are very similar and only differ from
translation table 11 by the UGA codon no longer encoding STOP, and instead 
coding Trp (table 4) and Gly (table 11). In order to compare translation tables
4 and 25, the called genes must be directly compared, bringing its own unique 
challenges and computational burden.
 
The decision was made to treat translation table 25 as equivalent to translation
table 4 out of necessity, due to the under-representation of organisms using 
translation table 25 (0.013% of NCBI release 86) and the single codon difference.

In a future release of Prodigal, a score for the stop codon will be available 
(`mscore`) which will hopefully lead to an obvious distinction between tables 4 
and 25.

### Prodigal gene statistics

Please refer to the [Prodigal documentation](https://github.com/hyattpd/prodigal/wiki/understanding-the-prodigal-output#gene-coordinates) 
for an explanation of the features used. 

Given that only translation tables 4 and 11 were considered, the model was built
as a comparison of these distributions. For each translation table, the median, 
standard deviation, and maximum were calculated for each of the following 
features: `conf`, `cscore`, `gc_cont`, `rscore`, `tscore`, and `uscore`. 
Additionally, the number of genes called and the coding density were included.


### Model performance
##### Validation (NCBI release 86)
The coding density is a very good indicator of which translation table is true, 
typically, a higher coding density indicates the correct translation table. 
The naive approach of selecting the translation table on the table which 
produced the highest coding density resulted in a classification accuracy of: 
99.9591% (n=117,482).

In addition to using the coding density as a predictor, the aforementioned 
statistics were included in a logistic classifier. Using the logistic classifier
model to select the highest probability translation table resulted in an 
accuracy of: 99.9881% (n=117,482, cv=10).

##### Testing (NCBI release 89)
The model was tested on genomes from NCBI release 89 which were not present or
updated from release 86, and which also passed quality control. The model
predicted the correct translation table with an accuracy of 99.9814% (n=16,101).
Three genomes were incorrectly predicted  as translation table 4 but were truly 
translation table 11 (GCA_003029925.1, GCA_900118515.1, GCA_900118825.1).

In comparison, the naive approach of selecting the highest coding density on the
same set of genomes results in an accuracy of 71.4055%, where 4604 translation 
table 11 genomes were misclassified as translation table 4.


### Caveats
If the number of nucleotides are less than 100,000 Prodigal is run using the 
`meta` parameter which produces an identical result for translation tables 4, 
and 11. In this case, genes will be called using translation table 11.

### Future work
A rolling model will be developed where the classifier is trained and validated
on all previous NCBI releases, with the exception unseen genomes in the most 
recent release which will be used as a test data set.
