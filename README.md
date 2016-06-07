# moloco
A tool for detecting shared genetic variants across multiple phenotypes using GWAS summary statistics. For details of the method, see ```moloco.pdf```.

### Getting started
```
git clone https://github.com/jimmyzliu/moloco
```
### Requirements
```
python 2.7.3
numpy 1.9.2
scipy 0.14.0
```

### Example
To run on a single region using the example files:
```
python moloco.py \
--stats lung_cancer.15.gwax.assoc.gz,heart_disease.15.gwax.assoc.gz,bronchitis.15.gwax.assoc.gz \
--chr 15 \
--from 78800000 \
--to 79000000 \
--priors 1e-4,1e-5,1e-6 \
--out test
```

To run on multiple regions specified in ```chr15.bed```:
```
python moloco.py \
--stats lung_cancer.15.gwax.assoc.gz,heart_disease.15.gwax.assoc.gz,bronchitis.15.gwax.assoc.gz \
--bed chr15.bed \
--priors 1e-4,1e-5,1e-6 \
--out test
```

The input GWAS summary files require the following columns:
```
CHR SNP BP OR SE P
```

```SE``` is the standard error of the log odds ratio. ```OR``` can be replaced with ```BETA``` for quantitative traits

### Output
The output file will look something like this:
```
config logBF PP
0 0 0.0259778813369
a -1.79404173024 0.00431977677255
b -3.99937376505 0.000476099550502
c -3.48692371991 0.000794789275465
ab -6.35561402013 4.51227573589e-05
ac 3.59547236141 0.946449693824
bc -6.68792696973 3.23648762171e-05
abc -1.79317822772 0.00432350852164
a,b -5.79341550099 7.91690338903e-05
a,c -5.28103706508 0.000132153440041
a,bc -8.48196950332 5.38184490108e-06
b,c -7.48629750719 1.45661920415e-05
ac,b -0.403901403893 0.0173456898908
ab,c -9.84254087183 1.380519483e-06
a,b,c -9.28033921699 2.42216439883e-06
```

The first column indicates the trait configuration. Traits are lettered a, b, c etc. corresponding to the order of the ```--stats``` input. When two (or more) letters are next to each other without a comma, this indicates a configuration where the two traits share a common causal variant. When two traits are separated by a comma, this indicates that there are signals in both traits, but that they do not colocalize. Config "0" denotes the situation where there is no evidence for association with any of the traits.

In the example above: ```a``` means there is only a signal for trait ```a```.

```ab``` means that there is a common signal for traits ```a``` and ```b```.

```a,bc``` means that there are signals for all three traits, one of which is shared between traits ```b``` and ```c``` and another one that is independent in trait ```a```, and so on.

The second column, ```logBF``` gives the log Bayes factor for each configuration. The final column, ```PP``` is that configuration's posterior probability.

### Contact
jliu at nygenome dot org
