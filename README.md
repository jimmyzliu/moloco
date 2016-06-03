# moloco
A tool for performing multiple-trait colocalization test using GWAS summary statistics

To get started:
```
git clone https://github.com/jimmyzliu/moloco
```
Requirements:
```
python 2.7.3
scipy 0.14.0
```

To run using the example files:
```
python moloco.py --stats lung_cancer.15.gwax.assoc.gz,heart_disease.15.gwax.assoc.gz,bronchitis.15.gwax.assoc.gz --chr 15 --from 78516053 --to 80860978 --priors 1e-4,1e-5,1e-6 --out test.out
```

The output file will look something like this:
```
config logBF PP
a 3.14680336111 0.0149183474663
b -0.657338372998 0.000332355564714
c 2.02583633886 0.00486285058827
ab -1.62508527007 0.000126274546479
ac 6.69962714556 0.520825791988
bc -2.88627046433 3.57758570759e-05
abc 2.09524681533 0.00521237331747
a,b 2.48946498795 0.00773111778662
a,c 5.17263965393 0.11311761542
a,bc 0.260532890523 0.000832203196364
b,c 1.36849796572 0.00252006938184
ac,b 6.04228877254 0.269906942029
ab,c 0.400751063351 0.000957470403343
a,b,c 4.51530132698 0.0586208124549
```

The first column indicates the trait configuration. Traits are lettered a, b, c etc. corresponding to the order of the ```--stats``` input. When two (or more) letters are next to each other without a comma, this indicates a configuration where the two traits share a common causal variant. When two traits are separated by a comma, this indicates that there are signals in both traits, but that they do not colocalize.

In the example above: ```a``` means there is only a signal for trait ```a```.

```ab``` means that there is a common signal for traits ```a``` and ```b```.

```a,bc``` means that there are signals for all three traits, one of which is shared between traits ```b``` and ```c``` and another one that is independent in trait ```a```, and so on.

The second column, ```logBF``` gives the log Bayes factor for each configuration. The final column, ```PP``` is that configuration's posterior probability.
