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
python moloco.py --stats lung_cancer.15.gwax.assoc.gz,heart_disease.15.gwax.assoc.gz,bronchitis.15.gwas.assoc.gz --chr 15 --from 78516053 --to 80860978 --priors 1e-4,1e-5,1e-6 --out test.out
```
