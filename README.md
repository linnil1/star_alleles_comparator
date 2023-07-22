# Star Allele Comparator (star_allele_comp)
The comparator can compare HLA or KIR alleles between cohort


## Install
``` bash
pip install git+https://github.com/linnil1/star_alleles_comparator
```


## Usage

### 1. Using command line

``` bash
star_allele_comp --csv hla_result1.csv hla_result2.csv --family hla --save tmp --plot -v
```

The input CSV should adhere to the following format:

Format 1: Separate Columns for Alleles

For this format, each allele is represented in separate columns:

``` csv
method,id,allele1,allele2,allele3,allele4
method1,id1,"A*01:02:03:04","A*01:02","B*01:01:01:01"
method1,id2,"A*01:02:03:04","A*01:02","B*01:01:01:01"
method1,id2,"C*03","C*03:03"
method2,id1,"A*01:02:03:04","A*01:02:03","B*01:02:02:01","B*04:01:02"
method2,id2,"A*01:02:03:04","A*01:02:03","B*01:02:02:01","B*04:01:02"
method2,id2,"C*03:03", "C*03:02"
```

Format 2: Using "alleles" Column with Underscore as Separator

In this format, the `alleles` column contains a single string with alleles separated by underscores:

``` csv
method,id,alleles
method1,id3,"KIR2DL1*0010203_KIR2DL1*001_KIR2DS1*0010101"
method1,id4,"KIR2DL1*0010203_KIR2DL1*00102_KIR2DS1*00101"
method1,id3,"KIR2DL1*03105_KIR2DL1*03:03"
method2,id3,"KIR2DL1*001_KIR2DL1*0030203_KIR2DS1*0010208_KIR2DS1*0040102"
method2,id4,"KIR2DL1*0010203_KIR2DL1*0010203_KIR2DS1*0010202_KIR2DS1*0040302"
method2,id4,"KIR2DL1*00303_KIR2DL1*03002"
```

Rules

* The input CSV should contain an `id` column for the sample ID.
* An optional `method` column can be included to indicate the method used for the sample.  If the method column is not specified, filename will be used.
* Columns starting with `allele` will be used to store the allele for each sample.


### 2. Using Python functions

``` python
from star_allele_comp import compare_method, print_all_summary, plot_summary
cohort = {
    "method1": { "sample_id1": [ "A*01:02:03:04", "A*01:02", "B*01:01:01:01", "B*03:01"] },
    "method2": { "sample_id1": [ "A*01:02:03:04", "A*01:02:03", "B*01:02:02:01", "B*04:01:02"] },
}
ground_truth_method = "method1"
result = compare_method(cohort, ground_truth_method, "hla")
# print result allele by allele
print(result)
# Method method2
# Sample sample_id1
# A*01:02:03:04    =4= A*01:02:03:04
# A*01:02          =2= A*01:02:03
# B*01:01:01:01    =1= B*01:02:02:01
# B*03:01          =0= B*04:01:02
# Note:
# Left hand side is alleles in reference method/cohort
# Right hand side is the allele in another method/cohort
# print summary (i.e. Accuracy vs Resolution, Confusion Matrix)
# details are in star_allele_comp/summary.py
df_cohort = result.to_dataframe()
print_all_summary(df_cohort)
# Accuracy summary
#            Accuracy                                num_match                     num_ref
# Resolution        0     1    2    3    4   FP   FN         0  1  2  3  4  FP  FN       0  1  2  3  4  FP  FN
# method
# method1         1.0  1.00  1.0  1.0  1.0  0.0  0.0         4  4  4  2  2   0   0       4  4  4  2  2   0   0
# method2         1.0  0.75  0.5  0.5  0.5  0.0  0.0         4  3  2  1  1   0   0       4  4  4  2  2   0   0
# Note In the accuracy summary table:
# * num_match represents the number of alleles that match the alleles in the ground truth method under the specific `Resolution`.
# * num_ref indicates the number of reference alleles with resolution >= `Resolution`
# * Accuracy is calculated as the ratio of num_match to num_ref.
# * Accuracy in FP is False Discovery Rate (FDR)
# * Accuracy in FN is False Negative Rate (FNR)
# Confusion matrix (not the same sample)
#             Count
#  match_res      -1  0  1  2  3
#  ref_res
# -1              2  0  0  0  0
#  1              1  1  0  0  0
#  2              1  0  2  6  0
#  3              0  0  0  0  1
#  4              0  0  0  0  1
# Accuracy summary per resolution per gene
#              Accuracy                               num_match                     num_ref
# Resolution          0    1    2    3    4   FP   FN         0  1  2  3  4  FP  FN       0  1  2  3  4  FP  FN
# method  gene
# method1 A         1.0  1.0  1.0  1.0  1.0  0.0  0.0         2  2  2  1  1   0   0       2  2  2  1  1   0   0
#         B         1.0  1.0  1.0  1.0  1.0  0.0  0.0         2  2  2  1  1   0   0       2  2  2  1  1   0   0
# Plot summary (i.e. Accuracy vs Resolution, gene, methods)
figs = plot_summary(df_cohort)
# You can use Dash to show it
from dash import dcc, html, Dash
app = Dash(__name__)
app.layout = html.Div(children=[dcc.Graph(figure=fig) for fig in figs])
app.run(debug=True)
```


## Develop

``` bash
pip install pdoc
pdoc star_allele_comp --docformat google
```