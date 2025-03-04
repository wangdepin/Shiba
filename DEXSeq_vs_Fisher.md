# Comparison of DEXSeq, Fisher's Exact Test, and Welch's t-test for Differential Splicing Analysis

## DEXSeq for DSE Analysis

### Statistical Model
- Uses a negative binomial (NB) distribution to model count data
- Accounts for biological variation and dispersion between replicates
- Tests for differential exon usage using generalized linear models

### Implementation Details
- Creates exonic parts for "included" and "excluded" forms of each event
- Processes junction counts into DEXSeq-compatible format
- Uses design formula `~ sample + exon + condition:exon`
- Estimates size factors and dispersion parameters
- Tests for differential exon usage (DEU) with `testForDEU()`
- Calculates per-gene q-values with `perGeneQValue()`

### Statistical Power
- Better for experiments with biological replicates
- Can handle within-group variation
- Requires more samples per condition (ideally 3+)

## Fisher's Exact Test (Shiba's Default)

### Statistical Model
- Creates 2x2 contingency tables of junction read counts
- Compares inclusion versus exclusion junctions between two conditions
- Performs Fisher's exact test on these tables
- Takes maximum p-value from multiple tests
- Corrects for multiple testing using Benjamini-Hochberg

### Implementation Details
- Considers a splicing event significant if:
  - Odds ratios indicate consistent direction of change
  - q-value < FDR threshold
  - Absolute dPSI ≥ threshold

### Statistical Power
- Works with small sample sizes, even with just one sample per condition
- Cannot account for biological variation within conditions
- Less robust to outliers and variation between replicates

## Welch's t-test (Shiba's Optional Method)

### Statistical Model
- Compares the means of PSI values between two groups
- Assumes samples come from normal distributions with unequal variances
- Tests the null hypothesis that two populations have equal means

### Implementation Details
- Applied directly to PSI (Percent Spliced In) values, not junction counts
- Uses SciPy's `stats.ttest_ind()` with `equal_var=False` parameter
- Results in a separate p-value ("p_ttest") for each splicing event

### Statistical Power
- Requires multiple samples per condition to be effective
- Can account for biological variation between replicates
- Less powerful when sample sizes are very small or data is not normally distributed
- Better when direct PSI values (rather than raw counts) are the focus of analysis

## Key Differences

1. **Data Input and Modeling**:
   - DEXSeq: Raw count data with negative binomial distribution
   - Fisher's: Junction counts in 2x2 contingency tables
   - Welch's: PSI values with assumed normal distribution

2. **Handling of Biological Variation**:
   - DEXSeq: Models sample-to-sample dispersion explicitly
   - Fisher's: Cannot model within-condition variation
   - Welch's: Accounts for different variances between conditions

3. **Sample Size Requirements**:
   - DEXSeq: Works best with 3+ replicates per condition
   - Fisher's: Can work with even a single sample per condition
   - Welch's: Requires multiple samples per condition (ideally 3+)

4. **Appropriate Use Cases**:
   - DEXSeq: Datasets with biological replicates, when specificity is critical
   - Fisher's: Limited or no replicates, pilot studies, when sensitivity is more important
   - Welch's: When PSI values are directly comparable and normally distributed

5. **Computational Complexity**:
   - DEXSeq: Most computationally intensive
   - Fisher's: Moderate computation
   - Welch's: Least computationally intensive

The choice between these methods should depend on experimental design (number of replicates), biological question (sensitivity vs. specificity needs), and computational constraints.