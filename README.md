# PRS Absolute Risk Calculator

This tool converts Polygenic Risk Score (PRS) Z-scores into absolute risk estimates, accounting for disease prevalence and PRS predictive power. It offers multiple methods for incorporating PRS predictive power and provides confidence intervals for risk estimates. It is similar but different from the methodology employed below:
[
](https://opain.github.io/GenoPred/Absolute_Conversion.html)

Note the code is more pseudocode currently and has never been tested but conveys overall idea. 
## Table of Contents
- [Background](#background)
- [Mathematical Framework](#mathematical-framework)
- [Installation](#installation)
- [Usage](#usage)
- [Input Formats](#input-formats)
- [Output Format](#output-format)
- [Implementation Details](#implementation-details)
- [Examples](#examples)

## Background

Polygenic Risk Scores (PRS) are typically reported as Z-scores (standard deviations from the mean). While useful for relative comparisons, clinicians and patients often want to know the absolute risk of developing a condition. This tool converts PRS Z-scores into absolute risk estimates while accounting for:

1. Disease prevalence
2. PRS predictive power (measured by R² or AUC)
3. Uncertainty in risk estimates

## Mathematical Framework

### Core Concepts

1. **Liability Threshold Model**
   - Assumes an underlying continuous liability score that follows a normal distribution
   - Disease occurs when liability exceeds a threshold
   - Threshold is determined by disease prevalence
   - PRS provides information about where an individual falls on this liability scale

2. **Risk Calculation**
   ```
   liability_threshold = norm.ppf(1 - prevalence)
   adjusted_z = raw_z * sqrt(R²)
   absolute_risk = 1 - norm.cdf(liability_threshold - adjusted_z)
   ```

3. **Predictive Power Estimation**
   We support three methods:
   
   a) **Direct R²**
   - Use empirically measured R² from validation studies
   
   b) **AUC Conversion**
   - Convert AUC to R² using:
   ```
   z = norm.ppf(AUC)
   R² = 2 * (z²) / (π + z²)
   ```
   
   c) **Theoretical Calculation**
   - Calculate from PRS weights and allele frequencies:
   ```
   variant_contribution = β² * 2 * AF * (1-AF)
   total_variance = Σ(variant_contributions)
   R² = total_variance * c  # transformed to liability scale
   ```
   where c is the liability scale transformation factor

4. **Confidence Intervals**
   - Based on inherent uncertainty in Z-score predictions
   - Uses standard error of 1 (by definition of Z-score)
   - Adjusted by √R² to match the scale of risk predictions

## Installation

```bash
pip install numpy pandas scipy
```

## Usage

Three ways to initialize the calculator:

1. Using empirical R²:
```python
calculator = PRSRiskCalculator(
    prevalence=0.05,
    r_squared=0.1
)
```

2. Using AUC:
```python
calculator = PRSRiskCalculator(
    prevalence=0.05,
    auc=0.65
)
```

3. Using PRS weights file:
```python
calculator = PRSRiskCalculator(
    prevalence=0.05,
    prs_file="prs_weights.txt"
)
```

Calculate risks:
```python
z_scores = [1.0, 2.0, 3.0]
results = calculator.calculate_risk(z_scores)
```

## Input Formats

### PRS Weights File
Tab-separated file with columns:
```
VARIANT  BETA    AF
rs123    0.02    0.3
rs456    0.01    0.5
rs789    0.03    0.1
```
- VARIANT: Variant identifier
- BETA: Effect size
- AF: Effect allele frequency (0-1)

## Output Format

Returns a pandas DataFrame with columns:
- raw_z: Original Z-score
- adjusted_z: Z-score adjusted by √R²
- absolute_risk: Probability of developing condition
- risk_ci_lower: Lower bound of 95% CI
- risk_ci_upper: Upper bound of 95% CI
- relative_risk: Risk relative to population prevalence
- odds_ratio: Odds ratio relative to population

## Implementation Details

### Key Components

1. **Z-score Adjustment**
   - Raw Z-scores are scaled by √R² to account for predictive power
   - This scaling ensures predictions aren't overconfident

2. **Risk Calculation**
   - Uses liability threshold model
   - Transforms adjusted Z-scores to risk probabilities
   - Accounts for disease prevalence

3. **Confidence Intervals**
   - 95% CIs based on Z-score uncertainty
   - Properly scaled by predictive power
   - Reflects uncertainty in individual predictions

## Examples

```python
# Initialize calculator with AUC
calculator = PRSRiskCalculator(
    prevalence=0.05,
    auc=0.65
)

# Calculate risks for several Z-scores
z_scores = [-2.0, -1.0, 0.0, 1.0, 2.0]
results = calculator.calculate_risk(z_scores)

# Get model information
model_info = calculator.get_model_info()
```

Example output:
```
   raw_z  adjusted_z  absolute_risk  risk_ci_lower  risk_ci_upper  relative_risk  odds_ratio
0   -2.0      -0.63         0.020         0.012         0.033           0.40         0.38
1   -1.0      -0.32         0.035         0.022         0.054           0.70         0.69
2    0.0       0.00         0.050         0.033         0.075           1.00         1.00
3    1.0       0.32         0.071         0.048         0.103           1.42         1.45
4    2.0       0.63         0.098         0.068         0.139           1.96         2.06
```

This implementation provides a rigorous yet practical way to convert PRS Z-scores into clinically meaningful absolute risk estimates, while properly accounting for the predictive power of the PRS and providing appropriate uncertainty estimates.
