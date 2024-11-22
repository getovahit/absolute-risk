# Mathematical Framework for PRS Absolute Risk Calculation

## Introduction
This document provides a comprehensive mathematical analysis of our methodology for converting Polygenic Risk Scores (PRS) to absolute risk estimates. We present the theoretical foundations, assumptions, derivations, and justifications for our approach.

## 1. Core Mathematical Framework

### 1.1 The Liability Threshold Model
The foundation of our approach rests on the liability threshold model, which assumes:
1. An underlying continuous liability score $$L$$ that follows a standard normal distribution in the population
2. Disease occurs when this liability exceeds a threshold $$T$$
3. The threshold $$T$$ is determined by the population prevalence $$k$$

Mathematically:
$$L \sim N(0,1)$$
$$P(\text{Disease}) = P(L > T) = k$$
$$T = \Phi^{-1}(1-k)$$

where $$\Phi^{-1}$$ is the inverse standard normal cumulative distribution function.

### 1.2 PRS Z-scores
The PRS Z-score represents how many standard deviations an individual's polygenic score deviates from the population mean. For an individual with Z-score $$z$$:
$$\text{PRS}_z \sim N(z, 1)$$

### 1.3 Predictive Power and Scale Adjustment
The key insight is that the raw Z-score must be adjusted by the predictive power of the PRS. We use $$R^2$$ (for continuous traits) or convert AUC to an equivalent $$R^2$$ (for binary traits).

For a PRS explaining variance $$R^2$$, the adjusted Z-score is:
$$z_{adj} = z\sqrt{R^2}$$

Justification:
- $$R^2$$ represents the proportion of variance explained
- $$\sqrt{R^2}$$ represents the correlation coefficient
- Multiplying by $$\sqrt{R^2}$$ appropriately scales the effect size

## 2. Converting to Absolute Risk

### 2.1 Binary Traits
For binary traits, the absolute risk calculation uses the liability threshold model:

$$P(\text{Disease}|z) = 1 - \Phi(\frac{T - z_{adj}}{\sqrt{1-R^2}})$$

where:
- $$T$$ is the liability threshold determined by prevalence
- $$z_{adj}$$ is the adjusted Z-score
- $$\sqrt{1-R^2}$$ accounts for residual variance

### 2.2 Confidence Intervals
We incorporate uncertainty in our predictions through confidence intervals that account for:
1. Inherent uncertainty in Z-score (sampling variation)
2. Uncertainty in $$R^2$$ estimation

For a given confidence level $$\alpha$$ (e.g., 95%):
$$z_{lower} = z_{adj} - z_{\alpha/2}\sigma_{adj}$$
$$z_{upper} = z_{adj} + z_{\alpha/2}\sigma_{adj}$$

where $$\sigma_{adj} = \sqrt{R^2}$$ (since Z-scores have unit variance by definition)

## 3. Estimation of Predictive Power

### 3.1 Direct $$R^2$$ or AUC
When empirical measures are available:
- Use $$R^2$$ directly for continuous traits
- For binary traits, convert AUC to $$R^2$$ using:
  $$R^2 = \frac{2(z_{AUC})^2}{\pi + (z_{AUC})^2}$$
  where $$z_{AUC} = \Phi^{-1}(AUC)$$

### 3.2 Theoretical Calculation
When empirical measures aren't available, we estimate $$R^2$$ from PRS weights and allele frequencies:

$$R^2_{theoretical} = \sum_{i=1}^n \beta_i^2(2AF_i(1-AF_i))$$

where:
- $$\beta_i$$ is the effect size for variant $$i$$
- $$AF_i$$ is the effect allele frequency
- The sum is transformed to liability scale for binary traits

## 4. Mathematical Justifications

### 4.1 Why Scale by $$\sqrt{R^2}$$?
1. Consider the linear model: $$y = \beta x + \epsilon$$
2. When $$x$$ is standardized (Z-score), $$\beta = r$$ where $$r$$ is correlation
3. $$R^2 = r^2$$ is the squared correlation
4. Therefore, $$\sqrt{R^2}$$ is the appropriate scaling factor

### 4.2 Liability Scale Transformation
For binary traits, we transform observed scale $$R^2$$ to liability scale using:

$$R^2_L = R^2_O \frac{k(1-k)}{z^2}$$

where:
- $$k$$ is prevalence
- $$z$$ is the normal density at threshold $$T$$
- This accounts for the non-linear relationship between liability and observed scales

## 5. Comparison with Alternative Approaches

### 5.1 Advantages of Our Approach
1. Computational Efficiency:
   - Direct calculations without numerical integration
   - No need for bivariate normal calculations

2. Individual-Level Predictions:
   - Works directly with Z-scores
   - Provides personalized confidence intervals

3. Flexibility:
   - Multiple methods for obtaining $$R^2$$
   - Handles both binary and continuous traits

### 5.2 Potential Limitations
1. Assumes normality of genetic effects
2. May not capture some complex correlation structures
3. Relies on accuracy of $$R^2$$ estimates

## 6. Extensions and Future Directions

### 6.1 Non-Normal Traits
For traits that deviate from normality:
1. Consider quantile normalization
2. Use empirical distribution functions
3. Apply appropriate transformations

### 6.2 Population Stratification
Account for population differences through:
1. Ancestry-specific allele frequencies
2. Population-specific prevalence estimates
3. Stratified $$R^2$$ calculations

## References
1. Dempster & Lerner (1950): Heritability of Threshold Characters
2. Lee et al. (2011): Estimating Missing Heritability
3. Wray et al. (2010): The Genetic Interpretation of AUC
