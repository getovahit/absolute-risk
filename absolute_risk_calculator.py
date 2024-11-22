# Import required libraries
import numpy as np
from scipy.stats import norm
from typing import List, Dict, Optional, Union
import pandas as pd
import json

class PRSRiskCalculator:
    """
    A class to calculate absolute risk from PRS Z-scores.
    
    This calculator can use three different methods to estimate PRS predictive power:
    1. Direct R² from validation studies
    2. Converted from AUC
    3. Theoretical calculation from PRS weights and allele frequencies
    
    The calculator uses the liability threshold model to convert Z-scores to absolute risks,
    accounting for disease prevalence and PRS predictive power.
    """
    
    def __init__(
        self,
        prevalence: float,
        r_squared: Optional[float] = None,
        auc: Optional[float] = None,
        prs_file: Optional[str] = None
    ):
        """
        Initialize the risk calculator with one of three methods.
        
        Args:
            prevalence: Disease prevalence in the population (0 to 1)
            r_squared: Optional liability scale R² from validation study
            auc: Optional AUC from validation study
            prs_file: Optional path to PRS weights file (tab-separated)
                     Must contain columns: VARIANT, BETA, AF
        
        Raises:
            ValueError: If inputs are invalid or missing required information
        """
        # Validate prevalence
        if not 0 < prevalence < 1:
            raise ValueError("Prevalence must be between 0 and 1")
            
        self.prevalence = prevalence
        
        # Determine which method to use for predictive power
        if r_squared is not None:
            # Method 1: Direct R²
            if not 0 <= r_squared <= 1:
                raise ValueError("R² must be between 0 and 1")
            self.r_squared = r_squared
            self.method_used = "empirical_r2"
            
        elif auc is not None:
            # Method 2: Convert from AUC
            if not 0.5 <= auc <= 1:
                raise ValueError("AUC must be between 0.5 and 1")
            self.r_squared = self._auc_to_r2(auc)
            self.method_used = "converted_auc"
            
        elif prs_file is not None:
            # Method 3: Calculate from PRS weights
            # Read and validate the weights file
            weights_df = pd.read_csv(prs_file, sep='\t')
            required_columns = {'VARIANT', 'BETA', 'AF'}
            if not required_columns.issubset(weights_df.columns):
                raise ValueError(f"PRS file must contain columns: {required_columns}")
            
            # Validate allele frequencies
            if not all((weights_df['AF'] >= 0) & (weights_df['AF'] <= 1)):
                raise ValueError("All allele frequencies must be between 0 and 1")
            
            # Calculate theoretical R²
            self.r_squared = self._calculate_theoretical_r2(
                weights_df['BETA'].values,
                weights_df['AF'].values
            )
            self.method_used = "theoretical"
            
        else:
            raise ValueError("Must provide either R², AUC, or PRS weights file")
        
        # Calculate corresponding AUC from R²
        self.auc = self._r2_to_auc(self.r_squared)
        
    def _auc_to_r2(self, auc: float) -> float:
        """
        Convert AUC to liability scale R².
        
        Uses the relationship between AUC and R² derived from the
        liability threshold model.
        
        Args:
            auc: Area Under the ROC Curve (0.5 to 1)
            
        Returns:
            Equivalent liability scale R²
        """
        z = norm.ppf(auc)
        return 2 * (z ** 2) / (np.pi + z ** 2)
    
    def _r2_to_auc(self, r2: float) -> float:
        """
        Convert R² to AUC.
        
        Inverse of the AUC to R² conversion.
        
        Args:
            r2: Liability scale R² (0 to 1)
            
        Returns:
            Equivalent AUC
        """
        return norm.cdf(np.sqrt(r2 / (2 - r2)))
    
    def _calculate_theoretical_r2(self, betas: np.ndarray, allele_freqs: np.ndarray) -> float:
        """
        Calculate theoretical R² from PRS weights and allele frequencies.
        
        Uses the variance explained by each variant:
        - Each variant contributes β² * 2AF(1-AF)
        - Sum is transformed to liability scale using prevalence
        
        Args:
            betas: Effect sizes (log odds ratios)
            allele_freqs: Effect allele frequencies
            
        Returns:
            Theoretical liability scale R²
        """
        # Calculate variance explained by each variant
        variant_contributions = (betas ** 2) * (2 * allele_freqs * (1 - allele_freqs))
        total_variance = np.sum(variant_contributions)
        
        # Transform to liability scale using disease prevalence
        k = self.prevalence
        i = norm.pdf(norm.ppf(1-k))
        c = k * (1-k) / (i ** 2)
        
        return total_variance * c

    def calculate_risk(self, z_scores: Union[float, List[float], np.ndarray]) -> pd.DataFrame:
        """
        Calculate absolute risk from PRS Z-score(s).
        
        Uses the liability threshold model:
        1. Adjusts Z-scores by √R² to account for predictive power
        2. Finds threshold based on prevalence
        3. Calculates probability of exceeding threshold
        4. Provides confidence intervals based on Z-score uncertainty
        
        Args:
            z_scores: Single Z-score or list/array of Z-scores
            
        Returns:
            DataFrame with columns:
            - raw_z: Original Z-score
            - adjusted_z: Z-score adjusted by √R²
            - absolute_risk: Probability of developing condition
            - risk_ci_lower: Lower bound of 95% CI
            - risk_ci_upper: Upper bound of 95% CI
            - relative_risk: Risk relative to population prevalence
            - odds_ratio: Odds ratio relative to population
        """
        # Convert input to numpy array
        if isinstance(z_scores, (float, int)):
            z_scores = [z_scores]
        z_scores = np.array(z_scores)
        
        # Calculate liability threshold from prevalence
        liability_threshold = norm.ppf(1 - self.prevalence)
        
        # Adjust Z-scores by √R² to account for predictive power
        adjusted_z = z_scores * np.sqrt(self.r_squared)
        
        # Calculate absolute risks
        absolute_risk = 1 - norm.cdf(liability_threshold - adjusted_z)
        
        # Calculate confidence intervals (95% CI)
        # Standard error of Z-score is 1 by definition
        z_se = np.ones_like(z_scores)
        adjusted_z_se = z_se * np.sqrt(self.r_squared)
        
        z_lower = adjusted_z - 1.96 * adjusted_z_se
        z_upper = adjusted_z + 1.96 * adjusted_z_se
        
        risk_lower = 1 - norm.cdf(liability_threshold - z_lower)
        risk_upper = 1 - norm.cdf(liability_threshold - z_upper)
        
        # Calculate relative risk and odds ratio
        relative_risk = absolute_risk / self.prevalence
        odds_ratio = (absolute_risk / (1 - absolute_risk)) / (self.prevalence / (1 - self.prevalence))
        
        # Create results DataFrame
        results = pd.DataFrame({
            'raw_z': z_scores,
            'adjusted_z': adjusted_z,
            'absolute_risk': absolute_risk,
            'risk_ci_lower': risk_lower,
            'risk_ci_upper': risk_upper,
            'relative_risk': relative_risk,
            'odds_ratio': odds_ratio
        })
        
        return results
    
    def get_model_info(self) -> Dict[str, float]:
        """
        Get summary of model parameters.
        
        Returns:
            Dictionary containing:
            - r_squared: Liability scale R²
            - auc: Area Under the ROC Curve
            - prevalence: Disease prevalence
            - method_used: Which method was used to get R²
        """
        return {
            'r_squared': self.r_squared,
            'auc': self.auc,
            'prevalence': self.prevalence,
            'method_used': self.method_used
        }
    
    def save_model_info(self, output_file: str):
        """
        Save model information to JSON file.
        
        Args:
            output_file: Path to save JSON file
        """
        with open(output_file, 'w') as f:
            json.dump(self.get_model_info(), f, indent=4)

def main():
    """Example usage demonstrating all three initialization methods."""
    # Example data
    z_scores = [-2.0, -1.0, 0.0, 1.0, 2.0]
    prevalence = 0.05
    
    # Example 1: Using empirical R²
    print("\n1. Using empirical R²:")
    calc_r2 = PRSRiskCalculator(prevalence=prevalence, r_squared=0.1)
    results_r2 = calc_r2.calculate_risk(z_scores)
    print("\nModel Info:")
    print(json.dumps(calc_r2.get_model_info(), indent=2))
    print("\nRisk Calculations:")
    print(results_r2.round(3))
    
    # Example 2: Using AUC
    print("\n2. Using AUC:")
    calc_auc = PRSRiskCalculator(prevalence=prevalence, auc=0.65)
    results_auc = calc_auc.calculate_risk(z_scores)
    print("\nModel Info:")
    print(json.dumps(calc_auc.get_model_info(), indent=2))
    print("\nRisk Calculations:")
    print(results_auc.round(3))
    
    # Example 3: Show PRS weights file format
    print("\n3. Using PRS weights file (example format):")
    print("VARIANT\tBETA\tAF")
    print("rs123\t0.02\t0.3")
    print("rs456\t0.01\t0.5")
    print("rs789\t0.03\t0.1")
    
if __name__ == "__main__":
    main()
