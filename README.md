# microsimulation-cerebralpalsy
############################################################################################
################# Microsimulation modeling using R: a tutorial #### 2025 ###################
############################################################################################
# This code forms the basis for the microsimulation model of the article: 
#
# Diana Marcela Nova Díaz, Sergio Aguilera Albesa, Eduardo Sánchez Iriso. 
# Microsimulation modeling for health decision sciences using R: A tutorial.

# Cost-Effectiveness of Complementary and Alternative Therapies for Children with Cerebral Palsy: 
# Evidence from Real-World Data and Microsimulation Modelling

# Please cite the article when using this code

This repository contains the R scripts used for the cost-effectiveness analysis of complementary and alternative therapies (CATs) in children with cerebral palsy (CP). The analysis combines real-world clinical data with short-term statistical modelling and long-term individual-level microsimulation.

## Project Description

- **Population**: 148 children with cerebral palsy stratified by Gross Motor Function Classification System (GMFCS) levels.
- **Short-term analysis**: Seemingly unrelated regression equations (SURE) to estimate incremental costs and QALYs based on EQ-5D-Y.
- **Long-term modelling**: First-order individual-level microsimulation projecting costs and QALYs over a 30-year time horizon.
- **Perspective**: Spanish public healthcare system.
- **Discount rate**: 3% annually for costs and QALYs.

## Files Included

- `MICROSIM_PHE.R`: Microsimulation model structure and parameters.
- `Microsim_AlterTherapies.R`: Analysis of alternative therapies (e.g., homeopathy, Peto Method).
- `Microsim_CompleTherapies.R`: Analysis of complementary therapies (e.g., Speech therapy, physiotherapy).

## Instructions

Each script can be run independently in RStudio. It is recommended to have the following R packages installed: `dplyr`, `data.table`, `ggplot2`. Scripts simulate patient-level transitions and calculate incremental cost-effectiveness ratios (ICERs) and cost-effectiveness acceptability curves (CEACs).

For any questions, please contact:  
**Diana Marcela Nova Díaz**  
Email: dianamarcela.nova@unavarra.es







