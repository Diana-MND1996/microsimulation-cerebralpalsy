**Microsimulation modeling in R for pediatric patients with cerebral palsy**

This GitHub repository provides the code explained in the paper: Nova Díaz Diana Marcela, Aguilar Albesa Sergio, Eduardo Sánchez-Iriso.
Cost-Effectiveness Analysis of Complementary and Alternative Therapies for Children with Cerebral Palsy: Evidence from Real-World Data and Microsimulation Modelling


This repository includes the R scripts used to perform a cost-effectiveness analysis of complementary and alternative therapies (CATs) in children with cerebral palsy (CP). The analysis uses real-world clinical data collected from our study population. It applies a microsimulation model at the individual level to evaluate the long-term outcomes of these therapies in patients with CP.

  ***Project Description***

- **Population**: 148 children with cerebral palsy stratified by Gross Motor Function Classification System (GMFCS) levels.
- **Long-term modelling**: First-order individual-level microsimulation projecting costs and QALYs over a 30-year time horizon.
- **Perspective**: Spanish public healthcare system.
- **Discount rate**: 3% annually for costs and QALYs.

  ***Sections***

- [MICROSIM_PHE.R](MICROSIM_PHE.R): Microsimulation model structure and parameters. Cost-effectiveness analysis of complementary therapies, alternative therapies, standard treatment versus no treatment.
- [Microsim_AlterTherapies.R](Microsim_AlterTherapies.R): Cost-effectiveness analysis of each of the alternative therapies versus standard treatment.
- [Microsim_CompleTherapies.R](Microsim_CompleTherapies.R): Cost-effectiveness analysis of each of the complementary therapies versus standard treatment.

***List of Contributors***:
- Diana Marcela Nova Díaz
- Eduardo Sánchez Iriso

For any questions, please contact:  
**Diana Marcela Nova Díaz**  
Email: dianamarcela.nova@unavarra.es
