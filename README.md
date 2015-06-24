# EfficientRegionalSeismicHazardAnalysis
Seismic hazard analysis method for spatially distributed infrastructure

This repository includes a set of Matlab and Python programs that can be used to develop a set of probabilistic ground motion maps for a region, together with their associated hazard-consistent annual occurrence probabilities. This output can be used to represent regional hazard in a way that is comprehensive, computationally efficient, replicable, and captures spatial correlation (i.e., is appropriate for loss or performance analysis of spatially distributed infrastructure). The program implements the method recommended in Han and Davidson (2012) with one modification, using input hazard information from the New Zealand National Seismic Hazard Model. 

The code is implemented for Christchurch, New Zealand in particular, so it would have to be modified to apply to another geographic region. In particular, the ground motion prediction equations and other details of the seismic hazard analysis would have to be changed. 

The repository includes:
* A user manual that describes the programs in detail.
* The code and input files for Christchurch.
* Select results for the Christchurch example that go with the results in Appendix C of the user manual.

Disclaimer: This code is provided as a service to help others who may be interested in applying the method in other regions, or improving it in different ways. We do not offer any guarantees or take any responsibility for anyone who may use it.

If you do use the code, please cite it properly and share what you have done and learned. 

Han, Y., and Davidson R. 2012. Probabilistic seismic hazard analysis for spatially distributed infrastructure. Earthquake Engineering and Structural Dynamics 41(15), 2141–2158.
