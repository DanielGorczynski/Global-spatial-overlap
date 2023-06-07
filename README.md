# Global-spatial-overlap
## Data and code for “Human density modulates spatial associations among tropical forest terrestrial mammal species”

### Data
1. Site.covs.2015.csv – Data file containing covariates used in the global SIF model. Covariates include NDVI over the camera trap sampling area (abbreviated NDVI.2015), average human density in a 5-km buffer surrounding the sampling points in persons per square kilometer (abbreviated Hum.Dens), and the Shannon diversity index of landcover types in the sampling area (abbreviated Hab.div.2015).
2. Full trait list bi 1 kg.csv- Trait data for all species used in this analysis. Variables include average body mass in kilograms, binomial (1 for present, 0 for absent) trait data related to diet, sociality, scan, scansoriality, and activity period, and average litter size in number of individuals. 
3. SIF.models.2015a.csv- SIF estimates from model outputs of two-species occupancy models, used in global analysis
4. Site Trait and Occ data.zip – Zip file containing folder with species lists for each protected area

### Code
1. Co-occurence modeling code GCB.R- Main R file with code necessary to run two-species occupancy models, global SIF model, and create all figures presented in the manuscript
2. Global regression sd cutoff GCB.R- Supplementary R file with code to run filtered global model, only including species pairs with low variance in their SIF estimates. Results presented in Fig. S2

