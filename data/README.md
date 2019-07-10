### Data files supporting Vesk & Yen (submitted) How many sprouts and how deep? Flexible modelling of multi-species experimental disturbances. 


#### Overview
Three data files are provided: `compiled-sprout-depths.csv`, `survival-data.csv`, and `sprouts-trait-data.csv`. These files contain the raw resprouting depths, binary resprouting success/failure, and species-level trait data, respectively. Each file is described below.


##### Compiled sprout depths
The `compiled-sprout-depths.csv` file has one row per replicate (individual within treatment within species). The file is in wide format, with 5 columns containing identifying information and 88 columns containing sprout depths measured in each replicate. The columns are:
- `spp`: species code (defined in `sprouts-trait-data.csv`).
- `plant`: plant ID within each treatment.
- `coll`: year of measurement, with month in years with multiple surveys.
- `yr`: binary year of measurement (2000 = 0; 2001 = 1).
- `treat`: treatment (clip = c, clip and burn = b).
- `depths.mm.`: all remaining columns are depths of sprouts in mm, with values for measured sprouts and NA in all columns to the right of the last measured sprout.


##### Resprouting success
The `survival-data.csv` file has one row per replicate (individual within treatment within species). The file is in long format, with 4 columns containing identifying information and a binary variable indicating resprouting success. The columns are:
- `spp`: as for compiled sprout depths.
- `treat`: as for compiled sprout depths.
- `yr`: as for compiled sprout depths.
- `survival`: binary variable indicating whether an individual had any sprouts (1) or no sprouts (0).


##### Trait information
The `sprouts-trait-data.csv` file has one row per species with 9 columns containing information on each species. The columns are:
- `CODE`: identical to the `spp` column in compiled sprout depths data.
- `MAXHT`: maximum species heights in mm.
- `GFORM`: species growth form (one of grass gr, chenopod subshrub ssh, shrub s, forb f, or tree tr).
- `STEMCATS`: ordinal variable, reflecting the number of stems of each species, measured on a half-log10 scale.
- `SLA`: specific leaf area for each species measured in mm2 / g.
- `LWC`: (not used in analyses) the leaf water content for each species.
- `BA`: (not used in analyses) the basal area of each species in mm.
- `CA`: (not used in analyses) the canopy area of each species in mm.
- `SPECIES`: the full species name.

