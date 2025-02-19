# chase-clustering

This repository contains data and R scripts for clustering analyses of Quaternary mammal communities, focusing on both a novel iterative “chase” algorithm and a spatially constrained Ward’s hierarchical approach. With these scripts, you can replicate the results of:

Brook, B.W., Alroy, J. _et al._ (2025). Title. _EcoEvoRxiv_.** (placeholder; forthcoming)

## Contents

- **data/**  
  CSV files storing the fossil mammal occurrence data and site coordinates:
  - `pleistocene_all.csv`  
  - `holocene_all.csv`  
  - `all_site_coords.csv`  
  - `pleistocene_paired.csv`  
  - `holocene_paired.csv`  
  - `paired_site_coords.csv`  
  
- **src/**  
  Three main R scripts that conduct all analyses described in the paper:
  1. `qe_biogeography_clustchase.r` – Implements the chase clustering algorithm, including data filtering and iterative improvement of cluster assignments based purely on compositional similarity.  
  2. `qe_biogeography_clustgeo.r` – Applies Ward’s hierarchical clustering with spatial constraints (via ClustGeo) to illustrate the role of geography in shaping community structure.  
  3. `qe_turnover_paired_sites.r` – Calculates site-level turnover for paired localities that have both Late Pleistocene and Holocene faunal data, displaying shifts in species composition through time.

## Requirements

- **R version:** 4.0 or higher (tested on 4.4.2)  
- **R packages:**  
  - `dplyr` (data wrangling)  
  - `magrittr` (pipe operator `%>%`)  
  - `tibble` (row/column manipulation)  
  - `ggplot2` (plotting)  
  - `sf` or `sp`/`rgdal` (for mapping)  
  - `maps` or equivalent (for quick world map polygons)  
  - `ClustGeo` (for Ward’s hierarchical + spatial constraints)  
  - Additional packages used may be loaded via `import::from` statements in the scripts

Make sure these packages are installed, for example:
```r
install.packages(c("dplyr", "magrittr", "tibble", "ggplot2", "sf", "maps", "ClustGeo"))
```
(This list may not be exhaustive; check the script headers for any others.)

## Usage

1. **Clone or download** this repository to your local machine.  
2. **Open R** (or RStudio) and set your working directory to the cloned folder:
   ```r
   setwd("path/to/chase-clustering")  # adapt as needed
   ```
3. **Run the scripts** in the `/src` folder in whichever order suits your interest:
   - `qe_biogeography_clustchase.r` for the chase clustering method.  
   - `qe_biogeography_clustgeo.r` for Ward’s constrained clustering.  
   - `qe_turnover_paired_sites.r` for paired-site turnover calculations.  
4. **Inspect outputs**: each script produces console output, plots (in the R graphics window), and optionally saves figures to disk (`.png` or `.pdf`).  

You may edit script parameters (e.g., `max_clusters`, `min_sp`, etc.) to explore different settings.

## Data Description

All CSV files in `data/` are read by the R scripts to replicate the analyses in the forthcoming paper:

- `pleistocene_all.csv` and `holocene_all.csv`: large presence/absence (or counts) matrices (rows = species, columns = sites).  
- `all_site_coords.csv`: coordinates (longitude, latitude) and other metadata for each site in the full dataset.  
- `pleistocene_paired.csv`, `holocene_paired.csv`, `paired_site_coords.csv`: smaller subset of 34 “paired” sites that have data in both periods, used for turnover analyses.

## Replicating the Study

1. **Data Filtering**: scripts automatically exclude sites below a certain species threshold (usually 5) to reduce noise.  
2. **Clustering**: 
   - _Chase_ method (purely compositional) is iterative and random, so results can vary slightly run to run.  
   - _Ward’s_ method merges sites by minimizing a composite distance of compositional and geographic information (from ClustGeo).  
3. **Turnover**: For the paired sites, a straightforward index of compositional change between Late Pleistocene and Holocene is computed, and a map is generated to illustrate which sites underwent the greatest faunal shifts.

Refer to the inline comments in each script for detailed steps and explanations.

## Citation

If you use this repository, code, or data in your own work, please cite:

- The forthcoming paper:  
  Brook, B.W., Alroy, J. _et al._ (2025). Title. _EcoEvoRxiv_. (Exact details will be added here once published.)

- This repository:  
  ```
  @misc{brook2025chaseclustering,
    author = {Brook, B.W. and others},
    title = {chase-clustering},
    howpublished = {\url{https://github.com/bwbrook/chase-clustering}},
    year = {2025}
  }
  ```

## License

This project is released under the MIT License. See [LICENSE](LICENSE) for details. Data and code are provided “as is,” without warranty of any kind.

## Contributing

Pull requests and issue reports are welcome. For major changes, please open an issue first to discuss proposed modifications.

## Contact

For questions or collaboration inquiries, please contact Barry Brook at barry.brook@utas.edu.au.

---  

Please explore, adapt, and share. We hope these methods facilitate further research into Quaternary mammal community patterns and beyond.
