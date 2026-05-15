# ML-Geo-Pixel-Poppers

**Automating Petrological Analysis: Mineral Phase Classification from EPMA Chemical Maps using Supervised Machine Learning**

ESS 469/569: Machine Learning in Geosciences, Autumn 2023
University of Washington, Earth and Space Sciences

**Project Lead:** Nicole Aikin (naikin@uw.edu) | [nicole-m-aikin](https://github.com/nicole-m-aikin)
**Contributors:** Aiman Shamsul (aiman188@uw.edu) | [AimanHS](https://github.com/AimanHS), Jonathan Lindenmann (linjon06@uw.edu)

---

## Scientific Motivation

Metamorphic petrology relies on identifying and characterizing mineral assemblages in thin sections using chemical maps collected by electron microprobe analysis (EMPA). The conventional workflow for extracting bulk mineral compositions from these datasets is labor-intensive and subjective: researchers manually compile raw TIFF elemental maps in Photoshop, overlay them in Illustrator, identify mineral phases by eye, and count pixels by hand before entering bulk compositions into thermodynamic modeling software such as Perple_X.

This project proposes an open-source, freely available alternative. We apply supervised pixel classification techniques to EMPA chemical map data from metamorphic thin sections, with the goal of automating mineral identification and bulk composition extraction. Our approach is motivated by Rubo et al. (2019), who demonstrated that random forest and convolutional neural network classifiers could automate mineralogy identification in sedimentary rocks. We extend this framework to metamorphic petrology, where mineral morphologies, chemistries, and structural relations are considerably more complex.

**Research question:** Can machine learning optimize and expand models for efficient bulk composition extraction in metamorphic petrology?

---

## Data

### Raw data
Chemical maps of four petrographic thin sections from two sampling locations in the Grand Canyon Upper Granite Gorge, Arizona. Maps were produced at the UMass Amherst Betis microprobe using electron microprobe X-ray analysis of ten elements: Al, Ca, Ce, Fe, K, Mg, Si, Ti, Y, and Zr. Each thin section is stored as a stack of single-element TIFF files.

| Thin section | Dimensions | Elements |
|---|---|---|
| 84.7-NA-2-1_Cold | 695 x 1152 px | 9 |
| 84.7-NA2-2_Cold | 695 x 1152 px | 9 |
| 78.7-10-1_Hot | 703 x 1100 px | 10 |
| 78.7-10-2_Hot | 575 x 1026 px | 10 |

### Labels
Mineral label masks were created manually in Photoshop by visual identification of mineral grains in RGB composite images. A custom interactive widget (PP1) allows users to stack any three elemental maps into a three-channel RGB image, which is then used as a visual reference for hand-drawing mineral masks. Masks are stored as binary TIFF files (1 = mineral present, 0 = absent) and projected into a single labeled image. Seven mineral phases were labeled for the primary thin section: albite, biotite, garnet, K-feldspar, quartz, sericite, and white mica.

### Labeled dataset
Label masks and raw chemical maps are combined into a labeled feature vector dataset stored in HDF5 format. Each pixel is a measurement instance described by an n-dimensional chemical feature vector, a mineral classification label, and (x, y) spatial coordinates. 57% of the 800,640 total pixels in the primary thin section were labeled.

---

## Project Workflow

### 1. Data preprocessing (`PP1`, `PP2`)
- Interactive widget for RGB composite image creation from elemental TIFF stacks
- Manual mineral mask creation in Photoshop
- HDF5 database construction: TIFF stacks + label masks → labeled pixel feature vectors
- Data cleaning: removal of unlabeled vectors, log transformation, normalization, flattening

### 2. Data characterization (`PP3`)
- Class balance evaluation (log plot of mineral label distribution)
- Per-element histograms within mineral label subsets
- Correlation matrices of elemental abundances within each mineral
- Covariance plots
- Principal component analysis (PCA) to assess feature redundancy and dimensionality

Key finding: No strong intra-mineral elemental correlations were identified, suggesting that compositional variance within mineral types is significant and dimensionality reduction may be difficult. The dataset is information-rich but highly unbalanced across mineral types.

### 3. Random forest classification (`PP6`)
- Supervised classification using scikit-learn RandomForestClassifier
- Training data: manually masked garnet and biotite pixels from thin section 78.7-10-1_Hot
- 75/25% train/test split with MinMaxScaler normalization
- Hyperparameter tuning via both grid search and random search cross-validation (10-fold CV)
  - Grid search parameter space: `n_estimators` [50, 100], `criterion` ['gini', 'entropy'], `bootstrap` [True, False]
  - Random search parameter space: `n_estimators` randint(10, 100), `max_depth` randint(1, 20), `min_samples_split` randint(2, 20), `min_samples_leaf` randint(1, 20)

---

## Results

Random forest classification produced spatially coherent mineral probability maps that visually match the known mineral distributions in the thin section.

**Garnet classification:** >96% accuracy. Grid search cross-validation outperformed random search, producing well-defined garnet probability maps that closely match the hand-drawn garnet mask. Signs of overfitting are present, consistent with a small and compositionally specific training dataset.

**Biotite classification:** ~80% accuracy. Both hyperparameter search approaches produced similar results, suggesting better generalization. The lower accuracy reflects the greater compositional overlap between biotite and other phases.

**Out-of-bag error:** Decreases toward 0 as the number of trees increases, with test accuracy stabilizing at ~97% for garnet. The constant accuracy from early in training is consistent with overfitting on the small dataset.

---

## Repository Structure

```
ML-Geo-Pixel-Poppers/
├── codebase/
│   ├── PP0_overview.ipynb              # Project overview and workflow summary
│   ├── PP1_Creating_labelled_data.ipynb # RGB widget + label mask creation
│   ├── PP2_data_preprocessing.ipynb    # HDF5 database construction
│   ├── PP3_Data_characterization.ipynb # Exploratory data analysis
│   ├── PP6_RandomForestTest.ipynb      # Random forest classification + results
│   └── PP90_parallel_coordinates_plot.ipynb # Supplementary visualization
├── Aikin_Data/
│   └── [thin section folders]/
│       ├── RAW/        # Raw elemental TIFF maps
│       ├── labels/     # Binary mineral mask TIFFs
│       ├── RGB/        # RGB composite images
│       └── single_element/ # Single-element visualizations
├── Presentation_images/
├── environment.yml
├── git-cheat-sheet-education.pdf
└── README.md
```

**Note on file paths:** Notebooks contain hardcoded local file paths from the original development environment. To run the code on your own system, update the file paths in each notebook to point to your local copy of the data directory.

---

## Environment and Dependencies

Install the conda environment from the provided file:

```bash
conda env create -f environment.yml
conda activate mlgeo
```

Key dependencies: `numpy`, `pandas`, `scikit-learn`, `matplotlib`, `h5py`, `PIL`, `tifffile`, `ipywidgets`

---

## Known Limitations and Future Directions

This project is a proof of concept developed within the scope of a graduate course. Known limitations include:

- **Overfitting:** Models trained on small, compositionally specific datasets generalize poorly across thin sections or mineral compositions.
- **Grayscale input:** The current random forest implementation converts RGB composites to grayscale, discarding most of the spectral information. Future work should use the full multi-element feature vectors.
- **No spatial features:** Pixel coordinates are preserved in the dataset but not yet incorporated as model features. Spatial context (grain boundaries, inclusion phases) is important for mineralogy that compositional data alone cannot resolve.

**Suggested next steps:**
- Incorporate spatial features using pixel coordinates and convolutional filters
- Train random forest on full PCA-reduced multi-element feature vectors rather than grayscale RGB
- Implement CNN for grain and subgrain boundary recognition
- Resample unbalanced mineral classes to improve generalization
- Extend to additional thin sections and rock types

---

## Citation and Contact

This project was developed as a final project for ESS 469/569: Machine Learning in Geosciences, Autumn 2023, University of Washington.

For questions about the data, instrumentation specifications, or the broader research context, contact Nicole Aikin at nmaikin27@gmail.com.

Reference: Rubo et al. (2019). Machine learning classification of mineral images in petrographic thin sections. *Journal of Petroleum Science and Engineering.*
