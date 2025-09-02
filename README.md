The analysis begins with the standard workflow including:

1. Meta data parsing

2. Clustering

3. Dimensionality reduction

4. Integration

5. Automated annotation with single R

6. Manual annotation

Pseudobulk analysis per cell type (‘Top Marker for Cell Typing’):

1. ‘pseudobulk_final.R’

2. ‘pseudobulk_functions.R’

3. SearchLight2 results for each cell type

Output:

1. Contains all cell type pseudobulk matrices and the corresponding DE tables as well as the sample sheet.

2. It also contains the top markers.tsv

Selecting regions of interest:

*Since there were two replicates for each fuso_pos sample, I separated all 28 samples using C1, C2, C3 identifiers in the format - sampleID_C1/2/3):

C1- fuso_neg

C2- fuso_pos

C3-fuso_pos

Files needed:

· centroids1.csv

· selector.py

· draw_circ.py


1. Run the draw_circ.py script.

Since the centroids represent the position of the center of mass of each segmented cell, the x_max and y_max are just the max values of X and Y and the center_x and center_y is computed by calculating the mean of X and Y as the mean represents the location of the center of mass of the entire distribution of cells.

2. The script will automatically save the cells in each circle excluding the previous to separate files per sample, so all that needs to be specified is the sample directory.

Example: C:/Users/AakankshaChoudhary/OneDrive/Desktop/cosmx/INC0046_C2/my_regions
