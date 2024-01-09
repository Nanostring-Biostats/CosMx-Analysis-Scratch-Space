## Visualize cellular neighborhood in gallery mode

A complete CosMx dataset will contain cell metadata, morphology/protein images and cell label results of cell segmentation. 
We've created a toolkit for visualizing the neighborhood of query cells in terms of protein staining, cell segmentation border, numeric and categorical metadata. (Note: we are *not* performing cell typing or cell segmetnation here, just drawing boundaries from the existing cell label/segmetnation results.)

You can find the package [here](code/NeighVizGallery). See the corresponding tutorial (code/NeighVizGallery/vignettes/tutorial.R) inside the package for more details. 

Plotting morphology images and cell borers of query cells's neighborhood 
![image](https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/assets/62775692/c3011b26-2bab-4005-bb75-da868199be8c)
Plotting numeric and categorical metadata of query cells' neighborhood
![image](https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/assets/62775692/4f43e0be-88f0-4a3b-9036-92e31b61d5aa)
