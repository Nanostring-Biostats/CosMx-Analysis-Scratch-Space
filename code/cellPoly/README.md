# cellPoly
An R package for quickly deriving and plotting cell polygons.

![image](https://user-images.githubusercontent.com/4357938/234404680-3fcb5b3a-d679-4946-b804-380bfc8570b4.png)

### System requirements
- R (>= 3.5.0)
- UNIX, Mac or Windows
- see DESCRIPTION for full dependencies

### Demo
See the "vignettes" folder. Vignettes should run in <2 minutes. 

### Instructions for use

The only input data you need is a transcript data frame, a data frame containing these 3 columns:
"cell_id", "x" and "y". 

The function \code{cellPolys} calculates polygon boundaries for given cell IDs, based on their transcript locations.

The function \code{drawPolys} draws polygon boundaries in the R graphics window (using base graphics). It calcualtes any polygons it needs along the way.

The workflow is designed so you only need to compute polygons for cells you plot. 
You start with an empty \code{polys} object, and as you use drawPolys to draw selected cells,
it computes their polygons and saves them to the ever-growing \code{polys} list.

### Installation
```

```
Installation should take < 1 min on a normal desktop computer. 


