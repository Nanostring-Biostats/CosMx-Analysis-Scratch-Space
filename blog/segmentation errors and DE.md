# Impact of segmentation errors on DE

Cell segmentation is always imperfect, leaving some cells' expression profiles contaminated with transcripts properly belonging to other cells. 
For many analyses, this is a largely ignorable source of noise. But for differential expression (DE) analyses, it's reliably confounding (both statistically and emotionally).

## The problem:

For an example, see the below cartoon, where a T-cell's expression profile is contaminated with transcripts from the tumor cells surrounding it:

![image](https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/assets/4357938/c7d4e91f-a225-4ed8-a935-f77c0d817670)

Now say we want to compare T-cells in the tumor bed vs. T-cells in the stroma. We'll find that T-cells in the tumor bed are enriched in
genes expressed by cancer cells (e.g. keratins), and similarly, T-cells in the stroma will be enriched in genes expressed by stroma cells 
(e.g. collagens). In practice, spurious findings like these are often the most significant genes emerging from a DE analysis. 

## How to think about segmentation error:

Segmentation error can be considered a missing term in a differential expression model. 
Say you're answering the above T-cell question by fitting the below model on T-cells:

*E(observed KRT9 expression) = B0 + B1 * (in_tumor)* 

When you do this, you're omitting an important term for contamination:

*E(observed KRT9 expression) = B0 + B1 * (in_tumor) + (**N_contaminating_transcripts**)*  

Your model is underspecified, so its results are biased. As you measure more cells and gain statistical power, you only gain more confidence in your biased results. 

## Countermeasures
