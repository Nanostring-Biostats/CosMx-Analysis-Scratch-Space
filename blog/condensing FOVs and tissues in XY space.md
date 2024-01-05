## Condensing cells in xy space for better plotting 

Minimizing whitespace while plotting cells in xy space is a constant challenge.
A single tissue will often have discontinuous FOVs, and aligning multiple tissues
in a sensible way can be onerous. 

Here, for example, are FOVs collected from core needle biopsies, where the cells can barely be seen against the vast expanse of white space. 

![image](https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/assets/4357938/8d7c1716-c3a6-4d85-9f5f-79cc9374fe0f)

As a partial solution, see the function consenseXY(), providing here: https://github.com/patrickjdanaher/Cosmx-Analysis-Scratch-Space/tree/main/code/condensing%20xy%20space

The main wrapper function contains an algorithm for pulling together FOVs from the same tissue, and an algorithm for tiling tissues across a plot. 

Here's a toy example of FOV groups from two tissues before and after the algorithm (color denotes tissue ID):

![image](https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/assets/4357938/3f897cc4-6b0c-4583-a825-d66b4aecec3a)

It's not perfect, but it's an improvement on the original spacing with no thought or manual labor. 

Warning: the FOV condensing code is inefficiently written and takes longer than it should, though it's still faster than working by hand. 
