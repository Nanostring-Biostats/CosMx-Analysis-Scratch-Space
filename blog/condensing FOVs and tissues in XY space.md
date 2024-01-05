## Condensing cells in xy space for better plotting 

Minimizing whitespace while plotting cells in xy space is a constant challenge.
A single tissue will often have discontinuous FOVs, and aligning multiple tissues
in a sensible way can be onerous. 

Here, for example, are FOVs collected from core needle biopsies, where the cells can barely be seen against the vast expanse of white space. 

![image](https://github.com/patrickjdanaher/Cosmx-Analysis-Scratch-Space/assets/4357938/e03a97fa-f819-4846-a5f0-45a625e66cfb)

As a partial solution, see the function consenseXY(), providing here: https://github.com/patrickjdanaher/Cosmx-Analysis-Scratch-Space/tree/main/code/condensing%20xy%20space

The main wrapper function contains an algorithm for pulling together FOVs from the same tissue, and an algorithm for tiling tissues across a plot. 

Here's a toy example of FOV groups from two tissues before and after the algorithm:

![image](https://github.com/patrickjdanaher/Cosmx-Analysis-Scratch-Space/assets/4357938/2e839ef1-e155-4d78-9768-5f1d3163e380)

It's not perfect, but it's an improvement on the original spacing with no thought or manual labor. 

Warning: the FOV condensing code is inefficiently written and takes longer than it should, though it's still faster than working by hand. 
