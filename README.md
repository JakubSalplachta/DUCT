# DUCT
Matlab pipeline for analysis of PV and BD systems
   - input: binary masks for each system, binary masks of main branches
   - output: numerical results of evaluated parameters
   
   - pipeline is based on skeletonization, followed by graph conversion and subsequent analysis
   - skeletozation based on: SKELETON3D by Philip Kollmannsberger (philipk@gmx.net) [1]
   - graph coversion based on: SKEL2GRAPH3D by Philip Kollmannsberger (philipk@gmx.net) [1]
   
The supplementary algorithm to DUCT pipeline was designed to specifically analyze morphological parameters of the BD and PV systems, and is compatible with the 3D binary masks. Using DUCT two separate masks for both BD and PV systems were generated and the analysis was divided in two independent parts. First, the analysis of the entire portal vein and biliary system and second, analysis of the corresponding main branch (= the longest branch) of each system. For detailed analysis and comparison of whole system versus only the main branch, two algorithms were developed. They differ in the input data and the evaluated parameters. For both algorithms, the first step is to create a 3D skeleton of the input binary mask. Any comments, corrections or suggestions are highly welcome.

If you include this in your own work, please cite our original publicaton [2].

 
 Jakub Salplachta (jakub.salplachta@ceitec.vutbr.cz
 
 References: 
[1] Kerschnitzki, Kollmannsberger et al., "Architecture of the osteocyte network correlates with bone material quality." Journal of Bone and Mineral Research, 28(8):1837-1845, 2013.
[2] ...
 
 
