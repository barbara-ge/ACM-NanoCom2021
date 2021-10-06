# ACM-NanoCom2021 - Astrocytic control in in-vitro and simulated neuron-astrocyte networks

In this repository you will find the MATLAB files and datasets used to obtain the results shown in the publication: 
- Barbara Genocchi, Annika Ahtiainen, Michael T. Barros, Jarno M.A. Tanskanen, Jari Hyttinen, and Kerstin Lenk. 2021. Astrocytic control in in vitro and simulated neuron-astrocyte networks. In Proceedings of the Eight Annual ACM International Conference on Nanoscale Computing and Communication (NANOCOM '21). Association for Computing Machinery, New York, NY, USA, Article 21, 1–7. DOI: https://doi.org/10.1145/3477206.3477458
 


In the folder **Experiments** you will find three subfolders:
- **analysis**: Already analysed data.
- **data**: Preprocessed datacells to feed in the code. Please mind that the script *correlation_analysis_experiments.m* might take up to one day to analyse the whole experimental dataset. So if you don't wish to modify any of the parameters, you can directly use the files in *Analysis*.
- **code**: Analysis script and functions and plotting scripts. There are step by step comments in the file where you need to modify the code to adapt it to your machine. 

The data contained in this folder has been obtained from MEA recordings as described in the Methods of the paper.


In the folder **Simulations** you will find two subfolders analogous in content to the previous:
- **code**
- **data**

The data in this folder has been simulated using INEXA (doi: https://doi.org/10.3389/fncom.2019.00092) as described in te Methods of the current paper.


For any further doubt or question, please do not hesitate to contact me at barbara.genocchi@tuni.fi. 
