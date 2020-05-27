# LEEPO
Liposome Encapsulation Efficiency Predictor &amp; Optimizer

This program was made based on the following research work: 

Xu, Xiaoming, Mansoor A. Khan, and Diane J. Burgess. "Predicting hydrophilic 
drug encapsulation inside unilamellar liposomes." 
International journal of pharmaceutics 423.2 (2012): 410-418.

DOI : https://doi.org/10.1016/j.ijpharm.2011.12.019



Python code to calculate Encapuslation efficiecny of Liposomes

(L)iposome

    (E)ncapsulation
    
        (E)fficiency
        
            (P)redictor &
            
                (O)ptimizer
            
Program Usage : python LEEPO.py mu sigma, sizerange ... 

Example command : python LEEPO.py 190.0 38.0 1000 0.23 5.2 0.45 140 5 60.04

Enter the system parameters in the order of:

1. Mean                         (μ - nm {float})
2. Standard deviation           (σ - nm {float})
3. Size range                   (eg: for 1-1000 nm enter "1000")
4. Anticipated defect           (% {float})
5. Bilayer thickness           (d - nm)
6. Molecular area of lipid      (a - nm^2)
7. Lipid molar concentration    (c - mMol/L)
8. Total volume                 (V - ml)
9. Expected Encapsulation efficiency (%) [optional] - provide only if defect value optimization is required.

'''
