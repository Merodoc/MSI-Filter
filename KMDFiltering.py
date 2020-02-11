# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from pyimzml.ImzMLParser import ImzMLParser
import numpy as np
from scipy import stats as st
import matplotlib.pyplot as plt
from pyimzml.ImzMLWriter import ImzMLWriter
import time
# Notes:
# New Plan - Matching to known lines with some level of variance
# completing regression on Defect vs Mass separately then mapping spectrum onto the regression line
# find values that match this line within certain tolerance and cut everything else 
# provide a large background dataset of known compounds and their defect ratios then classify them to these lines

#Need to filter mzlist by intensity


class DefectFilter:
    def __init__(self, filename):
        """ Initialize Filtering Framework from an imzml file """
        self.spectrum = ImzMLParser(filename)
        self.mzlist = []
        self.intensity_list = []
        self.filename = []
        self.KM2 = []
        self.KMD2 = []
        self.filtered_intens = []
        self.filtered_mass = []
        
        for idx, (x,y,z) in enumerate(self.spectrum.coordinates):
            self.mzs, self.intensities = self.spectrum.getspectrum(idx)
            self.mzlist.append(self.mzs)
            self.intensity_list.append(self.intensities)
        for i in range(len(self.intensity_list)):
            self.filterIntens(self.intensity_list[i], self.mzlist[i])
        
        self.filtered_mass = np.array(self.filtered_mass)
        self.filtered_intens = np.array(self.filtered_intens)
        
    def MSIFilter(self, coi, alpha):
        "Filter imzML file for complex of interest"
        if coi == "N-Glycan":
            self.glycanFilter()
            self.glycan_intens = []
            self.glycan_mzlist = self.filtered_mass
            for i in range(len(self.glycan_mzlist)):
                print("Iteration: "+str(i))
                if not np.all(self.glycan_mzlist[i]==0):
                    self.kendrickMass(self.glycan_mzlist[i])
                    self.glycan_intens.append(self.glycanProb(self.glycan_mzlist[i], self.KMD, alpha, self.filtered_intens[i]))
                else:
                    self.glycan_intens.append(self.glycan_mzlist[i])
            outname = 'glycanSpec'+str(np.random.randint(0,10**6))+'.imzML'
            mid = time.time()
            print("Iteration Time: "+str(mid-start))
            with ImzMLWriter(outname) as w:
                    for i in range(len(self.glycan_mzlist)):
                        w.addSpectrum(self.glycan_mzlist[i], self.glycan_intens[i], self.spectrum.coordinates[i])
            print("File Written to : " + outname)
        
    def glycanFilter(self, max_defect = 3):
        """create a line for the glycan filter based on ASMS 2019 poster """
        self.glycanMD = self.mzs*3.5*10**(-4) + 0.0039
        plt.plot(self.mzs, self.glycanMD)
        self.glycanDict = {}
        self.glycanSigma = 0.0173
        for i in range(len(self.mzs)):
            self.glycanDict[self.mzs[i]] = self.glycanMD[i]
            

    def glycanProb(self, KM, KMD, alpha, intensities, dist = 'Norm'):
        # Replace gylcanProb with a t-test or z-test from software
        """ Provide an intensity spectrum filtered for KMD values within alpha of known values """
        """ for single spectrum """
        glycanFilterInt = intensities.copy()
        for i in range(len(KM)):
            xbar = self.glycanDict[KM[i]]
            bestProb = 1
            for j in KMD[i]:
                if dist == 'Norm':
                    prob = st.norm.cdf(abs(xbar-j), loc = 0, scale = self.glycanSigma) - 0.5
                else:
                    break
                if prob < bestProb:
                    bestProb = prob
            if bestProb > alpha:
                glycanFilterInt[i] = 0
        return glycanFilterInt
                
    def kendrickMass(self, mzs, max_defect = 3):
        """ for single spectrum """
        # Start with KMs between 0 and 1:
        self.KMdict = {}
#        for i in range(max_defect + 1):
#            if i != 0:
#                self.KMdict[min(mzs)-i] = []
                
        for mz in mzs:
            self.KMdict[mz] = []

        for mz in mzs:
            defect, mass = np.modf(mz)
            for i in range(max_defect+1):
                    if mz-i in self.KMdict:
                        self.KMdict[mz - i].append(defect+i)
                    else:
                        self.KMdict[mz-i] = [defect+i]
        tbr = []
        for i in self.KMdict.keys():
            if i < min(mzs):
                tbr.append(i)
        
        for i in tbr:
            del self.KMdict[i]
            
        self.KM = list(self.KMdict.keys())
        self.KMD = list(self.KMdict.values())
            
    def kendrickMassList(self, mzs):
        """ for single spectrum """
        KM = mzs*14/14.01565
        self.KM2.append(KM)
        KMD = np.floor(KM) - KM
        self.KMD2.append(KMD)
    
    def KMDplot(self):
        axes = plt.axes()
        axes.set_ylim([-1,0])
        for i in range(len(self.filtered_mass)):
            plt.scatter(self.filtered_mass[i], self.KMD2[i])
        plt.show()
            
    def filterIntens(self, intens_list, mzlist,thresh = 0):
        print("iteration")
        intens = []
        mz = []
        if  np.all(intens_list <= thresh):
            np.append(self.filtered_intens, intens_list)
            np.append(self.filtered_mass, mzlist)
        else:
            for i in range(len(intens_list)):
                if intens_list[i] > thresh:
                    intens.append(intens_list[i])
                    mz.append(mzlist[i])
                else:
                    intens.append(0)
                    mz.append(0)
        np.append(self.filtered_intens, intens)
        np.append(self.filtered_mass, mz)
    
    def kendrickFilter(self, thresh,intens_list, mzlist):
        """ Takes full spectrum lists not single spectrum """
        for i in range(len(intens_list)):
            self.filterIntens(thresh, intens_list[i], mzlist[i])
        for i in self.filtered_mass:
            self.kendrickMassList(i)
                
                
                
        
        
            
        
                
            


start = time.time()    
test = DefectFilter('2019-110.imzML')
end = time.time()
print(end-start)
#test.glycanFilter()

start = time.time()
#test.KMDplot()
test.MSIFilter("N-Glycan",0.05)
end = time.time()
print(end-start)

"""
### 1. imzML data import
from pyimzml.ImzMLParser import ImzMLParser
import numpy as np

p = ImzMLParser('anenome.imzML')
mzlist = []
intensity_list = []
for idx, (x,y,z) in enumerate(p.coordinates):
    mzs, intensities = p.getspectrum(idx)
    mzlist.append(mzs)
    intensity_list.append(intensities)
### 2. Extratc mass spectrum from (x,y) position

### 3. Calculate KMD for each m/z

### 4. Kendrick Algorithm
Mass_list = []
Defect_list = []
MDR_list = []
for spectrum in mzlist:
    MDR = []
    Masses = []
    Defects = []
    masses = []
    defects = []
    for mzs in spectrum:
        Defect, Mass = np.modf(mzs)
        Masses.append(Mass)
        Defects.append(Defect)
        
        masses.append(Mass)
        defects.append(Defect)
        
        MDR.append(Defect/Mass)
    MDR_list.append(MDR)
    Mass_list.append(Masses)
    Defect_list.append(Defects)
    
import matplotlib.pyplot as plt
from sklearn import linear_model, metrics

X = masses
y = defects

from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X,y,test_size = 0.4, random_state = 1)

reg = linear_model.LinearRegression()

reg.fit(X_train, y_train)

print('Coefficients: \n', reg.coef_)

print('Variance score: {}'.format(reg.score(X_test, y_test)))

plt.plot(X_test, y_test, 'ro')
"""
""" Where Nominal - Mkr is the nominal mass and Mkr is the exact mass of the kendrick reference used for the amu definition """


### 5. i) if more than one chemically related compound 
    ### Mass Difference Clustering Algorithm
""" Use a standard clustering algorithm? """
### 5. ii) if "Targeted method?"
    ### Mass spectrum filtered
        ### generate image
### Back to step 3 with different KMD and m/z ranges 
"""split MDR into multiple spectra """
"""size = len(MDR)
idx_list = [idx + 1 for idx, val in enumerate(MDR) if val == 0]
res = [MDR[i:j] for i,j in zip([0]+idx_list,idx_list + ([size] if idx_list[-1] != size else []))]
"""

#
#from pyimzml.ImzMLWriter import ImzMLWriter
#with ImzMLWriter('output.imzML') as w:
#        for i in range(len(MDR_list)):
#            w.addSpectrum(MDR_list[i]*10**4, intensity_list[i], p.coordinates[i])
#            
#p = ImzMLParser('output.imzML')
#for idx, (x,y,z) in enumerate(p.coordinates):
#    mzs, intensities = p.getspectrum(idx)
#    print(len(mzs))
