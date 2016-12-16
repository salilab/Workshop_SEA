#!/bin/python

import os,sys,string,numpy
import scipy.stats,math

weightsFile=sys.argv[1] # format is one row per cluster and one column per run
wtsArray=numpy.transpose(numpy.loadtxt(weightsFile))

#check 1: make sure order of magnitude of elements is in range 10's i.e 10-100's (i.e. 2-3) 
print "Chi-square expected counts should be in the order of magnitude of 10's. Otherwise chi-square estimates are artificially blown up."
print "Maximum element is",wtsArray.max()

if wtsArray.max()>200.0:
        scalingFactor=math.pow(10,int(math.floor(math.log10(wtsArray.max())))-1)
        wtsArray=wtsArray/scalingFactor

print wtsArray

# check 2: make sure each element is greater than 10 
print "Minimum element is ",wtsArray.min()
print "All expected counts should be 10 or greater. Cochran (1952, 1954)"

if wtsArray.min()<10.0:
	scalingFactor=10.0/wtsArray.min()
	wtsArray=wtsArray*scalingFactor

print wtsArray

[chisquare,pvalue,dof,expected]=scipy.stats.chi2_contingency(wtsArray)

#print "When one set of marginal totals--the rows, say--is fixed by the sampling scheme, the hypothesis of no association is called homogeneity of proportions. It says the proportion of individuals in a particular column the same for all rows. If this does not hold for your data, the p-value will be significant i.e. less than 0.05"
#print "Chi-square value",chisquare
#print "Expected dof",dof 
print "If p value below is significant, i.e less than 0.05, then the two runs are not from the same distribution (and sampling hasn't converged). Else the two runs are from the same distribution (and you did a good job with sampling!)"
print "P-value",pvalue
#print "Expected distribution",expected


