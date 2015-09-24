__author__ = 'Zhiyang_Lin'
import sys
import numpy as np
import re
import os
import matplotlib.pyplot as plt
from Imagetool import *

def main():
    ####################################################
    # Read the input file
    ####################################################
    datafilename = 'output_it00001500'      #Modify the name to change the Inputfile 
    gbfilename = 'voronoi_7575_15.dat'      #Modify the name to change the Grain boundary file
    data = open(datafilename,'r')       #Open the input data file
    points = data.readlines()           #Read the input data file line by line
    gbdata = open(gbfilename,'r')       #Open the grain boundary file
    gbpoints = gbdata.readlines()       #Read the grain boundary file 
    lastpoint = points[-1] # Read the last line find out the maximum x, y and z
    print(lastpoint)
    lp = lastpoint.split()
    xmax = int(lp[0])
    print(xmax)
    ymax = int(lp[1])
    zmax = int(lp[2])
    print('zmax is' + str(zmax))
    gbzset = []             #Get the right z value for gb
    for z in range(zmax+1):
        gbzset.append(z*4)
    image = []
    c = 0
    ####################################################
    # Create image objects for all layers of data
    ####################################################
    for i in range(zmax+1):
        print('create image' + str(c))
        image.append(Image(xmax,ymax))
        c += 1
    ####################################################
    # For each image object collect pixel objects 
    ####################################################
    for p in points:
        data = p.split()
        x = int(data[0])
        y = int(data[1])
        z = int(data[2])
        pf = float(data[9])
        if z in range(zmax+1):
        	image[z].pixelcollect(x,y,z,pf)
    ####################################################
    # Transform the pixel sets to matrices
    ####################################################
    c = 0
    for i in image:
        print('form matrix of image' + str(c))
        i.formmatrix()
        c += 1
    ####################################################
    # Mark all grain boundary points
    ####################################################
    for gbp in gbpoints:
        gbdata = gbp.split()
        x = int(gbdata[0])
        y = int(gbdata[1])
        z = int(gbdata[2])
        gb = int(gbdata[3])
        if(gb == 1 and z in gbzset):
            image[z/4].markgb(x,y)
    c = 0
    ####################################################
    # Process every image
    ####################################################
    delgb = 0 # change it to 1 if you want to get rid off grain boudary
    clusters = {}
    for i in image:
        if delgb == 1:
            i.blendgb()
        else:
            i.enlargegb()
        i.getdlset()
        i.holefill()
        i.hrd()
        
        for p in i.edgesp:
        	if i.edgedetect(p) == 5 or i.edgedetect(p) == 10:
        		print('pass')
        		pass
        	else:
        		i.edges.append(i.edgecollect(p)) #find all edge points for every cluster
        xset = []
        yset = []
        ####################################################
	    # output a image file and a text file for all the
	    # edge points
	    ####################################################
        opfilename = 'plane' + str(c) + 'edge.txt'
        edgeout = open(opfilename,'w')
        for n in i.edges:
        	edgeop = 'cluster' + str(n[0].label) + 'edge:'
        	for p in n:
        		edgeop = edgeop + '(' + str(p.x) + ',' + str(p.y) + ') '
        		xset.append(p.x)
        		yset.append(p.y)
        	edgeout.write(edgeop)
        	edgeout.write('\n')
    	plt.figure()
        plt.figure(figsize=(10,10))
        plt.scatter(xset,yset)
        plt.ylim(0, ymax)
        plt.xlim(0, xmax)
        figname = 'plane'+ str(c) + 'edges' +'.png'
        plt.savefig(figname)
        plt.close()

        ####################################################
	    # output a image file and text file for all the clusters 
	    ####################################################

        print('plot image'+ str(c))
        xset = []
        yset = []
        zset = []
        for x in range(xmax + 1):
            for y in range(ymax + 1):
                if i.imagematrix[x][y].sets != 0:
                    xset.append(x)
                    yset.append(y)
                    zset.append(i.imagematrix[x][y].label)
        plt.figure()
        plt.figure(figsize=(10,10))
        plt.scatter(xset,yset,c=zset)
        plt.ylim(0, ymax)
        plt.xlim(0, xmax)
        figname = 'plane'+ str(c) + 'clusters' +'gb.png'
        plt.savefig(figname)
        plt.close()
        c += 1
        
        for l in range(i.maxlabel+1):
            if l != 0:
                ttpoints = 0
                for x in range(xmax+1):
                    for y in range(ymax+1):
                        if i.imagematrix[x][y].label == l:
                            ttpoints += 1
                            cset = i.imagematrix[x][y].sets
            	if cset == 1:
            		ttpoints = float(ttpoints)
                if ttpoints in clusters:
                    clusters[ttpoints] += 1
                else:
                    clusters[ttpoints] = 1

    ckeys = clusters.keys()
    ckeys.sort()
    fout = open('output.txt','w')
    print('output data')
    for i in ckeys:
        info = []
        info.append(str(i))
        info.append(str(clusters[i]))
        fout.write('{0[0]:<15}{0[1]:<15}'.format(info))
        fout.write('\n')

    num = []
    for i in ckeys:
        num.append(clusters[i])
    plt.figure()
    plt.figure(figsize=(10,10))
    plt.ylabel('number of dislocation with same area')
    plt.xlabel('area of dislocation')
    plt.loglog(ckeys,num,'bo')
    plt.savefig('logspace.png')
    plt.close()
    

if __name__ == '__main__':
    main()
