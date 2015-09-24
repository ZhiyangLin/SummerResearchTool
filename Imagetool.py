__author__ = 'Zhiyang_Lin'
import copy
import math
import operator

####################################################
# Pixel class: a pixel object contains all the info
# of a data point.
####################################################

class Pixel:
    def __init__(self, x=0, y=0, z=0, pf=0.0):
        self.x = x
        self.y = y
        self.z = z
        self.pf = pf
        self.sets = 0 # sets is the value of the initial cutoff division
        self.selected = 0 # indicates whether or not a point has been selected or not
        self.label = 0 # cluster label
        self.gb = 0 # indicates whether or not a point is a grain boundary

    def __str__(self):
        return "({},{}, pf = {}, label = {}/n)".format(self.x, self.y, self.pf, self.label)
####################################################
# Image class: a image class contains all the pixels
####################################################
class Image:
    def __init__(self,xmax,ymax):
        self.xmax = xmax
        self.ymax = ymax
        self.pixels = []
        self.threshold = 0.03
        self.gbpoints = []
        self.clusters = {}
        self.edges = []

    def __str__(self):
        return "(xmax = {}, ymax = {})".format(self.xmax, self.ymax)
    ####################################################
    # pixelcollect: collect all the pixels into a list 
    # belongs to this image according to their z value
    ####################################################
    def pixelcollect(self,x,y,z,pf):
            self.pixels.append(Pixel(x,y,z,pf))
    ######################################################
    # formmatrix: form the pixel list into a x by y matrix
    ######################################################
    def formmatrix(self):
        self.imagematrix = [] #create a empty list
        temppixel = copy.deepcopy(self.pixels) # make a deep copy of the pixel list 
        while temppixel != []: # while the copy of pixel list is not empty
            self.imagematrix.append(temppixel[:self.xmax+1]) # get the first x pixels into the matrix list
            temppixel = temppixel[self.xmax+1:] #delete the first x pixels in the copy

    ######################################################
    # markgb: mark the grain boundary points by change the
    # gb value to 1
    ######################################################
    def markgb(self,x,y):
        self.imagematrix[x][y].gb = 1
        self.gbpoints.append([x,y]) # record all the grain boundary points in a gbpoints list

    ######################################################
    # blendgb: get rid of the gb points
    ######################################################
    def blendgb(self):
        xmax = self.xmax
        ymax = self.ymax
        t = self.threshold
        while self.gbpoints: # while the gbpoints list is not empty
            x = self.gbpoints[0][0] # read the x value of the first gb point in the list
            y = self.gbpoints[0][1] # read the y value of the first gb point in the list
            del self.gbpoints[0] # delete the first gbpoint in the list
            tempx = x
            tempy = y
            if x == 0 and y == 0: #left bottom
                c = 0
                while tempx < xmax and c < 2:  # find 2 none-gb point along the x axis direction
                    tempx = tempx + 1
                    if self.imagematrix[tempx][y].gb == 0:
                        c = c + 1
                c = 0
                while tempy < ymax and c < 2:  # find 2 none-gb point along the y axis direction
                    tempy = tempy + 1
                    if self.imagematrix[x][tempy].gb == 0:
                        c = c + 1
                if tempx < tempy:   # if the gb points go along the y direction, get rid of the grain boundary alone 
                    i = 0           # the x direction
                    while x + i < tempx:
                        self.imagematrix[x + i][y].pf = self.imagematrix[tempx][y].pf
                        i = i + 1
                else:               # if the gb points go along the x direction, get rid of the grain boundary alone 
                    i = 0           # the y direction 
                    while y + i < tempy:
                        self.imagematrix[x][y + i].pf = self.imagematrix[x][tempy].pf
                        i = i + 1

            if (0 < x < xmax) and y == 0:
                c = 0
                while tempx < xmax and c < 2:
                    tempx = tempx + 1
                    if self.imagematrix[tempx][y].gb == 0:
                        c = c + 1
                if c > 0 and abs(self.imagematrix[tempx][y].pf - self.imagematrix[x-1][y].pf) < t:
                    i = 0
                    while 0 < x - i and x + i < tempx:
                        self.imagematrix[x + i][y].pf = self.imagematrix[tempx][y].pf
                        self.imagematrix[x - i][y].pf = self.imagematrix[tempx][y].pf
                        i = i + 1
                else:
                    c = 0
                    while tempy < ymax and c < 2:
                        tempy = tempy + 1
                        if self.imagematrix[x][tempy].gb == 0:
                            c = c + 1
                    i = 0
                    while y + i < tempy:
                        self.imagematrix[x][y + i].pf = self.imagematrix[x][tempy].pf
                        i = i + 1

            if x == xmax and y == 0:
                self.imagematrix[x][y].pf = self.imagematrix[x-2][y].pf
                self.imagematrix[x-1][y].pf = self.imagematrix[x-2][y].pf


            if x == 0 and (0 < y < ymax):
                c = 0
                while tempy < ymax and c < 2:
                    tempy = tempy + 1
                    if self.imagematrix[x][tempy].gb == 0:
                        c = c + 1
                if c > 0 and abs(self.imagematrix[x][tempy].pf - self.imagematrix[x][y-1].pf) < t:
                    i = 0
                    while 0 < y - i and y + i < tempy:
                        self.imagematrix[x][y + i].pf = self.imagematrix[x][tempy].pf
                        self.imagematrix[x][y - i].pf = self.imagematrix[x][tempy].pf
                        i = i + 1
                else:
                    c = 0
                    while tempx < xmax and c < 2:
                        tempx = tempx + 1
                        if self.imagematrix[x][tempy].gb == 0:
                            c = c + 1
                    i = 0
                    while x + i < tempy:
                        self.imagematrix[x + 1][y].pf = self.imagematrix[tempx][y].pf
                        i = i + 1

            if  0 < x < xmax and 0 < y < ymax:
                c = 0
                while tempx < xmax and c < 2:
                    tempx = tempx + 1
                    if self.imagematrix[tempx][y].gb == 0:
                        c = c + 1
                        if c == 1:
                            tempx1 = tempx
                if c > 0 and abs(self.imagematrix[tempx][y].pf - self.imagematrix[x-1][y].pf) < t:
                    i = 0
                    while 0 < x - i and x + i < tempx:
                        self.imagematrix[x + i][y].pf = self.imagematrix[tempx][y].pf
                        self.imagematrix[x - i][y].pf = self.imagematrix[tempx][y].pf
                        i = i + 1
                else:
                    c = 0
                    while tempy < ymax and c < 2:
                        tempy = tempy + 1
                        if self.imagematrix[x][tempy].gb == 0:
                            c = c + 1
                            if c == 1:
                                tempy1 = tempy
                    if c > 0 and abs(self.imagematrix[x][tempy].pf - self.imagematrix[x][y-2].pf) < t:
                        i = 0
                        while 0 < y - i and y + i < tempy:
                            self.imagematrix[x][y + i].pf = self.imagematrix[x][tempy].pf
                            self.imagematrix[x][y - i].pf = self.imagematrix[x][tempy].pf
                            i = i + 1

                    else:
                        if tempx - x < tempy - y:
                            self.imagematrix[x][y].pf = self.imagematrix[tempx1][y].pf
                        else:
                            self.imagematrix[x][y].pf = self.imagematrix[x][tempy1].pf


            if x == xmax and 0 < y < ymax:
                c = 0
                while tempy < ymax and c < 2:
                    tempy = tempy + 1
                    if self.imagematrix[x][tempy].gb == 0:
                        c = c + 1
                if c > 0 and abs(self.imagematrix[x][tempy].pf - self.imagematrix[x][y-1].pf) < t:
                    i = 0
                    while 0 < y - i and y + i < tempy:
                        self.imagematrix[x][y + i].pf = self.imagematrix[x][tempy].pf
                        self.imagematrix[x][y - i].pf = self.imagematrix[x][tempy].pf
                        i = i + 1
                else:
                    self.imagematrix[x][y].pf = self.imagematrix[x-2][y].pf
                    self.imagematrix[x-1][y].pf = self.imagematrix[x-2][y].pf

            if x == 0 and y == ymax:
                self.imagematrix[x][y].pf = self.imagematrix[x][y-2].pf
                self.imagematrix[x][y-1].pf = self.imagematrix[x][y-2].pf

            if (0 < x < xmax) and y == ymax:
                c = 0
                while tempx < xmax and c < 2:
                    tempx = tempx + 1
                    if self.imagematrix[tempx][y].gb == 0:
                        c = c + 1
                if c > 0 and abs(self.imagematrix[tempx][y].pf - self.imagematrix[x-1][y].pf) < t:
                    i = 0
                    while 0 < x - i and x + i < tempx:
                        self.imagematrix[x + i][y].pf = self.imagematrix[tempx][y].pf
                        self.imagematrix[x - i][y].pf = self.imagematrix[tempx][y].pf
                        i = i + 1
                else:
                    self.imagematrix[x][y].pf = self.imagematrix[x][y-2].pf
                    self.imagematrix[x][y-1].pf = self.imagematrix[x][y-2].pf

            if x == xmax and y == ymax:
                self.imagematrix[x][y].pf = self.imagematrix[x][y-2].pf
                self.imagematrix[x][y-1].pf = self.imagematrix[x][y-2].pf

	######################################################
	# enlargegb: enlarge grain boundary separate clusters
	######################################################
    def enlargegb(self):
        xmax = self.xmax
        ymax = self.ymax
        t = self.threshold
        while self.gbpoints: # while the gbpoints list is not empty
            x = self.gbpoints[0][0] # read the x value of the first gb point in the list
            y = self.gbpoints[0][1] # read the y value of the first gb point in the list
            del self.gbpoints[0] # delete the first gbpoint in the list
            tempx = x
            tempy = y
            if x == 0 and y == 0: # left bottom
                c = 0
                while tempx < xmax and c < 2: # find 2 none-gb point along the x axis direction
                    tempx = tempx + 1
                    if self.imagematrix[tempx][y].gb == 0:
                        c = c + 1
                c = 0
                while tempy < ymax and c < 2: # find 2 none-gb point along the y axis direction
                    tempy = tempy + 1
                    if self.imagematrix[x][tempy].gb == 0:
                        c = c + 1
                if tempx < tempy: # if the gb points go along the y direction, enlarge the grain boundary
                    i = 0         # along x direction
                    while x + i < tempx:
                        self.imagematrix[x + i][y].gb = 1
                        i = i + 1
                else:               # if the gb points go along the y direction, enlarge the grain boundary
                    i = 0           # along y direction
                    while y + i < tempy:
                        self.imagematrix[x][y + i].gb = 1
                        i = i + 1
            if 0 < x < xmax and y == 0:
                c = 0
                while tempx < xmax and c < 2:
                    tempx = tempx + 1
                    if self.imagematrix[tempx][y].gb == 0:
                        c = c + 1
                if c > 0: 
                    i = 0                                                                           
                    while 0 < x - i and x + i < tempx:                                             
                        self.imagematrix[x + i][y].gb = 1                   
                        self.imagematrix[x - i][y].gb = 1
                        i = i + 1
                else:
                    c = 0
                    while tempy < ymax and c < 2:
                        tempy = tempy + 1
                        if self.imagematrix[x][tempy].gb == 0:
                            c = c + 1
                    i = 0
                    while y + i < tempy:
                        self.imagematrix[x][y + i].gb = 1
                        i = i + 1
            if x == xmax and y == 0:
                self.imagematrix[x-1][y].gb = 1

            if x == 0 and (0 < y < ymax):
                c = 0
                while tempy < ymax and c < 2:
                    tempy = tempy + 1
                    if self.imagematrix[x][tempy].gb == 0:
                        c = c + 1
                if c > 0:
                    i = 0
                    while 0 < y - i and y + i < tempy:
                        self.imagematrix[x][y + i].gb = 1
                        self.imagematrix[x][y - i].gb = 1
                        i = i + 1
                else:
                    c = 0
                    while tempx < xmax and c < 2:
                        tempx = tempx + 1
                        if self.imagematrix[x][tempy].gb == 0:
                            c = c + 1
                    i = 0
                    while x + i < tempy:
                        self.imagematrix[x + 1][y].gb = 1
                        i = i + 1
            if  0 < x < xmax and 0 < y < ymax:
                c = 0
                while tempx < xmax and c < 2:
                    tempx = tempx + 1
                    if self.imagematrix[tempx][y].gb == 0:
                        c = c + 1
                if c > 0:
                    i = 0
                    while 0 < x - i and x + i < tempx:
                        self.imagematrix[x + i][y].gb = 1
                        self.imagematrix[x - i][y].gb = 1
                        i = i + 1
                else:
                    c = 0
                    while tempy < ymax and c < 2:
                        tempy = tempy + 1
                        if self.imagematrix[x][tempy].gb == 0:
                            c = c + 1
                    if c > 0:
                        i = 0
                        while 0 < y - i and y + i < tempy:
                            self.imagematrix[x][y + i].gb = 1
                            self.imagematrix[x][y - i].gb = 1
                            i = i + 1
                    else:
                        if tempx - x < tempy - y:
                            self.imagematrix[x][y].gb = 1
                        else:
                            self.imagematrix[x][y].gb = 1
            if x == xmax and 0 < y < ymax:
                c = 0
                while tempy < ymax and c < 2:
                    tempy = tempy + 1
                    if self.imagematrix[x][tempy].gb == 0:
                        c = c + 1
                if c > 0:
                    i = 0
                    while 0 < y - i and y + i < tempy:
                        self.imagematrix[x][y + i].gb = 1
                        self.imagematrix[x][y - i].gb = 1
                        i = i + 1
                else:
                    self.imagematrix[x-1][y].gb = 1
            if x == 0 and y == ymax:
                self.imagematrix[x][y-1].gb = 1
            if (0 < x < xmax) and y == ymax:
                c = 0
                while tempx < xmax and c < 2:
                    tempx = tempx + 1
                    if self.imagematrix[tempx][y].gb == 0:
                        c = c + 1
                if c > 0:
                    i = 0
                    while 0 < x - i and x + i < tempx:
                        self.imagematrix[x + i][y].gb = 1
                        self.imagematrix[x - i][y].gb = 1
                        i = i + 1
                else:
                    self.imagematrix[x][y-1].gb = 1
            if x == xmax and y == ymax:
                self.imagematrix[x][y].gb = 1
                self.imagematrix[x][y-1].gb = 1
######################################################
# getdlset: separate all the data point into different
# sets before clustering them
######################################################
    def getdlset(self): 
        self.dlsets = []
        self.dlsets.append([0, 0.4])
        self.dlsets.append([0.4, 0.9])
        self.dlsets.append([0.9, 1.1]) # append as many sets as you want
        for x in range(self.xmax+1):
            for y in range(self.ymax+1):
                for i,s in enumerate(self.dlsets):
                    if s[0] <= self.imagematrix[x][y].pf < s[1]:
                        self.imagematrix[x][y].sets = i
                if self.imagematrix[x][y].gb == 1: # if the point is a grain boundary point, classify it as sets 0
                        self.imagematrix[x][y].sets = 0
######################################################
# hrd: cluster all the data points
######################################################
    def hrd(self):
        xmax = self.xmax
        ymax = self.ymax
        maxlabel = 0
        edgesp = [] # A list to record a begining edge point of all clusters 
        for x in range(xmax + 1):
            for y in range(ymax + 1):
                points = [] # a list to add points that need to be processed
                if self.imagematrix[x][y].sets == 0: # if this sets belong to sets 0, then label it as cluster 0.
                    self.imagematrix[x][y].label = 0
                    self.imagematrix[x][y].selected = 1 # mark this point as selected after it has been processed
                if self.imagematrix[x][y].selected == 1:# if this point selected then skip it
                    pass
                else:
                    maxlabel += 1 # While maxlabel increases, it means it is processing a new cluster
                    self.imagematrix[x][y].label = maxlabel
                    edgesp.append(self.imagematrix[x][y]) #append this to the edge list
                    points.append([x,y])                  #append the [x,y] to the point list
                    while points: # while this points list is not empty.
                        tempx = points[0][0] # get the x value of the current point
                        tempy = points[0][1] # get the y value of the current point
                        del points[0] # delete the first point in the list
                        sets = self.imagematrix[tempx][tempy].sets
                        ######################################################
                        # check left point, right point, top and bottom point.
                        # to see if it has the same sets value as current point
                        #. If yes, append it to the points list.
                        ######################################################
                        if tempx == 0 and tempy == 0:
                            if self.imagematrix[xmax][tempy].selected == 0 and self.imagematrix[xmax][tempy].sets == sets:
                                self.imagematrix[xmax][tempy].selected = 1
                                self.imagematrix[xmax][tempy].label = maxlabel
                                points.append([xmax,tempy])
                            if self.imagematrix[tempx + 1][tempy].selected == 0 and self.imagematrix[tempx + 1][tempy].sets == sets:
                                self.imagematrix[tempx + 1][tempy].selected = 1
                                self.imagematrix[tempx + 1][tempy].label = maxlabel
                                points.append([tempx + 1,tempy])
                            if self.imagematrix[tempx][tempy + 1].selected == 0 and self.imagematrix[tempx][tempy + 1].sets == sets:
                                self.imagematrix[tempx][tempy + 1].selected = 1
                                self.imagematrix[tempx][tempy + 1].label = maxlabel
                                points.append([tempx,tempy + 1])
                            if self.imagematrix[tempx][ymax].selected == 0 and self.imagematrix[tempx][ymax].sets == sets:
                                self.imagematrix[tempx][ymax].selected = 1
                                self.imagematrix[tempx][ymax].label = maxlabel
                                points.append([tempx, ymax])
                        elif 0 < tempx < xmax and tempy == 0:
                            if self.imagematrix[tempx - 1][tempy].selected == 0 and self.imagematrix[tempx - 1][tempy].sets == sets:
                                self.imagematrix[tempx - 1][tempy].selected = 1
                                self.imagematrix[tempx - 1][tempy].label = maxlabel
                                points.append([tempx - 1, tempy])
                            if self.imagematrix[tempx + 1][tempy].selected == 0 and self.imagematrix[tempx + 1][tempy].sets == sets:
                                self.imagematrix[tempx + 1][tempy].selected = 1
                                self.imagematrix[tempx + 1][tempy].label = maxlabel
                                points.append([tempx + 1,tempy])
                            if self.imagematrix[tempx][tempy + 1].selected == 0 and self.imagematrix[tempx][tempy + 1].sets == sets:
                                self.imagematrix[tempx][tempy + 1].selected = 1
                                self.imagematrix[tempx][tempy + 1].label = maxlabel
                                points.append([tempx,tempy + 1])
                            if self.imagematrix[tempx][ymax].selected == 0 and self.imagematrix[tempx][ymax].sets == sets:
                                self.imagematrix[tempx][ymax].selected = 1
                                self.imagematrix[tempx][ymax].label = maxlabel
                                points.append([tempx,ymax])
                        elif tempx == xmax and tempy == 0:
                            if self.imagematrix[tempx - 1][tempy].selected == 0 and self.imagematrix[tempx - 1][tempy].sets == sets:
                                self.imagematrix[tempx - 1][tempy].selected = 1
                                self.imagematrix[tempx - 1][tempy].label = maxlabel
                                points.append([tempx - 1, tempy])
                            if self.imagematrix[0][tempy].selected == 0 and self.imagematrix[0][tempy].sets == sets:
                                self.imagematrix[0][tempy].selected = 1
                                self.imagematrix[0][tempy].label = maxlabel
                                points.append([0,tempy])
                            if self.imagematrix[tempx][tempy + 1].selected == 0 and self.imagematrix[tempx][tempy + 1].sets == sets:
                                self.imagematrix[tempx][tempy + 1].selected = 1
                                self.imagematrix[tempx][tempy + 1].label = maxlabel
                                points.append([tempx,tempy + 1])
                            if self.imagematrix[tempx][ymax].selected == 0 and self.imagematrix[tempx][ymax].sets == sets:
                                self.imagematrix[tempx][ymax].selected = 1
                                self.imagematrix[tempx][ymax].label = maxlabel
                                points.append([tempx,ymax])
                        elif tempx == 0 and 0 < tempy < ymax:
                            if self.imagematrix[xmax][tempy].selected == 0 and self.imagematrix[xmax][tempy].sets == sets:
                                self.imagematrix[xmax][tempy].selected = 1
                                self.imagematrix[xmax][tempy].label = maxlabel
                                points.append([xmax, tempy])
                            if self.imagematrix[tempx + 1][tempy].selected == 0 and self.imagematrix[tempx + 1][tempy].sets == sets:
                                self.imagematrix[tempx + 1][tempy].selected = 1
                                self.imagematrix[tempx + 1][tempy].label = maxlabel
                                points.append([tempx + 1,tempy])
                            if self.imagematrix[tempx][tempy + 1].selected == 0 and self.imagematrix[tempx][tempy + 1].sets == sets:
                                self.imagematrix[tempx][tempy + 1].selected = 1
                                self.imagematrix[tempx][tempy + 1].label = maxlabel
                                points.append([tempx,tempy + 1])
                            if self.imagematrix[tempx][tempy - 1].selected == 0 and self.imagematrix[tempx][tempy - 1].sets == sets:
                                self.imagematrix[tempx][tempy - 1].selected = 1
                                self.imagematrix[tempx][tempy - 1].label = maxlabel
                                points.append([tempx,tempy - 1])
                        elif 0 < tempx < xmax and 0 < tempy < ymax:
                            if self.imagematrix[tempx - 1][tempy].selected == 0 and self.imagematrix[tempx - 1][tempy].sets == sets:
                                self.imagematrix[tempx - 1][tempy].selected = 1
                                self.imagematrix[tempx - 1][tempy].label = maxlabel
                                points.append([tempx - 1, tempy])
                            if self.imagematrix[tempx + 1][tempy].selected == 0 and self.imagematrix[tempx + 1][tempy].sets == sets:
                                self.imagematrix[tempx + 1][tempy].selected = 1
                                self.imagematrix[tempx + 1][tempy].label = maxlabel
                                points.append([tempx + 1,tempy])
                            if self.imagematrix[tempx][tempy + 1].selected == 0 and self.imagematrix[tempx][tempy + 1].sets == sets:
                                self.imagematrix[tempx][tempy + 1].selected = 1
                                self.imagematrix[tempx][tempy + 1].label = maxlabel
                                points.append([tempx,tempy + 1])
                            if self.imagematrix[tempx][tempy - 1].selected == 0 and self.imagematrix[tempx][tempy - 1].sets == sets:
                                self.imagematrix[tempx][tempy - 1].selected = 1
                                self.imagematrix[tempx][tempy - 1].label = maxlabel
                                points.append([tempx,tempy - 1])
                        elif tempx == xmax and 0 < tempy < ymax:
                            if self.imagematrix[tempx - 1][tempy].selected == 0 and self.imagematrix[tempx - 1][tempy].sets == sets:
                                self.imagematrix[tempx - 1][tempy].selected = 1
                                self.imagematrix[tempx - 1][tempy].label = maxlabel
                                points.append([tempx - 1, tempy])
                            if self.imagematrix[0][tempy].selected == 0 and self.imagematrix[0][tempy].sets == sets:
                                self.imagematrix[0][tempy].selected = 1
                                self.imagematrix[0][tempy].label = maxlabel
                                points.append([0,tempy])
                            if self.imagematrix[tempx][tempy + 1].selected == 0 and self.imagematrix[tempx][tempy + 1].sets == sets:
                                self.imagematrix[tempx][tempy + 1].selected = 1
                                self.imagematrix[tempx][tempy + 1].label = maxlabel
                                points.append([tempx,tempy + 1])
                            if self.imagematrix[tempx][tempy - 1].selected == 0 and self.imagematrix[tempx][tempy - 1].sets == sets:
                                self.imagematrix[tempx][tempy - 1].selected = 1
                                self.imagematrix[tempx][tempy - 1].label = maxlabel
                                points.append([tempx,tempy - 1])
                        elif tempx == 0 and tempy == ymax:
                           if self.imagematrix[xmax][tempy].selected == 0 and self.imagematrix[xmax][tempy].sets == sets:
                               self.imagematrix[xmax][tempy].selected = 1
                               self.imagematrix[xmax][tempy].label = maxlabel
                               points.append([xmax, tempy])
                           if self.imagematrix[tempx + 1][tempy].selected == 0 and self.imagematrix[tempx + 1][tempy].sets == sets:
                               self.imagematrix[tempx + 1][tempy].selected = 1
                               self.imagematrix[tempx + 1][tempy].label = maxlabel
                               points.append([tempx + 1,tempy])
                           if self.imagematrix[tempx][0].selected == 0 and self.imagematrix[tempx][0].sets == sets:
                               self.imagematrix[tempx][0].selected = 1
                               self.imagematrix[tempx][0].label = maxlabel
                               points.append([tempx,0])
                           if self.imagematrix[tempx][tempy - 1].selected == 0 and self.imagematrix[tempx][tempy - 1].sets == sets:
                               self.imagematrix[tempx][tempy - 1].selected = 1
                               self.imagematrix[tempx][tempy - 1].label = maxlabel
                               points.append([tempx,tempy - 1])
                        elif 0 < tempx < self.xmax and tempy == self.ymax:
                            if self.imagematrix[tempx - 1][tempy].selected == 0 and self.imagematrix[tempx - 1][tempy].sets == sets:
                                self.imagematrix[tempx - 1][tempy].selected = 1
                                self.imagematrix[tempx - 1][tempy].label = maxlabel
                                points.append([tempx - 1, tempy])
                            if self.imagematrix[tempx + 1][tempy].selected == 0 and self.imagematrix[tempx + 1][tempy].sets == sets:
                                self.imagematrix[tempx + 1][tempy].selected = 1
                                self.imagematrix[tempx + 1][tempy].label = maxlabel
                                points.append([tempx + 1,tempy])
                            if self.imagematrix[tempx][0].selected == 0 and self.imagematrix[tempx][0].sets == sets:
                                self.imagematrix[tempx][0].selected = 1
                                self.imagematrix[tempx][0].label = maxlabel
                                points.append([tempx,0])
                            if self.imagematrix[tempx][tempy - 1].selected == 0 and self.imagematrix[tempx][tempy - 1].sets == sets:
                                self.imagematrix[tempx][tempy - 1].selected = 1
                                self.imagematrix[tempx][tempy - 1].label = maxlabel
                                points.append([tempx,tempy - 1])
                        elif tempx == xmax and tempy == ymax:
                            if self.imagematrix[tempx - 1][tempy].selected == 0 and self.imagematrix[tempx - 1][tempy].sets == sets:
                                self.imagematrix[tempx - 1][tempy].selected = 1
                                self.imagematrix[tempx - 1][tempy].label = maxlabel
                                points.append([tempx - 1, tempy])
                            if self.imagematrix[0][tempy].selected == 0 and self.imagematrix[0][tempy].sets == sets:
                                self.imagematrix[0][tempy].selected = 1
                                self.imagematrix[0][tempy].label = maxlabel
                                points.append([0,tempy])
                            if self.imagematrix[tempx][0].selected == 0 and self.imagematrix[tempx][0].sets == sets:
                                self.imagematrix[tempx][0].selected = 1
                                self.imagematrix[tempx][0].label = maxlabel
                                points.append([tempx,0])
                            if self.imagematrix[tempx][tempy - 1].selected == 0 and self.imagematrix[tempx][tempy - 1].sets == sets:
                                self.imagematrix[tempx][tempy - 1].selected = 1
                                self.imagematrix[tempx][tempy - 1].label = maxlabel
                                point = [tempx,tempy - 1]
                                points.append(point)
            self.maxlabel = maxlabel
            self.edgesp = edgesp
######################################################
# holefill: find the clusters that has only one point,
# and get rid off it. If the sets value of a point is 
# different from the left right, and top bottom points 
# to it, It is a cluster with only one point.
######################################################
    def holefill(self):
        xmax = self.xmax
        ymax = self.ymax
        for x in range(xmax + 1):
            for y in range(ymax + 1):
                point = self.imagematrix[x][y]
                leftc = 0
                rightc = 0
                upc = 0
                downc = 0
                if x == 0 and y == 0:
                    ll = self.imagematrix[xmax][y].sets
                    rl = self.imagematrix[x + 1][y].sets
                    ul = self.imagematrix[x][y + 1].sets
                    dl = self.imagematrix[x][ymax].sets
                    if ll != point.sets:
                        leftc = 1
                    if rl != point.sets:
                        rightc = 1
                    if ul != point.sets:
                        upc = 1
                    if dl != point.sets:
                        downc = 1
                if 0 < x < xmax and y == 0:
                    ll = self.imagematrix[x - 1][y].sets
                    rl = self.imagematrix[x + 1][y].sets
                    ul = self.imagematrix[x][y + 1].sets
                    dl = self.imagematrix[x][ymax].sets
                    if ll != point.sets:
                        leftc = 1
                    if rl != point.sets:
                        rightc = 1
                    if ul != point.sets:
                        upc = 1
                    if dl != point.sets:
                        downc = 1
                if x == xmax and y == 0:
                    ll = self.imagematrix[x - 1][y].sets
                    rl = self.imagematrix[0][y].sets
                    ul = self.imagematrix[x][y + 1].sets
                    dl = self.imagematrix[x][ymax].sets
                    if ll != point.sets:
                        leftc = 1
                    if rl != point.sets:
                        rightc = 1
                    if ul != point.sets:
                        upc = 1
                    if dl != point.sets:
                        downc = 1
                if x == 0 and 0 < y < ymax:
                    ll = self.imagematrix[xmax][y].sets
                    rl = self.imagematrix[x + 1][y].sets
                    ul = self.imagematrix[x][y + 1].sets
                    dl = self.imagematrix[x][y - 1].sets
                    if ll != point.sets:
                        leftc = 1
                    if rl != point.sets:
                        rightc = 1
                    if ul != point.sets:
                        upc = 1
                    if dl != point.sets:
                        downc = 1
                if 0 < x < xmax and 0 < y < ymax:
                    ll = self.imagematrix[x - 1][y].sets
                    rl = self.imagematrix[x + 1][y].sets
                    ul = self.imagematrix[x][y + 1].sets
                    dl = self.imagematrix[x][y - 1].sets
                    if ll != point.sets:
                        leftc = 1
                    if rl != point.sets:
                        rightc = 1
                    if ul != point.sets:
                        upc = 1
                    if dl != point.sets:
                        downc = 1
                if x == xmax and 0 < y < ymax:
                    ll = self.imagematrix[x - 1][y].sets
                    rl = self.imagematrix[0][y].sets
                    ul = self.imagematrix[x][y + 1].sets
                    dl = self.imagematrix[x][y - 1].sets
                    if ll != point.sets:
                        leftc = 1
                    if rl != point.sets:
                        rightc = 1
                    if ul != point.sets:
                        upc = 1
                    if dl != point.sets:
                        downc = 1
                if x == 0 and y == ymax:
                    ll = self.imagematrix[xmax][y].sets
                    rl = self.imagematrix[x + 1][y].sets
                    ul = self.imagematrix[x][0].sets
                    dl = self.imagematrix[x][y - 1].sets
                    if ll != point.sets:
                        leftc = 1
                    if rl != point.sets:
                        rightc = 1
                    if ul != point.sets:
                        upc = 1
                    if dl != point.sets:
                        downc = 1
                if 0 < x < xmax and y == ymax:
                    ll = self.imagematrix[x - 1][y].sets
                    rl = self.imagematrix[x + 1][y].sets
                    ul = self.imagematrix[x][0].sets
                    dl = self.imagematrix[x][y - 1].sets
                    if ll != point.sets:
                        leftc = 1
                    if rl != point.sets:
                        rightc = 1
                    if ul != point.sets:
                        upc = 1
                    if dl != point.sets:
                        downc = 1
                if x == xmax and y == ymax:
                    ll = self.imagematrix[x - 1][y].sets
                    rl = self.imagematrix[0][y].sets
                    ul = self.imagematrix[x][0].sets
                    dl = self.imagematrix[x][y - 1].sets
                    if ll != point.sets:
                        leftc = 1
                    if rl != point.sets:
                        rightc = 1
                    if ul != point.sets:
                        upc = 1
                    if dl != point.sets:
                        downc = 1
                if leftc == 1 and upc == 1 and rightc == 1 and downc == 1:
                    if point.x == 0:
                        point.sets = self.imagematrix[x + 1][y].sets
                    else:
                        point.sets = self.imagematrix[x - 1][y].sets
######################################################
# edgecollect: collect the edge points for every cluster
######################################################
    def edgecollect(self,point):
        xmax = self.xmax
        ymax = self.ymax
        edge = []
        newpoint = point
        prep = 0
        while newpoint not in edge:
            edge.append(newpoint)
            x = newpoint.x
            y = newpoint.y
            position = self.edgedetect(newpoint)
            if position[1] == 5:
                if prep == 1:
                    position[1] = 13
                if prep == 2:
                    position[1] = 7
                if prep == 4:
                    position[1] = 7
                if prep == 8:
                    position[1] = 13
            if position[1] == 10:
                if prep == 1:
                    position[1] = 11
                if prep == 2:
                    position[1] = 11
                if prep == 4:
                    position[1] = 14
                if prep == 8:
                    position[1] = 14
            if position[1] == 1 or position[1] == 3 or position[1] == 7:
                prep = 1
                if position[0] == 1:
                    if newpoint.label == self.imagematrix[xmax][ymax].label:
                        newpoint = self.imagematrix[xmax][ymax]
                    else:
                        newpoint = self.imagematrix[xmax][y]
                if position[0] == 2 or position[0] == 3:
                    if newpoint.label == self.imagematrix[x - 1][ymax].label:
                        newpoint = self.imagematrix[x - 1][ymax]
                    else:
                        newpoint = self.imagematrix[x - 1][y]
                if position[0] == 4 or position[0] == 7:
                    if newpoint.label == self.imagematrix[xmax][y - 1].label:
                        newpoint = self.imagematrix[xmax][y - 1]
                    else:
                        newpoint = self.imagematrix[xmax][y]
                if position[0] == 5 or position[0] == 6 or position[0] == 8 or position[0] == 9:
                    if newpoint.label == self.imagematrix[x - 1][y - 1].label:
                        newpoint = self.imagematrix[x - 1][y - 1]
                    else:
                        newpoint = self.imagematrix[x - 1][y]
            if position[1] == 2 or position[1] == 6 or position[1] == 14:
                prep = 2
                if position[0] == 1 or position[0] == 2:
                    if newpoint.label == self.imagematrix[x + 1][ymax].label:
                        newpoint = self.imagematrix[x + 1][ymax]
                    else:
                        newpoint = self.imagematrix[x][ymax]
                if position[0] == 3:
                    if newpoint.label == self.imagematrix[0][ymax].label:
                        newpoint = self.imagematrix[0][ymax]
                    else:
                        newpoint = self.imagematrix[xmax][ymax]
                if position[0] == 4 or position[0] == 5 or position[0] == 7 or position[0] == 8:
                    if newpoint.label == self.imagematrix[x + 1][y - 1].label:
                        newpoint = self.imagematrix[x + 1][y - 1]
                    else:
                        newpoint = self.imagematrix[x][y - 1]
                if position[0] == 6 or position[0] == 9:
                    if newpoint.label == self.imagematrix[0][y - 1].label:
                        newpoint = self.imagematrix[0][y - 1]
                    else:
                        newpoint = self.imagematrix[x][y - 1]
            if position[1] == 4 or position[1] == 12 or position[1] == 13:
                prep = 4
                if position[0] == 1 or position[0] == 2 or position[0] == 4 or position[0] == 5:
                    if newpoint.label == self.imagematrix[x + 1][y + 1].label:
                        newpoint = self.imagematrix[x + 1][y + 1]
                    else:
                        newpoint = self.imagematrix[x + 1][y]
                if position[0] == 3 or position[0] == 6:
                    if newpoint.label == self.imagematrix[0][y + 1].label:
                        newpoint = self.imagematrix[0][y + 1]
                    else:
                        newpoint = self.imagematrix[0][y]
                if position[0] == 7 or position[0] == 8:
                    if newpoint.label == self.imagematrix[x + 1][0].label:
                        newpoint = self.imagematrix[x + 1][0]
                    else:
                        newpoint = self.imagematrix[x + 1][y]
                if position[0] == 9:
                    if newpoint.label == self.imagematrix[0][0].label:
                        newpoint = self.imagematrix[0][0]
                    else:
                        newpoint = self.imagematrix[0][y]

            if position[1] == 8 or position[1] == 9 or position[1] == 11:
                prep = 8
                if position[0] == 1 or position[0] == 4:
                    if newpoint.label == self.imagematrix[xmax][y + 1].label:
                        newpoint = self.imagematrix[xmax][y + 1]
                    else:
                        newpoint = self.imagematrix[x][y + 1]
                if position[0] == 2 or position[0] == 3 or position[0] == 5 or position[0] == 6:
                    if newpoint.label == self.imagematrix[x - 1][y + 1].label:
                        newpoint = self.imagematrix[x - 1][y + 1]
                    else:
                        newpoint = self.imagematrix[x][y + 1]
                if position[0] == 7:
                    if newpoint.label == self.imagematrix[xmax][0].label:
                        newpoint = self.imagematrix[xmax][0]
                    else:
                        newpoint = self.imagematrix[x][0]
                if position[0] == 8 or position[0] == 9:
                    if newpoint.label == self.imagematrix[x - 1][0].label:
                        newpoint = self.imagematrix[x - 1][0]
                    else:
                        newpoint = self.imagematrix[x][0]
        return edge
######################################################
# edgedetect: detect edges
######################################################
    def edgedetect(self,point):
        xmax = self.xmax
        ymax = self.ymax
        x = point.x
        y = point.y
        pp = 0
        leftc = 0
        rightc = 0
        upc = 0
        downc = 0
        position = 0
        if x == 0 and y == 0:
            position = 1
            ll = self.imagematrix[xmax][y].label
            rl = self.imagematrix[x + 1][y].label
            ul = self.imagematrix[x][y + 1].label
            dl = self.imagematrix[x][ymax].label
            if ll != point.label:
                leftc = 1
            if rl != point.label:
                rightc = 1
            if ul != point.label:
                upc = 1
            if dl != point.label:
                downc = 1
        if 0 < x < xmax and y == 0:
            position = 2
            ll = self.imagematrix[x - 1][y].label
            rl = self.imagematrix[x + 1][y].label
            ul = self.imagematrix[x][y + 1].label
            dl = self.imagematrix[x][ymax].label
            if ll != point.label:
                leftc = 1
            if rl != point.label:
                rightc = 1
            if ul != point.label:
                upc = 1
            if dl != point.label:
                downc = 1
        if x == xmax and y == 0:
            position = 3
            ll = self.imagematrix[x - 1][y].label
            rl = self.imagematrix[0][y].label
            ul = self.imagematrix[x][y + 1].label
            dl = self.imagematrix[x][ymax].label
            if ll != point.label:
                leftc = 1
            if rl != point.label:
                rightc = 1
            if ul != point.label:
                upc = 1
            if dl != point.label:
                downc = 1
        if x == 0 and 0 < y < ymax:
            position = 4
            ll = self.imagematrix[xmax][y].label
            rl = self.imagematrix[x + 1][y].label
            ul = self.imagematrix[x][y + 1].label
            dl = self.imagematrix[x][y - 1].label
            if ll != point.label:
                leftc = 1
            if rl != point.label:
                rightc = 1
            if ul != point.label:
                upc = 1
            if dl != point.label:
                downc = 1
        if 0 < x < xmax and 0 < y < ymax:
            position = 5
            ll = self.imagematrix[x - 1][y].label
            rl = self.imagematrix[x + 1][y].label
            ul = self.imagematrix[x][y + 1].label
            dl = self.imagematrix[x][y - 1].label
            if ll != point.label:
                leftc = 1
            if rl != point.label:
                rightc = 1
            if ul != point.label:
                upc = 1
            if dl != point.label:
                downc = 1
        if x == xmax and 0 < y < ymax:
            position = 6
            ll = self.imagematrix[x - 1][y].label
            rl = self.imagematrix[0][y].label
            ul = self.imagematrix[x][y + 1].label
            dl = self.imagematrix[x][y - 1].label
            if ll != point.label:
                leftc = 1
            if rl != point.label:
                rightc = 1
            if ul != point.label:
                upc = 1
            if dl != point.label:
                downc = 1
        if x == 0 and y == ymax:
            position = 7
            ll = self.imagematrix[xmax][y].label
            rl = self.imagematrix[x + 1][y].label
            ul = self.imagematrix[x][0].label
            dl = self.imagematrix[x][y - 1].label
            if ll != point.label:
                leftc = 1
            if rl != point.label:
                rightc = 1
            if ul != point.label:
                upc = 1
            if dl != point.label:
                downc = 1
        if 0 < x < xmax and y == ymax:
            position = 8
            ll = self.imagematrix[x - 1][y].label
            rl = self.imagematrix[x + 1][y].label
            ul = self.imagematrix[x][0].label
            dl = self.imagematrix[x][y - 1].label
            if ll != point.label:
                leftc = 1
            if rl != point.label:
                rightc = 1
            if ul != point.label:
                upc = 1
            if dl != point.label:
                downc = 1
        if x == xmax and y == ymax:
            position = 9
            ll = self.imagematrix[x - 1][y].label
            rl = self.imagematrix[0][y].label
            ul = self.imagematrix[x][0].label
            dl = self.imagematrix[x][y - 1].label
            if ll != point.label:
                leftc = 1
            if rl != point.label:
                rightc = 1
            if ul != point.label:
                upc = 1
            if dl != point.label:
                downc = 1
        if leftc == 0 and upc == 0 and rightc == 0 and downc == 0:	#0	0000
            pp = 0
        elif leftc == 0 and upc == 0 and rightc == 0 and downc == 1:#1	1101
            pp = 1
        elif leftc == 0 and upc == 0 and rightc == 1 and downc == 0:#2	1000
            pp = 2
        elif leftc == 0 and upc == 0 and rightc == 1 and downc == 1:#3	1100
            pp = 3
        elif leftc == 0 and upc == 1 and rightc == 0 and downc == 0:#4	1110
            pp = 4
        elif leftc == 0 and upc == 1 and rightc == 0 and downc == 1:#5	0100
            pp = 5
        elif leftc == 0 and upc == 1 and rightc == 1 and downc == 0:#6	0110
            pp = 6
        elif leftc == 0 and upc == 1 and rightc == 1 and downc == 1:#7	0111
            pp = 7
        elif leftc == 1 and upc == 0 and rightc == 0 and downc == 0:#8	0010
            pp = 8
        elif leftc == 1 and upc == 0 and rightc == 0 and downc == 1:#9	0011
            pp = 9
        elif leftc == 1 and upc == 0 and rightc == 1 and downc == 0:#10	1011
            pp = 10
        elif leftc == 1 and upc == 0 and rightc == 1 and downc == 1:#11	0001
            pp = 11
        elif leftc == 1 and upc == 1 and rightc == 0 and downc == 0:#12	1001
            pp = 12
        elif leftc == 1 and upc == 1 and rightc == 0 and downc == 1:#13	1111
            pp = 13
        elif leftc == 1 and upc == 1 and rightc == 1 and downc == 0:
            pp = 14
        elif leftc == 1 and upc == 1 and rightc == 1 and downc == 1:
            pp = 15
        return [position,pp]