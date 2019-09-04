'''
Created on 21 Oct 2018

@author: thomasgumbricht
'''
#imports
import numpy as np
from os import path, makedirs
import array as arr
import subprocess
from sys import exit
import cmd
from geoimagine.support.webcolors import ReadColorHexXML

class MovieClock:
    """
    Constructor
    """
    def __init__(self, process, paramsD, locus):
        """
        Expects procParams from the kartturXML package.
        """
        booleans = ['textatclock','textinvideo']
        for b in booleans:
            if paramsD[b] in ['','N']:
                paramsD[b] = False
        self.params = lambda:None
        for item in paramsD:
            setattr(self.params, item, paramsD[item])

        self.overwrite = process.process.overwrite
        self.period = process.process.dstperiod
        self.dstLayers = process.process.dstLayerD
 
        self._SetItems()

        if self.timeline and self.params.clframecolor == self.params.clhandcolor:
            exit('The color coding for clframecolor and clhandcolor must differ')
        if self.clock and self.params.tlbordercolor == self.params.tlcolor:
            exit('The color coding for clframecolor and clhandcolor must differ')
   
        self._SetSizes()
        
        self._SetPositions()

        self._GetTiming()

        if self.timeline:
            self._SetYearTextGeometry()

        if self.clock:
            self._clBorder()

        self.colorDict, self.colorList = ReadColorHexXML()

        self._ImageFrame()

        self._TimeLoop(locus)
       
    def _SetItems(self):
        self.timeline = True
        self.clock = True

        if self.period.timestep[0:8] in ['seasonal']:
            self.timeline = False
            self.params.tlmargin = [0,0,0,0]
            self.params.tlticks = 0
        if self.period.pdTS.annualperiods == 1:
            self.clock = False
            self.params.clmargin = [0,0,0,0]
            
    def _clBorder(self):
        """
        Creates a mask circle - used as the rim of the annual cl.
        """
        if self.params.clradius-self.params.clborder <= 10:
            exit('The clborder is too thick comapared to the clradius, either increase clradius or decrease clborder')
        cx, cy = self.params.clradius,self.params.clradius
        x,y = np.ogrid[:self.params.clradius*2,:self.params.clradius*2]
        tmin,tmax = np.deg2rad((0,360))
        # convert cartesian --> polar coordinates
        r2 = (x-cx)*(x-cx) + (y-cy)*(y-cy)
        theta = np.arctan2(x-cx,y-cy) - tmin 
        # wrap angles between 0 and 2*pi
        theta %= (2*np.pi)
        # circular outer mask
        circmask = r2 <= self.params.clradius**2
        # angular mask
        anglemask = theta <= (tmax-tmin)
        # circular inner mask
        self.innercircmask = r2 <= (self.params.clradius-self.params.clborder)**2
        self.clrim = circmask*anglemask*(~self.innercircmask)

    def _SetSizes(self):
        """
        Assigns the columns and rows of the cl and the tl
        """
        if self.clock:
            #margins left, right, top, bottom 
            self.clWidthMargins = self.params.clradius*2+self.params.clmargin[0]+self.params.clmargin[1]
            self.clHeightMargins = self.params.clradius*2+self.params.clmargin[2]+self.params.clmargin[3]
        else:
            self.clWidthMargins = self.clHeightMargins =  0
        if self.timeline:
            self.tlWidthMargins = self.params.width - self.params.tlmargin[0]+self.params.tlmargin[1]
            self.tlHeightMargins = self.params.tlheight+self.params.tlmargin[2]+self.params.tlmargin[3]
        else:
            self.tlWidthMargins = self.tlHeightMargins =  0
        self.imgCols = max(self.clWidthMargins,self.tlWidthMargins)
        if self.timeline:
            self.imgRows = self.clHeightMargins + self.tlHeightMargins +  self.params.fontsize+1
        else:
            self.imgRows = self.clHeightMargins + self.tlHeightMargins 
        #For the video codec to work, the hight and width should be even 
        if not self.imgCols % 2 == 0:
            self.imgCols -= 1
        if not self.imgRows % 2 == 0:
            self.imgRows += 1
            
    def _SetPositions(self): 
        """
        Assigns the upper left and lower right corners of the cl and the tl.
        """
        if self.params.position == 'ul':
            if self.timeline:    
                self.tlul = (self.params.tlmargin[0]+self.params.tlborder,
                             self.params.tlmargin[2]+self.params.tlborder)
                
                self.tllr = (self.imgCols-self.params.tlmargin[1]-self.params.tlborder,  
                             self.params.tlmargin[2]+self.params.tlheight-self.params.tlborder)
            
            self.clCenter = [self.params.clradius+self.params.clmargin[0] , 
                             self.params.clradius + self.tlHeightMargins+self.params.clmargin[2]]  
            if self.params.textatclock:
                self.clCenter[1] += self.params.fontsize
                self.clCenter[1] -= self.params.tlmargin[3]
                
        elif self.params.position == 'ur':
            if self.timeline:      
                self.tlul = (self.params.tlmargin[0]+self.params.tlborder,
                             self.params.tlmargin[2]+self.params.tlborder)
                
                self.tllr = (self.imgCols-self.params.tlmargin[1]-self.params.tlborder,  
                             self.params.tlmargin[2]+self.params.tlheight-self.params.tlborder)
            if self.clock:    
                self.clCenter = [self.imgCols-self.params.clradius-self.params.clmargin[1] , 
                                 self.params.clradius + self.tlHeightMargins+self.params.clmargin[2]] 
                if self.params.textatclock:
                    self.clCenter[1] += self.params.fontsize
                    self.clCenter[1] -= self.params.tlmargin[3]
                
        if self.params.position == 'll':
            if self.timeline: 
                self.tlul = (self.params.tlmargin[0]+self.params.tlborder, 
                             self.imgRows - self.params.tlmargin[3]-self.params.tlheight+self.params.tlborder)
                
                self.tllr = (self.imgCols-self.params.tlmargin[1]-self.params.tlborder, 
                             self.imgRows - self.params.tlmargin[3] -self.params.tlborder)
            else:
                self.tlul = (0, self.imgRows)

            if self.clock:
                self.clCenter = [self.params.clradius+self.params.clmargin[0],
                                 self.tlul[1]-self.params.clradius-self.params.tlmargin[2]-self.params.clmargin[3]] 
                if self.params.textatclock:
                    self.clCenter[1] -= self.params.fontsize
                    self.clCenter[1] += self.params.tlmargin[2]
                
        elif self.params.position == 'lr':
            if self.timeline: 
                self.tlul = (self.params.tlmargin[0]+self.params.tlborder, 
                             self.imgRows - self.params.tlmargin[3]-self.params.tlheight+self.params.tlborder)
                
                self.tllr = (self.imgCols-self.params.tlmargin[1]-self.params.tlborder, 
                             self.imgRows - self.params.tlmargin[3] -self.params.tlborder)
            if self.clock:
                self.clCenter = [self.imgCols-self.params.clradius-self.params.clmargin[1], 
                                 self.tlul[1]-self.params.clradius-self.params.tlmargin[2]-self.params.clmargin[3]] 
                if self.params.textatclock:
                    self.clCenter[1] -= self.params.fontsize
                    self.clCenter[1] += self.params.tlmargin[2]
            
        if self.timeline:        
            self.tlheight = self.params.tlheight-2*self.params.tlborder 
    
            if self.tlheight <= 0:
                exit('The tlborder is thicker than your tl, either increase tlheight or decrease tlborder')
        if self.clock:
            self.clul = (self.clCenter[0]-self.params.clradius, self.clCenter[1]-self.params.clradius)
            self.cllr = (self.clCenter[0]+self.params.clradius, self.clCenter[1]+self.params.clradius)   
            self.clSize = (self.params.clradius*2,self.params.clradius*2)
        
    def _ImageFrame(self):
        '''
        Creates the initial time-block frame, with cl, clframe, clrim and tl holder (optional).
        '''
        self.frame = np.ones([self.imgCols,self.imgRows])        
        self.frame*= self.colorDict[self.params.bgcolor] 
        if self.clock:

            self.frame[ self.clul[0]:self.cllr[0], self.clul[1]:self.cllr[1] ] = self.colorDict[self.params.clframecolor]
            self.frame[ self.clul[0]:self.cllr[0], self.clul[1]:self.cllr[1] ][self.clrim] = self.colorDict[self.params.clbordercolor]
            self.frame[ self.clul[0]:self.cllr[0], self.clul[1]:self.cllr[1] ][self.innercircmask] = self.colorDict[self.params.boettcolor]

        if self.timeline and self.params.tlborder:

            ul = (self.tlul[0]-self.params.tlborder, self.tlul[1]-self.params.tlborder )
            lr = (self.tllr[0]+self.params.tlborder, self.tllr[1]+self.params.tlborder )
            self.frame[ ul[0]:lr[0], ul[1]:lr[1] ] = self.colorDict[self.params.tlbordercolor]

    def _GetTiming(self):
        '''
        Sets the clockhand (in degrees) adn the step in the timepline for each timestep
        '''
        if self.clock:
            self.sectorAngle = 360/self.period.pdTS.annualperiods
        if self.timeline:
            self.stepLength = (self.tllr[0]-self.tlul[0])/float(len(self.period.datumL))
            self.tllength = int(round(self.stepLength*len(self.period.datumL)))
        
    def _SetClock(self,angleRange):
        '''
        Sets the annual cl sector in each timestep.
        '''
        x,y = np.ogrid[:self.params.clradius*2,:self.params.clradius*2] 
        cx,cy = (self.params.clradius,self.params.clradius)
        tmin,tmax = np.deg2rad(angleRange)
        tmin -= np.pi/2
        tmax -= np.pi/2
        # ensure stop angle > start angle
        if tmax < tmin:
                tmax += 2*np.pi
        # convert cartesian --> polar coordinates
        r2 = (x-cx)*(x-cx) + (y-cy)*(y-cy)
        theta = np.arctan2(x-cx,y-cy) - tmin
        # wrap angles between 0 and 2*pi
        theta %= (2*np.pi)
        # circular mask
        circmask = r2 <= (self.params.clradius-self.params.clborder)**2
        # angular mask
        anglemask = theta <= (tmax-tmin)
        self.clhand = circmask*anglemask
        
    def _SetYearTextGeometry(self):
        from math import cos,radians,sqrt
        #Set the years

        self.period.pdTS._SetYear()
        if len(self.period.pdTS.yearL) < 2:
            self.yearlength = self.yearheight = 0 
            return

        self.rotationfactor = 1/cos(radians(self.params.rotate))
        self.yearlength = min(int(self.stepLength*self.period.pdTS.annualperiods*self.params.tlticks*4/5*self.rotationfactor ), self.params.fontsize*2)
        self.yearheight = min( self.params.fontsize,int(round(self.yearlength/2)) )
        self.textyshift = self.rotationfactor*(sqrt( (self.yearlength/2)**2 + (self.yearheight/2)**2) )
        self.textyshift -= self.yearheight/2  

        if self.params.textatclock:
            pass
           
    def _WriteYearCmd(self, labelD): 
        if self.params.rotate:
            cmd = " \( -size %(w)dx%(h)d! -background %(fontbg)s \
             -font %(font)s -fill %(fill)s caption:'%(y)d' \
             -trim -rotate '%(r)d' \)  -gravity northwest -geometry \
             +%(ulx)d+%(uly)d -composite " %labelD
        else:
            cmd = " \( -size %(w)dx%(h)d! -background %(fontbg)s \
             -font %(font)s -fill %(fill)s caption:'%(y)d' \
             -trim -rotate '%(r)d' \)  -gravity northwest -geometry \
             +%(ulx)d+%(uly)d -composite " %labelD
        return cmd

    def _WriteYearCmdlabel(self, labelD): 
        cmd = " \( -size %(w)dx%(h)d -background %(fontbg)s \
         -font %(font)s -fill %(fill)s label:'%(y)d'  \
         -trim -gravity center -extent %(w)dx%(h)d -rotate '%(r)d' \) \
         -geometry +%(ulx)d+%(uly)d -composite " %labelD

        return cmd

    def _TimeLoop(self,locus):
        '''
        Loop over all timesteps and create one png image per timestep
        '''
        #only one key for movieclock
        #for locus in self.dstLayers:
        for x,datum in enumerate(self.dstLayers[locus]):
            for comp in self.dstLayers[locus][datum]:
                if self.dstLayers[locus][datum][comp]._Exists() and not self.overwrite:
                    pass
                else:
                    tmpclockFPN = path.join(self.dstLayers[locus][datum][comp].FP,'tmp.png')
                    self._DoIt(tmpclockFPN,x,self.dstLayers[locus][datum][comp].FPN)    
                                        
    def _DoIt(self,tmpclockFPN,x, mcFPN):
        self._ImageFrame()
          
        angleRange = ((self.period.datumD[self.period.datumL[x]]['season']-1)*self.sectorAngle,self.period.datumD[self.period.datumL[x]]['season']*self.sectorAngle)
        #set the cl
        if self.clock:
            self._SetClock(angleRange)
            #transpose required to align cl properly with frame
            trans = np.transpose(self.clhand)
        
            #Set the clhand in the frame
            self.frame[ self.clul[0]:self.cllr[0], self.clul[1]:self.cllr[1] ][trans] = self.colorDict[self.params.clhandcolor]

        if self.timeline:
            #Set the tl
            ul = (self.tlul[0], self.tlul[1] )
            lr = (self.tlul[0]+int(round((x+1)*self.stepLength)),self.tllr[1])
            self.frame[ ul[0]:lr[0], ul[1]:lr[1] ] = self.colorDict[self.params.tlcolor] 
            
            #Write the years
            firstseason = self.period.datumD[self.period.datumL[0]]['season']
            if firstseason == 1:
                self.firsttick = 0
            else:
                self.firsttick = (self.period.pdTS.annualperiods-firstseason+1)*self.stepLength
            if self.params.tlticks:
                yearD = {}
                y = 0
                magickCmd = 'convert %(src)s ' %{'src':tmpclockFPN}
                while True:
                    ystep = y*self.params.tlticks
                    xc = int(round(ystep*self.period.pdTS.annualperiods*self.stepLength+self.firsttick))

                    if self.params.position in ['ul','ur']:
                        ul = [self.tlul[0]+xc, self.tlul[1]+4+self.params.tlmargin[3]]
                        #depress if rotate
                        if self.params.rotate:
                            ul[1] += self.textyshift
                            pass
                    else:
                        ul = [self.tlul[0]+xc, self.tlul[1]-4-self.yearheight]
                        #lift if rotate
                        if self.params.rotate:
                            ul[1] -= self.textyshift
                            pass
                    lr = (ul[0] + self.yearlength,ul[1]+self.yearheight)
                    if lr[0] >= self.tllength:
                        break
                    y += 1
                    if ystep >= len(self.period.pdTS.yearL):
                        break
                    if not self.params.textatclock:
                        if self.clul[0] < ul[0] < self.cllr[0] or self.clul[0] < lr[0] < self.cllr[0]:
                            pass
                        else:
                            yearD= {'y':self.period.pdTS.yearL[ystep], 'font':self.params.font, 'fill': self.params.fontcolor, 
                                                       'fontbg':self.params.fontbackground , 'w':self.yearlength, 'h':self.yearheight, 
                                                       'ulx':ul[0],'uly':ul[1], 'r':self.params.rotate}
                            cmd = self._WriteYearCmd(yearD)
                            magickCmd = '%(a)s %(b)s ' %{'a':magickCmd, 'b':cmd}
                    else:
                        yearD= {'y':self.period.pdTS.yearL[ystep], 'font':self.params.font, 'fill': self.params.fontcolor,
                                                    'fontbg':self.params.fontbackground , 'w':self.yearlength, 'h':self.yearheight, 
                                                    'ulx':ul[0],'uly':ul[1], 'r':self.params.rotate}
    
                        #Add the text if inimage
                        if not self.params.textinvideo:
                            cmd = self._WriteYearCmd(yearD)
                            magickCmd = '%(a)s %(b)s ' %{'a':magickCmd, 'b':cmd}
                                            
            magickCmd = '%(a)s -transparent %(trans)s %(dst)s ' %{'a':magickCmd, 'trans':self.params.transparent, 'dst':mcFPN}

            if self.params.tlticks:
                y = 0
                while True:
                    ystep = y*self.params.tlticks
                    xc = int(round(ystep*self.period.pdTS.annualperiods*self.stepLength+self.firsttick))
 
                    ul = (self.tlul[0]+xc-self.params.tltickwidth,self.tlul[1])
                    if ul[0] >= self.tllength:
                        break
                    lr = (self.tlul[0]+xc+self.params.tltickwidth,self.tllr[1])
                    if lr[0] > self.tllength+self.params.tlborder:
                        break
                    self.frame[ ul[0]:lr[0], ul[1]:lr[1] ] = self.colorDict[self.params.tickcolor]
                    y += 1

        flatA = self.frame.flatten('F')
        INDEX = flatA.tolist()
                    
        if not self.imgCols*self.imgRows == len(INDEX):

            exit('internal error in image size definition')

        if self.params.tlticks:
            if self.params.textinvideo:
                WritePng(INDEX, mcFPN, self.imgCols, self.imgRows, self.colorList)      
            else:
                WritePng(INDEX,tmpclockFPN,self.imgCols,self.imgRows,self.colorList)
                cr = subprocess.call('/usr/local/bin/' + magickCmd, shell=True)
                if cr:
                    print ('shell', subprocess.call('/usr/local/bin/' + magickCmd, shell=True) )
        else:
            WritePng(INDEX, mcFPN, self.imgCols, self.imgRows, self.colorList)
                        
def WritePng(INDEX,imageFPN,cols,lins,paletteL):
    import png as img_png
    #write the image file
    RGB = arr.array('B',(int(i) for i in INDEX))
    f = open(imageFPN, 'wb')      # binary mode is important    
    w = img_png.Writer(width=cols, height=lins, greyscale=False, bitdepth=8, palette=paletteL, chunk_limit = cols)
    w.write_array(f, RGB)
    f.close()
