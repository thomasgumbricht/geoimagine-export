'''
Created on 21 Oct 2018

@author: thomasgumbricht
'''
import os
from sys import exit
import numpy as np
import subprocess
from geoimagine.gdalutilities import GDALstuff
from geoimagine.ktgraphics import MapPlot

class ProcessExport():
    def __init__(self, process, session, verbose):
        #self.session = session
        self.verbose = verbose
        self.process = process   
        self.session = session 
        printstr = '    Starting ProcessExport: %(s)s' %{'s':self.process.proc.processid}
        print (printstr)
                
        #direct to subprocess
        if 'movieclock' in self.process.proc.processid.lower():
            self._MovieClock()
            if self.process.params.asscript and self.framecriptF:
                print ('        run',self.framescriptFPN)
                print ('        run',self.moviescriptFPN)
                exit('Exiting - you now have to run the shell scripts above')
            return
        #Set scriptF to False, if asscript == True, slef.scriptF gets a value
        self.scriptF = False
        if 'movieframe' in self.process.proc.processid.lower() and self.process.params.asscript:
            print ('self.process.srcIdDict',self.process.srcIdDict)

            self._IniMovieFrameScript()
        if 'movieoverlayframe' in self.process.proc.processid.lower():
            self._IniMovieFrameScript()

        self._LoopAllLayers()
        
        if 'movieframe' in self.process.proc.processid.lower() and self.scriptF:
            self.scriptF.close()
            print ('        run',self.scriptFPN)
        if 'movieoverlayframe' in self.process.proc.processid.lower() and self.scriptF:
            self.scriptF.close()
            print ('        run',self.scriptFPN)
            
        printstr = '    Finished ProcessExport: %(s)s' %{'s':self.process.proc.processid}
        print (printstr)
      

    def _LoopAllLayers(self):
        for srclocus in self.process.srcLayerD:
            if len(self.process.srcLayerD) == 1:
                dstlocus = list(self.process.dstLayerD.keys())[0]
            else:
                dstlocus = srclocus
            print ('        Processing locus',srclocus)
            for datum in self.process.srcLayerD[srclocus]:
                srcCompL = []
                for comp in self.process.srcLayerD[srclocus][datum]:
                    #Check if the file is there

                    if os.path.exists(self.process.srcLayerD[srclocus][datum][comp].FPN):
                        srcCompL.append(comp)
                if len(srcCompL) == 0:
                    SNULLEBULLE
                #Always only one dst comp
                dstcomp = (list(self.process.dstLayerD[dstlocus][datum].keys())[0])

                if 'archive' in self.process.proc.processid.lower():
                    self._ArchiveIni(srclocus,datum,srcCompL)
                elif 'exporttobyte' in self.process.proc.processid.lower(): 
                    self._ExportToByte(srclocus,datum,srcCompL)
                elif not self.process.dstLayerD[dstlocus][datum][dstcomp]._Exists() or self.process.overwrite:
                    if 'append' in self.process.proc.processid.lower():
                        self._MovieAppendFrames(srclocus,datum,comp) 
                    elif 'movieframe' in self.process.proc.processid.lower():
                        self._MovieFrames(srclocus,datum,comp) 
                    elif 'movieoverlayframe' in self.process.proc.processid.lower():
                        self._MovieOverlayFrames(srclocus,datum)
                    elif 'exporttosvg' in self.process.proc.processid.lower():
                        self._ExportSVG(srclocus,dstlocus,datum,srcCompL[0],dstcomp) 
                    elif 'exportmap' in self.process.proc.processid.lower():
                        self._ExportMap(srclocus,datum,srcCompL)                                  
                    else:
                        print (self.process.proc.processid)
                        ADDPROCESS
                else:
                    pass
                    #print ('exists',self.process.dstLayerD[locus][datum][dstcomp].FPN)
         
    def _ArchiveIni(self,locus,datum,srcCompL):
        '''
        '''
        for comp in srcCompL:
            print (self.process.dstLayerD[locus][datum][comp].FPN)
            SNULLE
            if not self.process.dstLayerD[locus][datum][comp]._Exists() or self.process.overwrite:
                if self.process.params.archive == 'zip':
                    self._ArchiveZip(locus,datum,comp)
                elif self.process.params.archive == 'gz':
                    self._archiveGunZip(locus,datum,comp)
                elif self.process.params.archive == 'tar':
                    self._archiveTar(locus,datum,comp)
                elif self.process.params.archive == 'gz.tar':
                    self._archiveGunZipTar(locus,datum,comp)
            
    def _ArchiveZip(self,locus,datum,comp):   
        import zipfile
        try:
            import zlib
            compression = zipfile.ZIP_DEFLATED
        except:
            compression = zipfile.ZIP_STORED
        zf = zipfile.ZipFile(self.process.dstLayerD[locus][datum][comp].FPN, mode='w')
        try:
            zf.write(self.process.srcLayerD[locus][datum][comp].FPN, compress_type=compression, arcname=self.process.srcLayerD[locus][datum][comp].FN)
        finally:
            zf.close()

                  
    def _ExportToByte(self,locus,datum,srcCompL):
        '''
        '''
        #self.SetMask()
        for comp in srcCompL:
            if not self.process.dstLayerD[locus][datum][comp]._Exists() or self.process.overwrite:
                self._SelectScaling(comp)
                self._SetPalette(comp)
                self._ExportLayer(locus,datum,comp)
                self._CreateLayout(locus,datum,comp)
                
    def _ExportMap(self,locus,datum,srcCompL):
        '''
        '''
        #self.SetMask()
        
        for comp in srcCompL:
            if not self.process.dstLayerD[locus][datum][comp]._Exists() or self.process.overwrite:
                self._SelectScaling(comp)

                self._SelectLegend(comp)
                measure = self.process.proc.dstcompD[comp]['measure']
                self._SetPalette(comp)
                mapArr = self._ExportLayer(locus,datum,comp, True)

                print (self.process.dstLayerD[locus][datum][comp].FPN)
                print (self.process.srcLayerD[locus][datum][comp].FPN)

                MapPlot(mapArr, self.process, self.process.dstCompD[comp].colorRamp, self.legend, self.scaling, measure, self.process.srcLayerD[locus][datum][comp], self.process.dstLayerD[locus][datum][comp])

                
    def _SelectLegend(self,comp):
        '''THIS IS A DUPLICATE FROM layout.legend
        Select legend from database
        '''
        legendD = self.session.IniSelectLegend(self.process.proc.dstcompD[comp])

        self.legend = lambda: None
        for key, value in legendD.items():
            setattr(self.legend, key, value)
        self.legend.frame = int(self.legend.framestrokewidth+0.99)
                
    def _ExportSVG(self,srclocus,dstlocus,datum,srccomp,dstcomp):

        src = self.process.srcLayerD[srclocus][datum][srccomp].FPN
        dst = self.process.dstLayerD[dstlocus][datum][dstcomp].FPN
        tempShp = dst.replace('.svg','.shp')
        #Get the bounding box for the region
        regExt = self.session._SelectRegionLonLatExtent(dstlocus,'D')

        #clipping done in ogr2ogr, better control and more stable
        clipVec = GDALstuff(src,tempShp,self.process.params)
        clipVec.SetClipBoxLLminmax(regExt)
        clipVec.ClipVector()

        srcbasename = os.path.splitext( self.process.dstLayerD[dstlocus][datum][dstcomp].FN )[0]
        styleD = {'srcbasename':srcbasename,'f':self.process.params.fill ,'s':self.process.params.stroke,'sw':self.process.params.strokewidth,'sd':self.process.params.strokedasharray}
        style = '.%(srcbasename)s { fill: %(f)s; stroke: %(s)s; stroke-width: %(sw)s; ' %styleD

        if len(self.process.params.strokedasharray) > 2:
            style += 'stroke-dasharray: %(sd)s; }' %styleD
        else:
            style += '}' %styleD
        print ('style',style)

        cmdD ={'epsg':self.process.params.t_epsg, 'src':tempShp, 'dst':dst, 'simp':self.process.params.simplify,
               'scale':self.process.params.scale,'padd':self.process.params.padding,'style':style}
        cmd =  "svgis draw %(src)s --crs EPSG:%(epsg)d -p %(padd)s --style '%(style)s' " %cmdD
        if self.process.params.simplify:
            cmd += "-s %(simp)d " %cmdD
        if self.process.params.scale:
            cmd += "-f %(scale)d " %cmdD
        cmd += "-o %(dst)s" %cmdD
        print (cmd)
        
        subprocess.call('/Users/thomasgumbricht/anaconda3/bin/' + cmd, shell=True)
        #os.system(cmd)

        pngFPN  = dst.replace('.svg','.png')
        
        
        src = dst
        #dst = pngFPN
        if self.process.params.crop != '':
            cropL = self.process.params.crop.split(',')
            cropL = [int(item) for item in cropL]
        else:
            cropL = False

        width, height = self.process.params.pngsize.split(',')
        params = {'w':int(width),'h':int(height), 'src':src, 'dst':pngFPN}

        magickCmd = 'convert -background none -density 72 -resize %(w)dx%(h)d! %(src)s %(dst)s' %params 

        print (magickCmd)
        if self.process.params.asscript:
            magickCmd += '\n'
            self.scriptF.write(magickCmd)
        else:
            subprocess.call('/usr/local/bin/' + magickCmd, shell=True)

    def _IniMovieFrameScript(self):
        #print (layerId = self.process.srcLayerD[locus][datum][srccomp].comp.id)
        locus = list(self.process.dstLayerD.keys())[0]
        datum = list(self.process.dstLayerD[locus].keys())[0]
        comp = self.process.srcIdDict['base']
        if not os.path.exists(self.process.dstLayerD[locus][datum][comp].FP):
            os.makedirs(self.process.dstLayerD[locus][datum][comp].FP)
        if self.process.params.asscript:
            scriptFN = 'images_%(comp)s.sh' %{'comp':comp}
            self.scriptFPN = os.path.join(self.process.dstLayerD[locus][datum][comp].FP,scriptFN)
            self.scriptF = open(self.scriptFPN,'w')
            
        '''
        for locus in self.process.dstLayerD:
            for datum in self.process.dstLayerD[locus]:
                if not os.path.exists(self.process.dstLayerD[locus][datum][comp].FP):
                    os.makedirs(self.process.dstLayerD[locus][datum][comp].FP)
            if self.process.params.asscript:

                scriptFN = 'images_%(comp)s.sh' %{'comp':comp}
                self.scriptFPN = os.path.join(self.process.dstLayerD[locus][datum][comp].FP,scriptFN)

                self.scriptF = open(self.scriptFPN,'w')
        '''
    
    def _MovieFrames(self,locus,datum,comp):
        '''
        ''' 
        src = self.process.srcLayerD[locus][datum][comp].FPN
        dst = self.process.dstLayerD[locus][datum][comp].FPN
        if self.process.params.crop != '':
            cropL = self.process.params.crop.split(',')
            cropL = [int(item) for item in cropL]
        else:
            cropL = False
        border = self.process.params.border
        if len(self.process.params.emboss) > 0:
            emboss = True
        else:
            emboss = False
        if len(self.process.params.vectoroverlay) > 4:
            if os.path.isfile(self.process.params.vectoroverlay):
                overlay = self.process.params.vectoroverlay
            else:
                exit('The vectoroverlay for movieframes does not exist')
        else:
            overlay = False
        embdimL = self.process.params.embossdims.split(',')
        embdimL = [int(item) for item in embdimL]
        params = {'w':self.process.params.width,'cw':cropL[0], 'ch':cropL[1], 'cx':cropL[2], 'cy':cropL[3],
                  'o':overlay, 'b':border,'bc':self.process.params.bordercolor,
                  'emboss':self.process.params.emboss, 'ptsize':self.process.params.embossptsize,
                  'embx': embdimL[0], 'emby':embdimL[1],
                  'src':src, 'dst':dst}
        
        magickCmd = 'convert \( -resize %(w)dx ' %params

        if cropL:
            magickCmd += '-crop %(cw)dx%(ch)d+%(cx)d+%(cy)d ' %params
        magickCmd += ' %(src)s \) ' %params
        
        if overlay:
            if cropL:
                magickCmd += '\( -background none -resize %(cw)dx%(ch)d! %(o)s \) -composite ' %params
            else:
                magickCmd += '\( -background none -resize %(w)dx %(o)s \) -composite ' %params
        
        if emboss:
            magickCmd += '\( -size %(embx)dx%(emby)s xc:none -font Trebuchet -pointsize %(ptsize)s -gravity center -draw "fill silver text 1,1 ' %params
            magickCmd += "'%(emboss)s' fill whitesmoke text -1,-1 '%(emboss)s' fill grey text 0,0 '%(emboss)s' " %params
            magickCmd += '" -transparent grey -fuzz 90%% \) -composite ' %params

        if border:
            magickCmd += '\( -border %(b)dx%(b)d -bordercolor %(bc)s \) ' %params
        magickCmd += '%(dst)s' %params

        #print (magickCmd)
        if self.process.params.asscript:
            magickCmd += '\n'
            self.scriptF.write(magickCmd)
        else:
            subprocess.call('/usr/local/bin/' + magickCmd, shell=True)
                       
    def _MovieOverlayFrames(self,locus,datum):
        '''
        ''' 
        baseSrc = self.process.srcLayerD[locus][datum][self.process.srcIdDict['base']].FPN
        overSrc = self.process.srcLayerD[locus][datum][self.process.srcIdDict['over']].FPN
        dst = self.process.dstLayerD[locus][datum][self.process.srcIdDict['base']].FPN
        cropL = self.process.params.basecrop.split(',')
        cropL = [int(item) for item in cropL]
        embdimL = self.process.params.embossdims.split(',')
        embdimL = [int(item) for item in embdimL]
        params = {'bw':self.process.params.basewidth,'cw':cropL[0], 'ch':cropL[1], 'cx':cropL[2], 'cy':cropL[3],
                  'ow':self.process.params.overlaywidth,'g':self.process.params.gravity,'ox':self.process.params.geomx, 'oy':self.process.params.geomy,
                  'emboss':self.process.params.emboss, 'ptsize':self.process.params.embossptsize,
                  'embx': embdimL[0], 'emby':embdimL[1],
                  'bSrc':baseSrc, 'oSrc':overSrc,'dst':dst}
        magickCmd = 'convert \( -resize %(bw)dx -crop %(cw)dx%(ch)d+%(cx)d+%(cy)d %(bSrc)s \) ' %params
        
        magickCmd += '\( -resize %(ow)dx -border 1x1 -bordercolor black %(oSrc)s \) -gravity %(g)s -geometry +%(ox)d+%(oy)d -composite ' %params
        
        magickCmd += '\( -size %(embx)dx%(emby)s xc:none -font Trebuchet -pointsize %(ptsize)s -gravity center -draw "fill silver text 1,1 ' %params
        magickCmd += "'%(emboss)s' fill whitesmoke text -1,-1 '%(emboss)s' fill grey text 0,0 '%(emboss)s' " %params
        magickCmd += '" -transparent grey -fuzz 90%% \) -composite %(dst)s' %params

        if self.process.params.asscript:
            magickCmd += '\n'
            self.scriptF.write(magickCmd)
        else:
            subprocess.call('/usr/local/bin/' + magickCmd, shell=True)
                 
    def _MovieAppendFrames(self,locus,datum,comp):
        '''
        ''' 
        
        if 'top' in self.process.srcIdDict:
            dst = self.process.dstLayerD[locus][datum][self.process.srcIdDict['top']].FPN
            firstSrc = self.process.srcLayerD[locus][datum][self.process.srcIdDict['top']].FPN
            secondSrc = self.process.srcLayerD[locus][datum][self.process.srcIdDict['bottom']].FPN
            #iAppend = 'horisontal'
            params = {'firstSrc':firstSrc, 'secondSrc':secondSrc, 'dst':dst}
            magickCmd = 'convert -append %(firstSrc)s %(secondSrc)s %(dst)s' %params
        elif 'left' in self.process.srcIdDict:
            dst = self.process.dstLayerD[locus][datum][self.process.srcIdDict['left']].FPN
            firstSrc = self.process.srcLayerD[locus][datum][self.process.srcIdDict['left']].FPN
            secondSrc = self.process.srcLayerD[locus][datum][self.process.srcIdDict['right']].FPN
            #iAppend = 'vertical'
            params = {'firstSrc':firstSrc, 'secondSrc':secondSrc, 'dst':dst}
            magickCmd = 'convert +append %(firstSrc)s %(secondSrc)s %(dst)s' %params
        else:
            NOTYET

        if self.process.params.asscript:
            magickCmd += '\n'
            self.scriptF.write(magickCmd)
        else:
            subprocess.call('/usr/local/bin/' + magickCmd, shell=True)
                     
    def _SelectScaling(self,comp):
        '''
        '''
        # self.process.proc.dstcompD[comp] is the dictionary form of the compostion
        scalingD = self.session.IniSelectScaling(self.process.proc.dstcompD[comp])
        self.scaling = lambda: None 
        for key, value in scalingD.items():
            setattr(self.scaling, key, value)
        '''
        print ('self.scaling.power',self.scaling.power)
        print ('self.scaling.powerna',self.scaling.powerna)
        print ('self.scaling.mirror0',self.scaling.mirror0)
        print ('self.scaling.scalefac',self.scaling.scalefac)
        print ('self.scaling.offsetadd',self.scaling.offsetadd)
        '''
    
    def _SetPalette(self,comp):
        #Set the palette for this composition
        self.process.dstCompD[comp].palette = self.process.params.palette
        #Create the palette
        self.process.dstCompD[comp]._CreatePalette()
        
    def _ExportLayer(self,locus,datum,comp,returnmap=False):
        '''
        '''
        #Open the src layer
        self.process.srcLayerD[locus][datum][comp].ReadSrcLayer()
        srcBAND = self.process.srcLayerD[locus][datum][comp].layer.NPBAND
        #print (srcBAND)
        #print (srcBAND[20,750])
        srcCellNull = self.process.srcLayerD[locus][datum][comp].layer.cellnull
        if self.scaling.power:                 
            if self.scaling.mirror0:
                #BAND = srcBAND+ self.legendD['offsetadd']
                BAND = srcBAND
                dstBANDIni = np.power(BAND, self.scaling.power)
                dstBANDIni *=  self.scaling.scalefac
                dstBANDInv = np.power(-BAND, self.scaling.power)
                dstBANDInv *= - self.scaling.scalefac
                MASK = (BAND < 0)
                dstBAND = np.copy(dstBANDIni)
                dstBAND[MASK] = dstBANDInv[MASK]
                dstBAND += 125
            else:
                #changed to work with VWB deficit 20181125
                BAND = srcBAND
                dstBANDIni = np.power(BAND, self.scaling.power)
                dstBANDIni *=  self.scaling.scalefac
                #dstBANDInv = np.power(-BAND, self.scaling.power)
                #dstBANDInv *= - self.scaling.scalefac
                #MASK = (BAND < 0)
                dstBAND = np.copy(dstBANDIni)
                #dstBAND[MASK] = dstBANDInv[MASK]
                dstBAND += self.scaling.offsetadd
        elif self.scaling.mirror0 or self.scaling.offsetadd == 125: 
            dstBAND =  srcBAND * self.scaling.scalefac
            dstBAND += 125
        else:  
            #Only scale and offset
            dstBAND = srcBAND*self.scaling.scalefac
            dstBAND += self.scaling.offsetadd

        dstBAND[dstBAND < 0] = 0
            
        dstBAND[dstBAND > 250] = 250

        #set nulll to 255
        dstBAND[srcBAND == srcCellNull] = 255
        #set nan to 255
        dstBAND[np.isnan(srcBAND)] = 255
        
        if returnmap:
            return dstBAND
        #Create the dst layer
        self.process.dstLayerD[locus][datum][comp].layer = lambda:None
        #Set the np array as the band
        self.process.dstLayerD[locus][datum][comp].layer.NPBAND = dstBAND
        #copy the geoformat from the src layer
        self.process.dstLayerD[locus][datum][comp].CopyGeoformatFromSrcLayer(self.process.srcLayerD[locus][datum][comp].layer)
        #Create the dst layer
        self.process.dstLayerD[locus][datum][comp].CreateDSWriteRasterArray()

    def _CreateLayout(self,locus,datum,comp):
        '''Create a layout and export to PNG
        '''
        if len(self.process.params.vectoroverlay) > 4:
            if os.path.isfile(self.process.params.vectoroverlay):
                overlay = self.process.params.vectoroverlay
            else:
                exit('The vectoroverlay for movieframes does not exist')
        else:
            overlay = False
        border = self.process.params.border
        
        if len(self.process.params.emboss) > 0:
            emboss = True
        else:
            emboss = False
        
        #only create a layout if either border or overlay is set
        if border or overlay or emboss or self.process.params.width:
            FPN = self.process.dstLayerD[locus][datum][comp].FPN
            #tgto do, change from tif and png to hdrfiletype and datafiletype
            pngFPN = FPN.replace('.tif','.png')
            
            if self.process.params.crop != '':
                cropL = self.process.params.crop.split(',')
                cropL = [int(item) for item in cropL]
            else:
                cropL = False
            
            embdimL = self.process.params.embossdims.split(',')
            embdimL = [int(item) for item in embdimL]
            params = {'w':self.process.params.width,'cw':cropL[0], 'ch':cropL[1], 'cx':cropL[2], 'cy':cropL[3],
                      'o':overlay, 'b':border,'bc':self.process.params.bordercolor,
                      'emboss':self.process.params.emboss, 'ptsize':self.process.params.embossptsize,
                      'embx': embdimL[0], 'emby':embdimL[1],
                      'src':FPN, 'dst':pngFPN}
            
            magickCmd = 'convert \( -resize %(w)dx ' %params
    
            if cropL:
                magickCmd += '-crop %(cw)dx%(ch)d+%(cx)d+%(cy)d ' %params
            magickCmd += ' %(src)s \) ' %params
            
            if overlay:
                if cropL:
                    magickCmd += '\( -background none -resize %(cw)dx%(ch)d! %(o)s \) -composite ' %params
                else:
                    magickCmd += '\( -background none -resize %(w)dx %(o)s \) -composite ' %params
            
            if emboss:
                magickCmd += '\( -size %(embx)dx%(emby)s xc:none -font Trebuchet -pointsize %(ptsize)s -gravity center -draw "fill silver text 1,1 ' %params
                magickCmd += "'%(emboss)s' fill whitesmoke text -1,-1 '%(emboss)s' fill grey text 0,0 '%(emboss)s' " %params
                magickCmd += '" -transparent grey -fuzz 90%% \) -composite ' %params
    
            if border:
                magickCmd += '\( -border %(b)dx%(b)d -bordercolor %(bc)s \) ' %params
            magickCmd += '%(dst)s' %params
            subprocess.call('/usr/local/bin/' + magickCmd, shell=True)
            
            if self.process.params.jpg:
                jpgFPN = FPN.replace('.tif','.jpg')
                params = {'src':pngFPN, 'dst':jpgFPN, 'q':self.process.params.jpg}
                magickCmd = 'convert %(src)s -quality %(q)d %(dst)s ' %params
                subprocess.call('/usr/local/bin/' + magickCmd, shell=True)
                
    def _MakeMovieDirs(self,locus):
        
            for datum in self.process.dstLayerD[locus]:
                #print ('locus, datum',self.dstLayers[locus][datum])
                for comp in self.process.dstLayerD[locus][datum]:
                    #print ('locus, datum, comp',self.dstLayers[locus][datum][comp])
                    if not os.path.exists(self.process.dstLayerD[locus][datum][comp].FP):
                        os.makedirs(self.process.dstLayerD[locus][datum][comp].FP)
                    #step one up to create folders for ready frames and movie
                    FP = os.path.split(self.process.dstLayerD[locus][datum][comp].FP)[0]
                    self.frameFP = os.path.join(FP,'frames')
                    if not os.path.exists(self.frameFP):
                        os.makedirs(self.frameFP)
                    self.movieFP = os.path.join(FP,'movie')
                    if not os.path.exists(self.movieFP):
                        os.makedirs(self.movieFP)
                    movieFN = os.path.splitext(self.process.dstLayerD[locus][datum][comp].FN)[0]

                    #moviedataum is already set in timestep
                    movieFN = movieFN.replace(datum,self.process.dstperiod.moviedatum)
                    self.movieFN = '%(mf)s.mp4' %{'mf':movieFN}
                    
                    #self.movieFN = '%(c)s.mp4' %{'c':self.process.dstLayerD[locus][datum][comp].comp.compid}
                    self.movieFPN = os.path.join(self.movieFP,self.movieFN)
  
                    #Break the loop and return
                    if self.process.params.asscript:
                        moviescriptFN = 'movie_%(comp)s.sh' %{'comp':comp}
                        framescriptFN = 'frame_%(comp)s.sh' %{'comp':comp}
                        self.moviescriptFPN = os.path.join(self.movieFP,moviescriptFN)
                        self.framescriptFPN = os.path.join(self.frameFP,framescriptFN)
                        if not os.path.isfile(self.moviescriptFPN) or self.process.overwrite:
                            self.moviescriptF = open(self.moviescriptFPN,'w')
                        else:
                            self.moviescriptF = False
                        if not os.path.isfile(self.framescriptFPN) or self.process.overwrite:
                            self.framecriptF = open(self.framescriptFPN,'w')
                        else:
                            self.framecriptF = False   
                    return
                                            
    def _MovieClock(self):
        from geoimagine.export import movie_clock
        query = {'name':self.process.params.name}
        paramL = ['tlmargin','clmargin','position','bgcolor','tlborder','clborder','tlbordercolor',
        'clbordercolor','clcolor','clcolor','tlheight','tlticks','tltickwidth','tickcolor',
        'boettcolor','textatclock','tlcolor','clradius','clhandcolor','clframecolor','fontsize', 
        'fontcolor','fontbackground','rotate','font','textinvideo','transparent']
        paramD = self.session._SelectMovieClock(query, paramL)
        
        for locus in self.process.dstLayerD:
            
            self._MakeMovieDirs(locus)
            if os.path.exists(self.movieFPN) and not self.process.overwrite:
                continue
            elif os.path.exists(self.movieFPN):
                os.remove(self.movieFPN)
            paramD['width'] = self.process.params.width
            paramD['movie'] = self.movieFN
            
            print ('winding movieclock',locus)
            movie_clock.MovieClock(self, paramD,locus)

            gravityD = {'ll':'southwest','ul':'northwest','lr':'southeast','ur':'northeast'}
            #Write the command to overlay the clock and the image
            params = {'gravity':gravityD[paramD['position']],'movie':self.movieFN}
            

            if self.process.params.asscript and self.framecriptF:
                cmd = 'for i in ../images/*.png; do convert "../images/$(basename $i)" "../clock/$(basename $i)" -gravity %(gravity)s -composite "$(basename $i)"; done' %params
    
                cd = 'cd %(d)s;\n' %{'d':self.frameFP}
                self.framecriptF.write(cd)
                cmd += '\n'
                self.framecriptF.write(cmd)
                self.framecriptF.close()
            else:
                cmd = 'for i in ../images/*.png; do /usr/local/bin/convert "../images/$(basename $i)" "../clock/$(basename $i)" -gravity %(gravity)s -composite "$(basename $i)"; done' %params
                cmd = 'cd %(d)s; %(cmd)s' %{'d':self.frameFP, 'cmd':cmd}
    
                os.system(cmd)
    
            #cd to where the data are saved
            cd = 'cd %(d)s;' %{'d':self.movieFP}
            if self.process.params.asscript and self.moviescriptF:
                #The movie command
                cmd = "ffmpeg -r 3 -f image2 -pattern_type glob -i '../frames/*.png' -vcodec libx264 -b:v 1M -pix_fmt yuv420p -crf 35 %(movie)s" %params
                cd += '\n'
                self.moviescriptF.write(cd)
                cmd += '\n'
                self.moviescriptF.write(cmd)
                self.moviescriptF.close()
            else:
                cmd = "/usr/local/bin/ffmpeg -r 3 -f image2 -pattern_type glob -i '../frames/*.png' -vcodec libx264 -b:v 1M -pix_fmt yuv420p -crf 35 %(movie)s" %params
                cmd = '%(cd)s %(cmd)s' %{'cd':cd, 'cmd':cmd}
                print ('cmd',cmd)
                os.system(cmd)
            
         