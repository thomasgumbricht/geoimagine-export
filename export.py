'''
Created on 21 Oct 2018

@author: thomasgumbricht
'''
import os
from sys import exit
import numpy as np

class ProcessExport():
    def __init__(self, process, session, verbose):
        #self.session = session
        self.verbose = verbose
        self.process = process   
        self.session = session         
        #direct to subprocess
  
        #Just loop - order dows not matter, process layer by layer
        self._LoopAllLayers()
      
    def _LoopAllLayers(self):
        for locus in self.process.srcLayerD:
            #print ('locus',locus)
            for datum in self.process.srcLayerD[locus]:
                #print ('    datum',datum)
                for comp in self.process.dstLayerD[locus][datum]:
                    print ('      dstcomp',comp)
                for comp in self.process.srcLayerD[locus][datum]:
                    print ('        comp',comp)

                    #Check if the file is there
                    if os.path.exists(self.process.srcLayerD[locus][datum][comp].FPN):
                        print ('cool',self.process.srcLayerD[locus][datum][comp].FPN) 
                        print ('target',self.process.dstLayerD[locus][datum][comp].FPN)
                        if not self.process.dstLayerD[locus][datum][comp]._Exists() or self.process.overwrite:
                            printstr = 'Exporting: %(fpn)s' %{'fpn':self.process.dstLayerD[locus][datum][comp].FPN}

                            if self.process.proc.processid.lower() == 'exporttobyteancillary':
                                self._ExportToByte(locus,datum,comp)
                            else:
                                ADDPROCESS
                        
    def _ExportToByte(self,locus,datum,comp):
        '''
        '''
        #self.SetMask()
        
        self.SelectScaling(comp)
        self._SetPalette(comp)
        self._ExportLayer(locus,datum,comp)
                 
    def SelectScaling(self,comp):
        '''
        '''
        # self.process.proc.dstcompD[comp] is the dictionary form of the compostion
        scalingD = self.session.IniSelectScaling(self.process.proc.dstcompD[comp])
        self.scaling = lambda: None 
        for key, value in scalingD.items():
            setattr(self.scaling, key, value)
        print ('self.scaling.power',self.scaling.power)
        print ('self.scaling.powerna',self.scaling.powerna)
        print ('self.scaling.mirror0',self.scaling.mirror0)
        print ('self.scaling.scalefac',self.scaling.scalefac)
        print ('self.scaling.offsetadd',self.scaling.offsetadd)

    def _SetPalette(self,comp):
        #Set the palette for this composition
        self.process.dstCompD[comp].palette = self.process.params.palette
        #Create the palette
        self.process.dstCompD[comp]._CreatePalette()

    def _ExportLayer(self,locus,datum,comp):
        '''
        '''
        #Open the src layer
        self.process.srcLayerD[locus][datum][comp].ReadSrcLayer()
        srcBAND = self.process.srcLayerD[locus][datum][comp].layer.NPBAND
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
                BAND = srcBAND+self.scaling.offsetadd
                BAND[BAND < 0] = 0
                BAND = np.power(BAND,self.scaling.power)
                BAND *= self.scaling.scalefac
                dstBAND = BAND
        elif self.scaling.mirror0 or self.scaling.offsetadd == 125: 
            #no power involved 
            #This means that zero will end up at 125
            #dstBANDIni = srcBAND*1
            dstBAND =  srcBAND * self.scaling.scalefac
            print ('self.scaling.mirror0')
            #dstBANDInv = srcBAND*-1
            #dstBANDInv *= self.scaling.scalefac
            dstBAND += 125
        else:  
            NOTYET
            dstBAND = srcBAND*1.0 #Multiply to get a new array, not just a copy
            
            BAND = srcBAND+self.scaling.offsetadd
            BAND[BAND < 0] = 0
            BAND *= self.scaling.scalefac
            dstBAND = BAND 

        dstBAND[dstBAND > 250] = 250

        #BAND[BAND < 0] = 0
        if self.scaling.power and not self.scaling.mirror0:
            NOTYEAT
            dstBAND[srcBAND < 0] = self.scaling.powerna
        
        #set all below 0 to 0
        dstBAND[dstBAND < 0] = 0
        #set nulll to 255
        dstBAND[srcBAND == srcCellNull] = 255

        #Apply mask is applicable
        '''
        if len(self.process.xparamTagD['maskin']) == 1:
            if self.process.parameters.masknotnullin == 255:
                dstBAND[self.maskBAND == 0] = 255
            else:
                dstBAND[(self.maskBAND == 0) & (dstBAND != 255)] = self.process.parameters.masknotnullin
        '''

        #Create the dst layer
        self.process.dstLayerD[locus][datum][comp].layer = lambda:None
        

        #Set the np array as the band
        self.process.dstLayerD[locus][datum][comp].layer.NPBAND = dstBAND
        #copy the geoformat from the src layer
        self.process.dstLayerD[locus][datum][comp].CopyGeoformatFromSrcLayer(self.process.srcLayerD[locus][datum][comp].layer)
        #Create the dst layer

        self.process.dstLayerD[locus][datum][comp].CreateDSWriteRasterArray()
        
        
        '''
        def CreateOpenDstRasterFiles(self, locationid):
        self.dstDateL = []
        for x,dstKey in enumerate(self.layerD):
            for datum in self.layerD[dstKey][locationid]: 

                self.layerD[dstKey][locationid][datum].dstDS = mj_gis.RasterCreateWithFirstLayer(self.layerD[dstKey][locationid][datum].FPN, self.layerD[dstKey][locationid][datum])
                if x == 0:
                    self.dstDateL.append(datum)
        self.dstDateL.sort()
        '''
        
    def CloseDstRasterFiles(self, locationid):
        for dstKey in self.layerD:
            for datum in self.layerD[dstKey][locationid]:           
                self.layerD[dstKey][locationid][datum].dstDS.CloseDS()
        #TGTODO tw051, tw052, two53 etc
        self.dstData.layerD[key][locationid][datum].BAND = dstBAND
        self.dstData.layerD[key][locationid][datum].WriteRasterLayer(flatten = False)             