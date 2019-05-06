import arcpy, glob, os, shutil, math
from arcpy import env
from arcpy.sa import *
import numpy as num


phi=-23.7951;
H= 546;
Pa=285.8;
T=11.8;
RH=55.5
SH=8.35
S=SH/10.68
#dew point temperature calculation
#Td1=243.04*((math.log(RH/100))+((17.625*T)/(243.04+T)))         #numerator
#Td2=17.625-(math.log(RH/100))-((17.625*11.8)/(243.04+11.8))     #denominator
#Td=Td1/Td2;
arcpy.env.overwriteOutput = True
arcpy.env.workspace="D:\\mohana\\post_doc\\papers\\et\\gis\\prism_data1\\tmean\\clip"
dl=arcpy.ListRasters("*","TIF")

t_mean_1="D:\\mohana\\post_doc\\papers\\et\\gis\\prism_data1\\tmean\\clip\\tmean_200001.tif"


dsc=arcpy.Describe(t_mean_1)
sr=dsc.SpatialReference
ext=dsc.Extent
ll=arcpy.Point(ext.XMin,ext.YMin)

arcpy.CheckOutExtension("Spatial")
#net radiation calculation
#v = glob.glob('D:\\mohana\\post_doc\\papers\\et\\gis\\et_calc\\v\\*')
#delta = glob.glob('D:\\mohana\\post_doc\\papers\\et\\gis\\et_calc\\delta\\*tif')
#ws2 = glob.glob('D:/mohana/post_doc/papers/et/gis/narr_data/dswrf/resam_ind/*tif')
output1 = 'D:\\mohana\\post_doc\\papers\\et_regional\\GIS\\TP'
output2 = 'D:\\mohana\\post_doc\\papers\\et_regional\\GIS\\deltaP'
output3 = 'D:\\mohana\\post_doc\\papers\\et_regional\\GIS\\ETP'
output4 = 'D:\\mohana\\post_doc\\papers\\et_regional\\GIS\\ETW'
output5 = 'D:\\mohana\\post_doc\\papers\\et_regional\\GIS\\ET_new'
output6 = 'D:\\mohana\\post_doc\\papers\\et_regional\\GIS\\RTP'
output = 'D:\\mohana\\post_doc\\papers\\et_regional\\GIS\\test'
i=0;

#change negative values
n=len(dl)
for raster in dl:
    print raster
    #break
    #try:
        
    i=i+1
    bn=str(i)
    basename = os.path.basename(raster).split("_")[1]
    #print basename

    dir1='D:\\mohana\\post_doc\\papers\\et_regional\\GIS\\v'
    basename1=dir1+"\\"+"v_"+basename
    v=arcpy.sa.Raster(basename1)
    va=arcpy.RasterToNumPyArray(v)
    #[row,col]=va.shape
    #print va

    dir2='D:\\mohana\\post_doc\\papers\\et_regional\\GIS\\delta'
    basename2=dir2+"\\"+"delta_"+basename
    delta=arcpy.sa.Raster(basename2)
    deltaa=arcpy.RasterToNumPyArray(delta)
    #print deltaa

    dir3='D:\\mohana\\post_doc\\papers\\et_regional\\GIS\\vtc'
    basename3=dir3+"\\"+"vtc_"+basename
    vtc=arcpy.sa.Raster(basename3)
    vtca=arcpy.RasterToNumPyArray(vtc)

    dir4='D:\\mohana\\post_doc\\papers\\et_regional\\GIS\\vd'
    basename4=dir4+"\\"+"vd_"+basename
    vd=arcpy.sa.Raster(basename4)
    vda=arcpy.RasterToNumPyArray(vd)

    dir5='D:\\mohana\\post_doc\\papers\\et_regional\\GIS\\net_rad\\resample\\netrad_p'
    basename5=dir5+"\\"+"netr_"+basename
    rt=arcpy.sa.Raster(basename5)
    rta=arcpy.RasterToNumPyArray(rt)

    dir6='D:\\mohana\\post_doc\\papers\\et_regional\\GIS\\htc'
    basename6=dir6+"\\"+"htc_"+basename
    htc=arcpy.sa.Raster(basename6)
    htca=arcpy.RasterToNumPyArray(htc)

    dir7='D:\\mohana\\post_doc\\papers\\et_regional\\GIS\\alpha'
    basename7=dir7+"\\"+"alpha_"+basename
    alpha=arcpy.sa.Raster(basename7)
    alphaa=arcpy.RasterToNumPyArray(alpha)

    dir8='D:\\mohana\\post_doc\\papers\\et_regional\\GIS\\beta'
    basename8=dir8+"\\"+"beta_"+basename
    beta=arcpy.sa.Raster(basename8)
    betaa=arcpy.RasterToNumPyArray(beta)

    dir9='D:\\mohana\\post_doc\\papers\\et_regional\\GIS\\gammap'
    basename9=dir9+"\\"+"gammap_"+basename
    gammap=arcpy.sa.Raster(basename9)
    gammapa=arcpy.RasterToNumPyArray(gammap)

    dir10='D:\\mohana\\post_doc\\papers\\et_regional\\GIS\\lhv'
    basename10=dir10+"\\"+"lhv_"+basename
    lhv=arcpy.sa.Raster(basename10)
    lhva=arcpy.RasterToNumPyArray(lhv)

    #print lhva

    rastera=arcpy.RasterToNumPyArray(raster)
    [row,col]=rastera.shape
    cellsize=row*col
    
    #print cellsize

    TPp=raster
    vpp=v
    deltapp=delta

    #arcpy.gp.RasterCalculator_sa("Con(\raster\ >= 0, 28.0, (28.0*1.15))", output+"\\"+basename+".tif")



    deltaTP= (((rt/vtc)+vd-vpp+(htc*(Raster(raster)-TPp)))/(deltapp+htc))

    count=0
    value=0
    TP=0
    TH=0.01

    #print deltaTP
    


    TPpa=arcpy.RasterToNumPyArray(TPp)
    vppa=arcpy.RasterToNumPyArray(vpp)
    deltappa=arcpy.RasterToNumPyArray(deltapp)
    deltaTPa=arcpy.RasterToNumPyArray(deltaTP)

    TP2=num.zeros([row,col])
    deltaP2=num.zeros([row,col])

    
    for row1 in range(0,row-1):
        for col1 in range(0,col-1):
            
            TPp1=TPpa.item(row1,col1)
            deltaTP1=deltaTPa.item(row1,col1)
            alpha1=alphaa.item(row1,col1)
            beta1=betaa.item(row1,col1)
            rt1=rta.item(row1,col1)
            vt1=vtca.item(row1,col1)
            vd1=vda.item(row1,col1)
            vpp1=vppa.item(row1,col1)
            htc1=htca.item(row1,col1)
            raster1=rastera.item(row1,col1)
            count=0
            #print TPp1
            while count < count+1:
                TP1=TPp1+deltaTP1
                vP1=6.11*Exp((alpha1*TP1)/(TP1+beta1))
                deltaP1=((alpha1*beta1*vP1)/(TP1+beta1)**2)
                TPp1=TP1
                vpp1=vP1
                deltapp1=deltaP1
                deltaTP1=(((rt1/vt1)+vd1-vpp1+(htc1*(raster1-TPp1)))/(deltapp1+htc1))
                count=count+1
                #print count
                #, raster1, vP1, deltaTP1, TP1, deltaP1, alpha1, beta1

                if count>100:
                    print "didnot converge"
                    break


                if abs(deltaTP1)<=TH:
                    #print count
                    TP2[row1,col1]=TP1
                    deltaP2[row1,col1]=deltaP1
                    break
                
    #newRaster = arcpy.NumPyArrayToRaster(TP2)
         
    TP3 = arcpy.NumPyArrayToRaster(TP2,ll,dsc.meanCellWidth,dsc.meanCellHeight)
    arcpy.DefineProjection_management(TP3, sr)

    deltaP3 = arcpy.NumPyArrayToRaster(deltaP2,ll,dsc.meanCellWidth,dsc.meanCellHeight)
    arcpy.DefineProjection_management(deltaP3, sr)

    
    TP3.save(output1+"\\"+"TP_"+basename)
    deltaP3.save(output2+"\\"+"deltap_"+basename)
            
                
    #potential et        
    ETP = rt-(htc*vtc*(TP3-Raster(raster)))

           
    #net radiation at equil.temp
    RTP = ETP+(gammap*vtc*(TP3-Raster(raster)))

    #estimate wet environment et
    b1="D:\\mohana\\post_doc\\papers\\et\\gis\\et_calc1\\b1\\b1_1.tif"
    b2="D:\\mohana\\post_doc\\papers\\et\\gis\\et_calc1\\b2\\b2_1.tif"

    ETW= b1+(b2/(1+(gammap/deltaP3)))*RTP

   
    #ETW.save(output+"\\"+"ETW_te"+basename)

    cal1=Con((ETW<ETP/2),ETP/2, ETW)
    cal2=Con((cal1>ETP),ETP, cal1)

    ETWa=cal2

    #outputt=output+"\\"+"ETW_te1"+basename
    #arcpy.gp.RasterCalculator_sa("cal2",outputt)
                 

    #estimate the actual et (w/m2)
    ET=(2*ETWa)-ETP

    #latent heat of vapourization(mm/month)
    ETP1=(ETP/lhv)*31
    ETW1=(ETWa/lhv)*31
    ET1=(ET/lhv)*31

    ETP2=Con((ETP1<0),0, ETP1)
    ETW2=Con((ETW1<0),0, ETW1)
    ET2=Con((ET1<0),0, ET1)

    ETP2.save(output3+"\\"+"ETP_"+basename)
    ETW2.save(output4+"\\"+"ETW_"+basename)
    ET2.save(output5+"\\"+"ET_"+basename)
    RTP.save(output6+"\\"+"RTP_"+basename)

    print"completed", basename

        
                        
                   
    #except:
     #   print "some prob"



         


          


        
