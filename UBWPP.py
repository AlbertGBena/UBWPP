##Read binary data from PCL 1300 UHF by DEGREANE
##All thanks the help from Bernard Campistron


import os
import numpy as np
import glob
import struct
import datetime
import datetime as dt
import calendar
import copy
from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel
from netCDF4 import Dataset,num2date,date2num
import math
import sys
from scipy import interpolate as inter_sci

list1=['NPAR-numero de parametros en el encabezado','NWDREC-numero de palabras de 8 bits en el encabezado del conjunto + medidas (encabezado + valores)','NHTS-numero de ventanas de tiempo (alturas) = numero de medidas','NRX-numero de receptor utilizado para la medicion',
       'NPTSP-numero de puntos por espectro','NSPEC-numero de promedios = numero de integraciones de senyal inconsistentes',
       'NCI-numero de integraciones coherentes','IPP-periodo entre pulsos','PW-ancho de pulso (duracion pasada al transmisor)',
       'DLY1-retraso del sistema 1','SPAC1-espaciamiento del sistema 1','NSAM1-numero de medidas en el sistema 2',
       'DLY2-retraso del sistema 2','SPAC2-espaciamiento del sistema 2','NSAM2-numero de medidas en el sistema 2',
       'RECYR-year','RECDAY-day','RECHR-hora','RECMIN-min','RECSEC-sec','RECDAYMON-dia del mes','RECMON-mes',
       'DCOPT-opcion de filtro continuo','WDOPT-ventana de apodizacion','AZ-haz de azimut (0.1grados)','FREQ','BMCODE-numerod e haz',
       'ALT-altitud del sitio en m (MSL)','EL-Elevacion del haz (0.1 grados)','SYS_DELAY1-retardo del sistema para establecer el ancho de banda en nanosegundos',
       'SYS_DELAY2-retardo del sistema para establecer el ancho de banda en nanosegundos','RECERR-error en los datos','CLOCK-periodo de reloj en nanosegundos',
       'SITE','LONGI','LONGI-M','LONGI-S','LAT','LAT-M','LAT-S','CODE_NBR','CODE_NBR_IPP','CODE_NBR_MOM','CODE_NBR_LIGNE','DECODE_TRONC',
       'N_PIC-numero de picos calculados','MOD_ST- VHF/1 ou UHF/0','PUIS_TX-Potencia maxima del transmisor en vatios','VERS_DATA-version en formato de datos',
       'FORMAT-2 compresion,1 tabla calidad','VERS_DSP','ALGO_DSP','RECOUVREMENT-superposicion espectral (en puntos)',
       'PULSE_TX-duracion del pulso (puede ser diferente de la duracion transmitida al transmisor)en nanosegundos','PLUIE-Deteccion de lluvia en la vertical del ciclo anterior.',
       'GAIN_ANT','TX_LOSS_ANT','RX_LOSS_ANT','RECMSSEC-Datacion de los datos ms de la medida.','CONCAT_LIMIT_VENT','CONCAT_LIMIT_SON',
       'NWDREC_HIGH-Parte superior de la palabra numero + medida (32 bits)','NCI_SON-NCI para la parte SON (Si = 0, esto no es un dato RASS)',
       'NCI_HARD-NCI para la tarjeta integradora (Si = 0, cf medida anterior con solo 1 NCI -> Establecer NCI_HARD = NCI)',
       'NPTSP_SON-numero de puntos por espectro','NSPEC_SON-numero de promedios = numero de integraciones de senyal inconsistentes',
       'DCOPT_SON','N_PIC_SON','RECOUVREMENT_SON','FREQ_HIGH','RX_NF','ADC_BITS','ADC_SCALE','ADC_IMPEDANCE','TX_CURRENT','ATTENUATION',
       'ROSE_AURIA','VENT_AURIA','PLUIE_AURIA','TEMP_AURIA','HUMID_AURIA','PRESSION_AURIA','WMO_BLOCK','OWNER_COUNTRY','OWNER_AGENCY',
       'INSTRUMENT (6 es automatica)','ANTENNA_TYPE-4 = matriz planar, 5 = CoCo, 6 = Yaggi, 7 = strip, 14 otros, 15 faltantes','BEAMWIDTH-0.01 GRAUS','STATION_TYPE']


def ContinuityProfile(vector):
    Ones=[];nVect=[]
    for i in range(len(vector)):
        if np.isnan(vector[i]):
            Ones.append(np.nan)
        else:
            if len(nVect)>=2:
                if abs(vector[i]-nVect[-1])<5:
                    nVect.append(vector[i])
                    Ones.append(1)
                else:
                    Ones.append(np.nan)
            else:
                nVect.append(vector[i])
                Ones.append(1)

    Final_vector=np.multiply(vector,Ones)
    return Final_vector,Ones


def CleanProfile(vector2):
    vector=np.copy(vector2)

    Mean=np.nanmean(vector2);Int=2.*np.sqrt(np.nanvar(vector2))
    LimitS=Mean+Int;LimitI=Mean-Int
    for j in range(len(vector2)):
            
        if (vector2[j]>LimitS or vector2[j]<LimitI) and j>1 and j<len(vector2)-2:
            vector[j]=np.nanmean([vector[j+1],vector[j-1]])
    return vector


def AvoidBadSignal(matrix,mask):
    

    for i in range(len(matrix)):#entro en cada instant
        for j in range(len(matrix[i])):#entro en cada altura
            if i==0 or i==len(matrix)-1:

                if i==0 and mask[i]==1:
                    matrix[i][j]=matrix[i+1][j]
                if i==len(matrix[i])-1 and mask[i]==1:
                    matrix[i][j]=matrix[i-1][j]
            else:
                if mask[i]==1:
                    
                    matrix[i][j]=np.nanmean((matrix[i+1][j],matrix[i-1][j]))

                        #if np.isnan(matrix[i+1][j]):
                        #    matrix[i][j]=matrix[i-1][j]+np.nanmean(matrix[i+1][j],matrix[i-1][j])
                        #if np.isnan(matrix[i-1][j]):
                        #    matrix[i][j]=matrix[i+1][j]-np.nanmean(matrix[i+1][j],matrix[i-1][j])
                


    return matrix







def Continuity(vector):
    v1=np.copy(vector)
    for i in tange(len(vector)):
        if i>=3 and i!=len(vector)-1:
            if vector[i-1]-vector[i]>10:#aquí tenim un problema
                v1[i]=np.nan


def Extrapolate(v1,v2,p):#v1 is the height and v2 is the paaremetr to extrapolate the first 2 heights
    v3=np.copy(v2[3:7])
    v4=np.copy(v1[3:7])
    #v1_last=v1[0:7]
    v2_last=np.copy(v2)
    Para=np.count_nonzero(~np.isnan(v3))
    #print(len(v3),Para)
    if Para==len(v3):#condition to extrapolate
    
        #print(v4,v3)
        v1_cut=np.copy(v1[0:7])
        #print(v1_cut)

        f = inter_sci.interp1d(v4,v3,fill_value="extrapolate")
        ynew=f(v1_cut)
        if p==1 or p==2:
            if p==1:
                if abs(ynew[0])>12.:
                    v2_last[0]=np.nan
                else:

                    v2_last[0]=ynew[0]

                if abs(ynew[1])>12.:
                    v2_last[1]=np.nan
                else:

                    v2_last[1]=ynew[1]
                if abs(ynew[2])>12.:
                    v2_last[2]=np.nan
                else:

                    v2_last[2]=ynew[2]
            if p==2:
                if abs(ynew[0])>6.:
                    v2_last[0]=np.nan
                else:

                    v2_last[0]=ynew[0]

                if abs(ynew[1])>6.:
                    v2_last[1]=np.nan
                else:

                    v2_last[1]=ynew[1]
                if abs(ynew[2])>6.:
                    v2_last[2]=np.nan
                else:

                    v2_last[2]=ynew[2]

        else:
            v2_last[0]=ynew[0];v2_last[1]=ynew[1];v2_last[2]=ynew[2]



            
    return v2_last







def PeaksAll(Pot,alt,velX,Noise,att,cte_ra,TimeStamp):
    P_FINAL=[];Z_FINAL=[];C2N=[];VEL=[];VEL_M=[];SIG=[];SK=[];KUR=[];SNR=[];Mask=[];V_AIR=[];V_HIDRO=[]

    for j in range(len(Pot)):
        vector_O=np.copy(Pot[j])
        Noi_plot=np.copy(Noise[j])
        Att=att[j]
        eix_y=alt
        eix_x=velX

        
        v=[];z=[];p_final=[];Sig=[];Sk=[];Kur=[];snr=[];v_m=[];Max=[];Max2=[];v_air=[];v_hidro=[];bimodal=[]#;Min_p=[];c2n=[]
        for k in range(len(vector_O)):
            null_v=np.ones(len(vector_O[k]))*np.nan
            if k==0:#eliminem la primera altura
                snr.append(np.nan)
                v.append(np.nan)
                Sig.append(np.nan)
                Sk.append(np.nan)
                Kur.append(np.nan)
                #c2n.append(np.nan)
                v_m.append(np.nan)
                v_air.append(np.nan)
                v_hidro.append(np.nan)
                #r_w_std_h.append(np.nan)
                Max.append(np.nan)
                z.append(np.nan)
                bimodal.append(np.nan)
                #Min_p.append(0)
            else:
                
                
                #print('resultat',str(vector),str(round(eix_y[k],2)),str(round(Noi_plot[k],2)))
                
                if np.isnan(vector_O[k]).all():
                    snr.append(np.nan)
                    v.append(np.nan)
                    Sig.append(np.nan)
                    Sk.append(np.nan)
                    Kur.append(np.nan)
                    #c2n.append(np.nan)
                    v_m.append(np.nan)
                    v_air.append(np.nan)
                    v_hidro.append(np.nan)
                    #r_w_std_h.append(np.nan)
                    Max.append(np.nan)
                    z.append(np.nan)
                    bimodal.append(np.nan)
                    #Min_p.append(0)
                    
                    vector=vector_O[k]
                else:
                    
                    val_ratio=np.nanmax(vector_O[k])/(np.nanstd(vector_O[k]))
                    val_snr=10.*np.log10(np.nanmax(vector_O[k])/Noi_plot[k])
                    if val_snr<5:# or val_ratio>5:
                        snr.append(np.nan)
                        v.append(np.nan)
                        Sig.append(np.nan)
                        Sk.append(np.nan)
                        Kur.append(np.nan)
                        #c2n.append(np.nan)
                        v_m.append(np.nan)
                        v_air.append(np.nan)
                        v_hidro.append(np.nan)
                        vector=null_v
                        #r_w_std_h.append(Interval)
                        Max.append(np.nan)
                        z.append(np.nan)
                        bimodal.append(np.nan)
                        #Min_p.append(0)
                        
                        
                    else:
                        
                        vector,v_b,vec_air,vec_hidro,moments_m,m_air,m_hidro=peaksfinder(vector_O[k], Noi_plot[k], eix_x)
                        bimodal.append(v_b)
                        if np.nansum(bimodal)>0:#suposem que si hem detectat algun bimodal en el perfil, tenim precipitacio
                            if np.isnan(vec_hidro).all():
                                vec_hidro=np.copy(vec_air)
                        #print(round(np.nanmax(vector),2),round(Noi_plot[k],2),round(alt[k],2))
                        
                        val_ratio=np.nanmax(vector)/(np.nanstd(vector))
                        val_snr=10.*np.log10(np.nanmax(vector)/Noi_plot[k])
                        val_ratio=np.nanmax(vector)/(np.nanstd(vector))
                        val_snr=10.*np.log10(np.nanmax(vector)/Noi_plot[k])
                        ######PER LA COMPOSICIO AIRE+HIDRO
                        
                        W_med=moments_m[0];sig=moments_m[1];ske=moments_m[2];kur=moments_m[3]
                        #####PER EL VALOR DEL AIRE
                        
                        W_med_a=m_air[0];sig_a=m_air[1];ske_a=m_air[2];kur_a=m_air[3]

                        #####PER EL VALOR DEL AIRE
                        
                        W_med_hidro=m_hidro[0];sig_hidro=m_hidro[1];ske_hidro=m_hidro[2];kur_hidro=m_hidro[3]




                        if abs(sig)>=5:# or abs(W_med)>12.:#CONTION IN SPECTRAL WIDTH and vertical velocity
                        
                            snr.append(np.nan)
                            v.append(np.nan)
                            Sig.append(np.nan)
                            Sk.append(np.nan)
                            Kur.append(np.nan)
                            #c2n.append(np.nan)
                            v_m.append(np.nan)
                            v_air.append(np.nan)
                            v_hidro.append(np.nan)
                            vector=null_v
                            #r_w_std_h.append(Interval)
                            Max.append(np.nan)
                            z.append(np.nan)
                            #Min_p.append(0)
                        else:
                            snr.append(val_snr)
                            if np.isnan(vector).all():
                                v.append(np.nan)
                            else:
                                v.append(eix_x[np.nanargmax(vector)])
                            Sig.append(sig)
                            Sk.append(ske)
                            Kur.append(kur)
                            v_m.append(W_med)
                            v_air.append(W_med_a)
                            v_hidro.append(W_med_hidro)
                            #r_w_std_h.append(Interval)
                            Max.append(10*np.log10(np.nansum(vector)))
                            if k==1:
                                z.append(np.nan)#treiem la primera altura ja que dona senyals de potencia molt baixos
                            else:
                                if np.isnan(np.nansum(vector)) or np.nansum(vector)<=0:
                                    z.append(np.nan)
                                else:
                                    z.append(10*np.log10(np.nansum(vector))+20.*np.log10(eix_y[k])+Att-cte_ra)
                            #c2n.append(z[-1]*.1-13.4)
                            #print('minim',round(np.nanmin(vector),1), 'altura ',round(alt[k],1))
                            #if np.nanmin(vector)>100000:
                            #    Min_p.append(1)
                                #print('minim',round(np.nanmin(vector),1),'decision',Min_p[-1], 'altura ',round(alt[k],1))
                            #else:
                            #    Min_p.append(0)

                            
                p_final.append(vector)
        ##DELETE SPORIUS 
        z=Deletesporius(z);v=Deletesporius(v);v_m=Deletesporius(v_m);Sig=Deletesporius(Sig);Sk=Deletesporius(Sk);Kur=Deletesporius(Kur)
        v_m=CleanProfile(v_m)
        v_m,skull=ContinuityProfile(v_m)#analize the continuity on the profile
        Sig=np.multiply(Sig,skull);Sk=np.multiply(Sk,skull);Kur=np.multiply(Kur,skull)
        #c2n=Deletesporius(c2n)
        Min=[]
        for i in range(len(p_final)):
            Min.append(np.nanmin(p_final[i]))
        if j==0:

            if np.nanmin(Min)>100000:#nomes guardo en beam0
                Mask.append(1)

            else:
                Mask.append(0)

        ##EXTRAPOLATION THE VALUES TO the frist three gates--> En Joan ho desestima el 11/03/2022
        
        ####z=Extrapolate(alt,z,0);v=Extrapolate(alt,v,1);v_m=Extrapolate(alt,v_m,1);Sig=Extrapolate(alt,Sig,2);Sk=Extrapolate(alt,Sk,0);Kur=Extrapolate(alt,Kur,0)
        
        c2n=np.copy(z)*.1-13.4#equation from equality......

        #print('despres',np.around(z,decimals=1))
        ##plt.subplot(111)
        ##plt.scatter(v_m,alt)
        ##plt.xlabel('speed (m/s)')
        ##plt.ylabel('height (m)')
        ##plt.grid()
        ##plt.xlim(-12,12)
        ##plt.title('Beam '+str(j)+' at '+str(unix2date(TimeStamp)))
        ##plt.show()

        P_FINAL.append(p_final);Z_FINAL.append(z);C2N.append(c2n);VEL.append(v);VEL_M.append(v_m);SNR.append(snr)
        SIG.append(Sig);SK.append(Sk);KUR.append(Kur);V_AIR.append(v_air);V_HIDRO.append(v_hidro)
        #    P_FINAL=p_final;Z_FINAL=z;C2N=c2n;VEL=v;VEL_M=v_m;SNR=snr
        #    SIG=Sig;SK=Sk;KUR=Kur
        #else:
        #    P_FINAL=np.vstack((P_FINAL, p_final));Z_FINAL=np.vstack((Z_FINAL,z));C2N=np.vstack((C2N, c2n));VEL=np.vstack((VEL, v));VEL_M=np.vstack((VEL_M, v_m))
        #    SNR=np.vstack((SNR, snr));SIG=np.vstack((SIG, Sig));SK=np.vstack((SK, Sk));KUR=np.vstack((KUR, Kur))
    return P_FINAL,Z_FINAL,C2N,VEL,VEL_M,SIG,SK,KUR,SNR,Mask,V_AIR,V_HIDRO#no tenen la primera altura-->H[0]

def Moments(vector,eix_x):
    PT=np.nansum(vector)
    W_med=np.nansum(np.multiply(eix_x,vector))/PT
    sig=np.sqrt(np.nansum(np.multiply((eix_x-W_med)**2,vector))/PT)#spectral width
    ske=np.nansum(np.multiply((eix_x-W_med)**3,vector))/(PT*sig**3)#skewness
    kur=np.nansum(np.multiply((eix_x-W_med)**4,vector))/(PT*sig**4)-3.#kurtosis normalized
    return W_med,sig,ske,kur

def peaksfinder(vector, noise, eixX):#funcio per detectar el pic o pics de la senyal
    v1=np.copy(vector)

    
    v_final=np.ones(len(v1))*np.nan
    v_p1=np.ones(len(v1))*np.nan
    v_p2=np.ones(len(v1))*np.nan

    m_air=[np.nan,np.nan,np.nan,np.nan]
    m_hidro=[np.nan,np.nan,np.nan,np.nan]
    moments_final=[np.nan,np.nan,np.nan,np.nan]
    v_bimodal=0.
    if ~np.isnan(v1).all() and np.nanargmax(v1)>2 and np.nanargmax(v1)<len(v1)-2:
        
        box_kernel = Box1DKernel(4)
        box_kernel2 = Box1DKernel(8)
        v1 = convolve(v1, box_kernel)
        v1 = convolve(v1, box_kernel2)
        v2=np.copy(v1)
        

        if np.nanmax(v1)!=None:

            Ind1=np.nanargmax(v1)
            Max1=v1[Ind1]

            val1=np.nanmax(v1)
            
            
            for i in range(len(v1)-Ind1):
                if i!=0:

                    if v1[Ind1+i]<val1  and (v1[Ind1+i]/Max1)>=.1 and ~np.isnan(v1[Ind1+i]) and i+1<len(v1)-Ind1:

                        val1=v1[Ind1+i]
                        
                        v1[Ind1+i]=np.nan
                    else:
                        Indmax1=Ind1+i
                        break
            
            val1=np.nanmax(v1)
            v2=np.copy(v1)

            for i in range(Ind1):
                if i!=0:

                    if v1[Ind1-i]<val1 and (v1[Ind1-i]/Max1)>=.1 and ~np.isnan(v1[Ind1-i]) and Ind1>i+1:
                        
                        val1=v1[Ind1-i]
                        v1[Ind1-i]=np.nan
                    else:
                        Indmin1=Ind1-i+1
                        break
    
            Inter1=np.arange(Indmin1,Indmax1,1)
    
            for i in range(len(Inter1)):
    
                v_final[Inter1[i]]=vector[Inter1[i]]
                v_p1[Inter1[i]]=vector[Inter1[i]]
            w_1,sig_1,sk_1,ku_1=Moments(v_p1,eixX)
            v1[Ind1]=np.nan


        moments1=[w_1,sig_1,sk_1,ku_1]

        if np.nanmax(v1)!=None and np.nanargmax(v1)>2 and np.nanargmax(v1)<len(v1)-2:
            Ind2=np.nanargmax(v1)
            Max2=v1[Ind2]


            if abs(eixX[Ind1]-eixX[Ind2])>=1 and Max2>=.6*Max1:
                v_bimodal=1.

                val1=np.nanmax(v1)
                
                
                for i in range(len(v1)-Ind2):
                    if i!=0:

                        if v1[Ind2+i]<val1  and (v1[Ind2+i]/Max2)>=.1 and ~np.isnan(v1[Ind2+i]) and i+1<len(v1)-Ind2:
                            val1=v1[Ind2+i]
                            
                            v1[Ind2+i]=np.nan
                        else:
                            Indmax2=Ind2+i
                            break
                val1=np.nanmax(v1)
              
                Ind2=np.nanargmax(v1)
                for i in range(Ind2):
                    if i!=0:

                        if v1[Ind2-i]<val1 and (v1[Ind2-i]/Max2)>=.1 and ~np.isnan(v1[Ind2-i]) and Ind2>i+1:
                            val1=v1[Ind2-i]
                            
                            v1[Ind2-i]=np.nan
                        else:
                            Indmin2=Ind2-i+1
                            break
                Inter2=np.arange(Indmin2,Indmax2,1)
                for i in range(len(Inter2)):
                    v_final[Inter2[i]]=vector[Inter2[i]]
                    v_p2[Inter2[i]]=vector[Inter2[i]]
                w_2,sig_2,sk_2,ku_2=Moments(v_p2,eixX)
                moments2=[w_2,sig_2,sk_2,ku_2]
            else:
                v_binomal=0.
        else:
            v_binomal=0.

    if ~np.isnan(v_p1).all() and ~np.isnan(v_p2).all():
        if abs(eixX[np.nanargmax(v_p2)]-eixX[np.nanargmax(v_p1)])<10:
            if np.nanargmax(v_p1)>np.nanargmax(v_p2):
                v_air=np.copy(v_p2)
                v_hidro=np.copy(v_p1)
                m_air=moments2
                m_hidro=moments1
            else:
                v_air=np.copy(v_p1)
                v_hidro=np.copy(v_p2)
                m_air=moments1
                m_hidro=moments2
        else:
            v_air=np.ones(len(v1))*np.nan
            v_hidro=np.ones(len(v1))*np.nan
            v_final=np.ones(len(v1))*np.nan
            m_air=[np.nan,np.nan,np.nan,np.nan]
            m_hidro=[np.nan,np.nan,np.nan,np.nan]
    else:
        if np.isnan(v_p2).all() and np.isnan(v_p1).all():
            m_hidro=[np.nan,np.nan,np.nan,np.nan]
            m_air=[np.nan,np.nan,np.nan,np.nan]
        else:
            if np.isnan(v_p2).all():
                m_air=moments1
                m_hidro=[np.nan,np.nan,np.nan,np.nan]
        
            
        v_air=np.copy(v_p1)
        v_hidro=np.copy(v_p2)
        
    
    w_f,sig_f,sk_f,ku_f=Moments(v_final,eixX)
    moments_final=[w_f,sig_f,sk_f,ku_f]
##    plt.plot(v_final,'b',alpha=.5)
##    plt.plot(vector,'r',alpha=.5)
##    plt.show()
    

    return v_final,v_bimodal,v_air,v_hidro,moments_final,m_air,m_hidro

def SmoothMatrix(matrix):#the function analyze the matrix and if found differences between three values at the same height with different time genrate the mean of them
    n_matrix=np.copy(matrix)
    LimitValue=5.
    for i in range(len(matrix[0])):#insise 1 hegight
        for j in range(len(matrix)-1):#insed one time
            if j==0:
                #if j==0:#primer perfila  t=0
                Dif1=matrix[j+1][i]-matrix[j][i]
                Dif2=matrix[j+2][i]-matrix[j+1][i]
                if np.isnan(Dif1) and np.isnan(Dif2):#valor entre dos nans
                    n_matrix[j][i]=np.nan
                if ~np.isnan(Dif1):#el valor anterior es anna, pero el poserior no
                    if abs(Dif1)>LimitValue:
                        n_matrix[j][i]=np.nanmean([matrix[j+1][i],matrix[j][i]])
                    
                #else:#darrer perfil a t=len(matrix)-1
                    
                #    Dif1=matrix[j+1][i]-matrix[j][i]
                #    Dif2=matrix[j][i]-matrix[j-1][i]
                #    if np.isnan(Dif1) and np.isnan(Dif2):#valor entre dos nans
                #        n_matrix[j][i]=np.nan
                #    if ~np.isnan(Dif2):#el valor anterior es anna, pero el poserior no
                #        if abs(Dif2)>6.:
                #            n_matrix[j][i]=np.nanmean(matrix[j+2][i],matrix[j+1][i])
                
            else:
                Dif1=matrix[j][i]-matrix[j-1][i]
                Dif2=matrix[j+1][i]-matrix[j][i]
                if np.isnan(Dif1) and np.isnan(Dif2):#valor entre dos nans
                    n_matrix[j][i]=np.nan
                if np.isnan(Dif1) and ~np.isnan(Dif2):#el valor anterior es anna, pero el poserior no
                    if abs(Dif2)>LimitValue:
                        n_matrix[j][i]=np.nanmean([matrix[j+1][i],matrix[j][i]])
                if ~np.isnan(Dif1) and np.isnan(Dif2):#el valor anterior es anna, pero el poserior no
                    if abs(Dif1)>LimitValue:
                        n_matrix[j][i]=np.nanmean([matrix[j-1][i],matrix[j][i]])
                if ~np.isnan(Dif1) and ~np.isnan(Dif2):#el valor anterior es anna, pero el poserior no
                    if abs(Dif1)>LimitValue and abs(Dif2)>LimitValue:
                        n_matrix[j][i]=np.nanmean([matrix[j-1][i],matrix[j+1][i]])

    return n_matrix


                    




            




def ExtractNoise(matrix,vector,Nbeams):#matrix -->la entrada es una matriu, senyal, amb el nombre de beams *num Heights*num bins. El vector es el soroll amb Num beams*NumHeights
    M2=[]#sera la matriu de sortida
    for i in range(Nbeams):#entrem en cada beam
        v1=vector[i]
        m1=matrix[i]
        m2=[]
        for j in range(len(v1)):#entrem en cada altura
            m2.append(m1[j]-v1[j])#restem el valor del soroll
        M2.append(m2)
    return M2

def Inter1D_z(vector):
    Vec=[]
    for i in range(len(vector)):
        if i>0 and i<len(vector)-1:
            if ~np.isnan(vector[i-1]) and np.isnan(vector[i]) and ~np.isnan(vector[i+1]):
                Vec.append(np.nanmean([vector[i-1],vector[i+1]]))
            else:
                Vec.append(vector[i])
        else:
            Vec.append(vector[i])
    return Vec

def Inter1D_Mask(vector):#add value for 1 or 2 nan between real values
    Vec=[];Vec2=[]
    for i in range(len(vector)):
        if i>0 and i<len(vector)-1:
            if ~np.isnan(vector[i-1]) and np.isnan(vector[i]) and ~np.isnan(vector[i+1]):
                Vec.append(1.)
            else:
                Vec.append(vector[i])
        else:
            Vec.append(vector[i])
    for i in range(len(Vec)):
        if i>0 and i<len(Vec)-2:
            if ~np.isnan(Vec[i-1]) and np.isnan(Vec[i]) and np.isnan(Vec[i+1])and ~np.isnan(Vec[i+2]):
                Vec2.append(1.)
            else:
                Vec2.append(Vec[i])
        else:
            Vec2.append(Vec[i])
    return Vec2

def Inter1D_SNR(vector):
    Vec=[]
    LimitSnr=15.
    for i in range(len(vector)):
        if i>0 and i<len(vector)-1:
            if ~np.isnan(vector[i-1]) and vector[i]<=LimitSnr and ~np.isnan(vector[i+1]):
                Vec.append(np.nanmean([vector[i-1],vector[i+1]]))
            else:
                Vec.append(vector[i])
        else:
            Vec.append(vector[i])
    return Vec
def Inter1D_M(vector):
    Vec=[]
    for i in range(len(vector)):
        if i>0 and i<len(vector)-1:
            if ~np.isnan(vector[i-1]) and np.isnan(vector[i]) and ~np.isnan(vector[i+1]):
                Vec.append(1.)
            else:
                Vec.append(vector[i])
        else:
            Vec.append(vector[i])
    return Vec

def Inter1D(vector):
    y=np.asarray(vector)



    indx=np.argwhere(~np.isnan(y))



    if len(indx)>5:
        nou=[];noux=[];Nanx=[]

        for i in range(int(indx[-1]-indx[0])):
            if ~np.isnan(y[indx[0]+i]):
                
            
                nou.append(float(y[indx[0]+i]))
                noux.append(i)
            Nanx.append(i)
              

        y2= np.interp(Nanx, noux, nou)
        Inici=np.ones(int(indx[0]))*np.nan
        Fi=np.ones(len(y)-int(indx[-1]))*np.nan

        y3=np.concatenate((Inici,y2))

        y4=np.concatenate((y3,Fi))
    else:
        y4=np.copy(y)

    return y4



def Ave_NoiseAndPeaks(P,noi,eix_X):#P es potencia, on forma numAv*Nheight*Nbins / Noi es el soroll agupat en numAv*Nheights /eix_X es le vector velocitat té dimensions Nbins
    W_last_F=[];Sigma_F=[];Ske_F=[];Kur_F=[];P_cor_F=[];Snr_F=[]
    #print(len(P),NBeams)
    for k in range(len(P)):
        P_last=np.copy(P[k])
        #Ara el senyal ja no té soroll
        #for i in range(len(noi)):#entrem dins cada numAv
        #    P_withoutNoise=[]
        #P_last=np.copy(P)
        #for j in range(len(noi)):#entrem a cada altura
            #P_withoutNoise.append(P[j]-noi[j])
            
        #    P_last.append(P[j]-noi[j])
        #P_last=np.nanmean(P_C,axis=0)#potencia promig sense soroll
        W_last=[];Sigma=[];Ske=[];Kur=[];P_cor=[];Snr=[]
        for i in range(len(P_last)):#entrema  acda altura
            v=P_last[i]
            P_norm=np.copy(v)/np.nanmax(v)#normalitzada
            std_P=np.nanstd(P_norm)
            P_norm[P_norm<=std_P]=np.nan
            #Trobem el maxims de la senyal
            Av_vel=[];V_vel=[]
            w=[];M_w=[];IndP=[]
            for j in range(len(P_norm)):
                if j==len(P_norm)-2:
                    break
                if j!=0 or j!=len(P_norm)-1:
                    if P_norm[j-1]<=P_norm[j] and P_norm[j+1]<=P_norm[j]:
                        Av_vel.append(eix_X[j])
                        V_vel.append(P_norm[j])
                        IndP.append(j)
                        #print('velo del max',eix_X[j],v_med_C[j],std_med)
            n_peak=[]
            if not Av_vel or np.isnan(v).all():
                if np.isnan(v).all():
                    W=np.nan
                    snr=np.nan
                else:
                    W=0.
                    snr=0.
                #IndPeak=IndP[0]
                Signal=np.ones(len(v))*np.nan
                sig=np.nan
                ske=np.nan
                kur=np.nan
                
            else:
                for j in range(len(V_vel)):
                    if V_vel[j]<=.9 and V_vel[j]>=2.*std_P:
                        w.append(Av_vel[j])
                        M_w.append(V_vel[j])
                        n_peak.append(IndP[j])
                if not w:
                    W=np.nan
                    #IndPeak=IndP[0]
                    Signal=np.ones(len(v))*np.nan
                    sig=np.nan
                    ske=np.nan
                    kur=np.nan
                    snr=np.nan
                else:
                    if len(w)==1:
                        W=w[0]
                        Peak=n_peak[0]
                    else:
                        
                        W=w[np.nanargmax(M_w)]#es la velocitat final obtinguda
                        Peak=n_peak[np.nanargmax(M_w)]
                    ######refem el senyal    
                    Signal=[];i_Ind=0;f_Ind=len(P_norm)
                    for j in range(Peak):
                        if Peak-j-1==0:
                            break
                        if P_norm[Peak-j-1]>P_norm[Peak-j]:
                            i_Ind=Peak-j
                            break
                    for j in range(len(P_norm)-Peak):
                        if (Peak+j+1)==len(P_norm):
                            break
                        if P_norm[Peak+j+1]>P_norm[Peak+j]:
                            
                            f_Ind=Peak+j
                            
                            break
                    #print('index inici i final',i_Ind,f_Ind)
                    for j in range(len(v)):
                        if j>=i_Ind and j<=f_Ind:
                            Signal.append(v[j])
                        else:
                            Signal.append(np.nan)
                    PT=np.nansum(Signal)
                    W_med=np.nansum(np.multiply(eix_X,Signal))/PT
                    sig=np.sqrt(np.nansum(np.multiply((eix_X-W_med)**2,Signal))/PT)#spectral width
                    ske=np.nansum(np.multiply((eix_X-W_med)**3,Signal))/(PT*sig**3)#skewness
                    kur=np.nansum(np.multiply((eix_X-W_med)**4,Signal))/(PT*sig**4)-3.#kurtosis normalized
                    snr=10.*np.log10(PT/noi[k,i])
            W_last.append(W);Sigma.append(sig);Ske.append(ske);Kur.append(kur);P_cor.append(Signal);Snr.append(snr)
        W_last_F.append(W_last);Sigma_F.append(Sigma);Ske_F.append(Ske);Kur_F.append(Kur);P_cor_F.append(P_cor);Snr_F.append(Snr)

    return W_last_F,Sigma_F,Ske_F,Kur_F,P_cor_F,Snr_F# el Snr esta en dB



    

def VertInterp(m1,m2,speed):#in two matrix, out 1 matrix:
    ver=[]
    
    for i in range(len(m2)):#entry in one height gate
        if np.isnan(m2[i]).all():
            ver.append(np.nan)
            
        else:
            ver.append(np.nanmean(speed*(m2[i]/m2[i])))
            
    
    for i in range(len(m2)):#entry in one height gate
        if i>0 and i<len(m2)-2:
            dif1=ver[i]-ver[i-1]
            dif2=ver[i+1]-ver[i]
            dif3=ver[i+1]-ver[i-1]
            if dif3<=6. and dif1>=6.:
##                print('INTERPOLACIO NECESSARIA')
                ver[i]=((ver[i+1]-ver[i-1])/2.)+ver[i-1]
           
    return ver
                
        
        



def Direction (u,v):#determinate the direction of the horizontal wind
    if np.isnan(u) and np.isnan(v):
        Dout=np.nan
    else:
        
        if np.isnan(u) or np.isnan(v):
            if np.isnan(v) and ~np.isnan(u):
                if u>=0.:# wind from west to east
                    Dout=270.
                else:# wind from east to west
                    Dout=90.
            else:
                if v>=0.:# wind from south to north
                    Dout=180.
                else:# wind from south to north
                    Dout=0.
                
        else:
                      
            angle=180.*np.arctan(abs(u/v))/np.pi
            if u>=0 and v>=0:
                Dout=180.+angle
            if u<0 and v<0:
                Dout=angle
            if u>=0 and v<=0:
                Dout=360.-angle
            if u<=0 and v>=0:
                Dout=180.-angle
    return Dout
##type2 is from https://journals.ametsoc.org/view/journals/bams/76/10/1520-0477_1995_076_1717_usmdfn_2_0_co_2.xml  amb alguns canvis
def Type2(ZeVert,W,sigw,eixV,alt,ske):#entra una matriu amb temps i altures
    dv=[]#speed variation with height
    for j in range(len(alt)):

        dv.append(1+3.68*10**-5*alt[j]+1.71*10**-9*alt[j]**2)

    Type_F=[]
    for i in range(len(W)):#entrem en un temps, per tant tenim el perfil d'altura

        EstType=[]#here the values from type are recorded -10 for snow, 0 stratiform rain, 5 drizzle, 10 convective rain, 20 unkown
        Cfact=2#value from cover factor, is the number multiplicate to sigma, Initially I considered as 2.
    
        
        for j in range(len(alt)):#entry in one height gate
            if np.isnan(W[i][j]) or np.isnan(ZeVert[i][j]) or ZeVert[i][j]<=0:
                EstType.append(np.nan)
                
            else:

                if np.power(sigw[i][j],2)<1. and W[i][j]<2.*dv[j]:#cas de neu
                    EstType.append(-10.)
                else:
                    if W[i][j]>=3.*dv[j] :
                        if np.power(sigw[i][j],2)<=4:
                            EstType.append(10.)#AIXO ES STRATIFORME, ara li dono el mateix valor que al covectiu perque els fusiono
                        if W[i][j]>2.*dv[j] and W[i][j]<3.*dv[j]:
                            EstType.append(0)#AIXO ES mixed
            if len(EstType)-1<j:
                if ZeVert[i][j]<20:#limit oer neu
                    EstType.append(-10.)
                else:

                    EstType.append(20.)#AIXO ES UNKNOWN
        
        for m in range(len(EstType)):#to avoid sporadic values
            if m!=0 and m!=len(EstType)-1:
                s1=EstType[m-1]
                s2=EstType[m]
                s3=EstType[m+1]
                if s2==0 and s1==-10 and s3==-10:
                    EstType[m]=-10
                if s2==0 and s1==10 and s3==10:
                    EstType[m]=10
                if s2==20 and s1==10 and s3==10:
                    EstType[m]=10
                if s2==20 and s1==-10 and s3==-10:
                    EstType[m]=-10
                if s2==20 and s1==0 and s3==0:
                    EstType[m]=0
                if s2==20 and np.isnan(s1) and np.isnan(s3):
                    EstType[m]=np.nan
                if s2==10 and np.isnan(s1) and np.isnan(s3):
                    EstType[m]=np.nan
                if s2==0 and np.isnan(s1) and np.isnan(s3):
                    EstType[m]=np.nan
                if s2==-10 and np.isnan(s1) and np.isnan(s3):
                    EstType[m]=np.nan

        if i==0:
            Type_F=EstType
        #    print(EstType)
        else:
            
            Type_F=np.vstack((Type_F,EstType))
    return Type_F






def Type(ZeVert,W,sigw,eixV,alt):#entra una matriu amb temps i altures
    dv=[]#speed variation with height
    for j in range(len(alt)):

        dv.append(1+3.68*10**-5*alt[j]+1.71*10**-9*alt[j]**2)

    Type_F=[]
    for i in range(len(W)):#entrem en un temps, per tant tenim el perfil d'altura

        EstType=[]#here the values from type are recorded -10 for snow, 0 for mixed, and 10 for liquid, 20 unknown
        Cfact=2#value from cover factor, is the number multiplicate to sigma, Initially I considered as 2.
    
        
        for j in range(len(alt)):#entry in one height gate

            
           
            if np.isnan(ZeVert[i][j]) or ZeVert[i][j]<=0:# or raindetc==0.:
                EstType.append(np.nan)
            else:
                z=np.power(10.,ZeVert[i][j]/10.)#convert from dBZ to mm6 m-3
                w=W[i][j]
                s=sigw[i][j]
            

                #stratiform case (M-P)
                vwater=2.65*np.power(z,.114)#values a,b from Atlas et al. 1973
                vsnow=.817*np.power(z,.063)#values a,b from Atlas et al. 1973


                S=(w-(dv[j]*vsnow))
                L=(w-(dv[j]*vwater))



                if abs(S)<abs(L):#snow case

                    if abs(S)<=(Cfact*abs(s)) and abs(L)<=(Cfact*abs(s)):#case not liquid, possible snow

                        if w<=(dv[j]*vwater) and w>=(dv[j]*vsnow):
                            EstType.append(0)#mixed
                        else:
                            EstType.append(20)#unknow


                    else:
                        EstType.append(-10)#snow

                if abs(S)>=abs(L):#rain case
                            
                    if abs(S)<=(Cfact*abs(s)) and abs(L)<=(Cfact*abs(s)):#case liquid
                        

                        if w<=(dv[j]*vwater) and w>=(dv[j]*vsnow):
                            EstType.append(0)#mixed
                        else:
                            EstType.append(20)#unknow
                    else:
                        EstType.append(10)#rain
                if np.isnan(L) and np.isnan(S):
                    EstType.append(np.nan)
                if np.isnan(S) and ~np.isnan(L):#case liquid, but possible wrong election
                    EstType.append(10)#rain
                if ~np.isnan(S) and np.isnan(L):#case not liquid, but possible wrong election
                    EstType.append(-10)#snow
                #print('vels',round(w,2),round(vsnow,2),round(vwater,2),EstType[-1],round(s,2),'height',j)
                #print('S value',round(S,1),'L value',round(L,1))
                #print('height',round(alt[i],0),'rain ',round(dv[i]*vwater,1),'snow ',round(dv[i]*vsnow,1),'real ',round(w,2),'width',round(Cfact*s,2),'Decision',EstType[-1])

        for m in range(len(EstType)):#to avoid sporadic values
            if m!=0 and m!=len(EstType)-1:
                s1=EstType[m-1]
                s2=EstType[m]
                s3=EstType[m+1]
                if s2==0 and s1==-10 and s3==-10:
                    EstType[m]=-10
                if s2==0 and s1==10 and s3==10:
                    EstType[m]=10
                if s2==20 and s1==10 and s3==10:
                    EstType[m]=10
                if s2==20 and s1==-10 and s3==-10:
                    EstType[m]=-10
                if s2==20 and s1==0 and s3==0:
                    EstType[m]=0
                if s2==20 and np.isnan(s1) and np.isnan(s3):
                    EstType[m]=np.nan
                if s2==10 and np.isnan(s1) and np.isnan(s3):
                    EstType[m]=np.nan
                if s2==0 and np.isnan(s1) and np.isnan(s3):
                    EstType[m]=np.nan
                if s2==-10 and np.isnan(s1) and np.isnan(s3):
                    EstType[m]=np.nan
        if i==0:
            Type_F=EstType
        #    print(EstType)
        else:

            Type_F=np.vstack((Type_F,EstType))
        #    print(EstType)                


    return Type_F #return the vector with the estimation, where -10 is snow, 0 is mixed and 10 is liquid, 20 is unknown
            
            
                
                
            
            
            
    
    

def Read1FileNci(NameFile):#read the file to find the minium value for ncivrai and its frequency
    
    f=open(NameFile,'rb')

    byte=f.read() 
    f.close()
    
    LHeader=np.ndarray((1,),dtype='<i2',buffer=byte[0:2])
    Size=int(LHeader[0])
##    print('size',Size)

    Header=np.ndarray((Size,),dtype='<i2',buffer=byte[0:Size*2])
    LenHeader=int(Header[0])
    TotalDades=int(Header[1])/2
    NPTSP=int(Header[4])
    
    NHTS=int(Header[2])
    N_Pic=int(Header[45])
    NSPEC=int(Header[5])#number of cohorent integrations
    NCIVRAI=int(Header[6])*int(Header[63])*int(Header[41])
##    print('NCIVRAI',NCIVRAI)
##    print('intgr incohe',int(Header[5]))
    
    Type='High'
    Contador=0
    NewOffset=0
    TamanyFitxer=0
    Fri=[];NcI=[]
    
    while True:

                
        LenMoments=((4*N_Pic)+1)*NHTS#fiquem 4 perque tenim vel,sig,noise,snr,skew. he tret la quality
        OffsetValues=(LenMoments+Size)*2
        LonTrama=(Size+LenMoments+(NPTSP*NHTS))*2

        
                
        Contador=Contador+1
        TamanyFitxer=TamanyFitxer+(LonTrama/2)


        ##CALCULO ALGUNS PAARAMETRES
           
        pp=int(Header[7])*int(Header[32])*(10**-9)
        fri=1/pp
        Fri.append(fri)
        
        
        ncivrai=int(Header[6])*int(Header[63])*int(Header[41])#is possible that it is wrong
        NcI.append(ncivrai)


        if TamanyFitxer>=len(byte)/2:

            break

        LHeader=np.ndarray((1,),dtype='<i2',buffer=byte[NewOffset+(LonTrama*Contador):NewOffset+(LonTrama*Contador)+2])
        Size=int(LHeader[0])


        Header2=np.ndarray((Size,),dtype='<i2',buffer=byte[NewOffset+(LonTrama*Contador):NewOffset+(LonTrama*Contador)+(Size*2)])
        
    ######Loop to detect the height change
        if int(Header2[2])==int(Header[2]):
            Header=Header2
        else:
            #start the low mode
            NewOffset=LonTrama*Contador
            Contador=0
            Header=Header2
            
        

         
        LenHeader=int(Header[0])
        NPTSP=int(Header[4])
        NHTS=int(Header[2])
        N_Pic=int(Header[45])

    return NCIVRAI,NcI,Fri

def rainrate(al1,al2,al3,vel,ref,U,V):
    #first calculate the value from mu, then el value of ^, then N0 and R
    #constants from Bernard's presentation (pag 82 profilers-DQC_august-2013.pdf)
    

    #information extracted from http://tiwrm.haii.or.th/sharewater_download/books/Radar%20for%20Meteorological%20and%20Atmospheric%20Observations%20(2014).pdf (pag 173)
    c1=50.
    c2=1200.
    c3=3390.
    denwater=10.**6#density from water g m-3
    Mu=[];AA=[];N0=[];R=[];CW=[];HKEF=[];VKEF=[];
    for i in range(len(al1)):
        velH=np.sqrt(np.power(U[i],2)+np.power(V[i],2))
        v=vel[i];a1=al1[i];a2=al2[i];a3=al3[i]

        if np.isnan(ref[i]):
            
            mu=np.nan;A=np.nan;n0=np.nan;r=np.nan;cw=np.nan;hkef=np.nan;vkef=np.nan
        else:
            
            x=np.arange(0,10,.001)#will be the values for mu
            va1=(a1-v)/a2
            va2=1+(a3/(c1*(x**2)+c2*x+c3))
            if va1<=0:
                mu=0.
            else:
                
                f1=np.log(va1)*(-1/(x+7))
    ##            print(f1)

                f2=np.log(va2)
                mu=0.
    ##            plt.plot(f1,'b')
    ##            plt.plot(f2,'r')
    ##            plt.show()
                for j in range(len(f1)-1):
            
                    if f2[j]-f1[j]>=0 and f2[j+1]-f1[j+1]<=0:
                        mu=(x[i]+x[i+1])/2.
##                        print('mu',mu)
                        break
                
            A=(c1*np.power(mu,2))+(c2*mu)+c3# units m-1
                
                
            z=np.power(10.,ref[i]/10.)#convert th dBz to mm6 m-3
            z=float(z)/np.power(10.,18)#convert mm6 m-3 to m3 for calculate n0
            n0=z*np.power(A,mu+7.)/(math.gamma(mu+7.))# units m^-(mu+4)
            r=n0*(math.gamma(mu+4.))*np.pi*(1/6.)*((a1/np.power(A,mu+4.))-(a2/np.power(A+a3,mu+4.)))*(3600.*1000.)#multiply by 3600*1000 to change from m/s to mm/h
            hkef=denwater*(velH**2.)*r/(2.*1000.*3600.)#units g s-3
            vkef=denwater*n0*math.gamma(mu+4.)*(np.pi/12.)*(((a1**3.)/np.power(A,mu+4.))-((a2**3.)/np.power(A+(3.*a3),mu+4.))+((3.*a1*(a2**2.))/np.power(A+(2.*a3),mu+4.))-((3.*a2*(a1**2))/np.power(A+(a3),mu+4.)))#units g s-3
               
            cw=(np.pi*denwater*n0*math.gamma(mu+4.)/(6.*np.power(A,4.+mu)))#g m-3
            
        Mu.append(mu)
        AA.append(A)
        HKEF.append(hkef)
        VKEF.append(vkef)
        if n0<=0 or np.isnan(n0):
            N0.append(np.nan)
        else:
            N0.append(np.log10(n0))
        R.append(r)
        CW.append(cw)
    
    Rout=Deletesporius(R)
    N0out=Deletesporius(N0)
    CWout=Deletesporius(CW)
    HKEFout=Deletesporius(HKEF)
    VKEFout=Deletesporius(VKEF)
##    delete aberrant values
    for j in range(len(Rout)):
        if Rout[j]>100. and ~np.isnan(Rout[j]):
            Rout[j]=np.nan
        if CWout[j]>10. and ~np.isnan(CWout[j]):
            CWout[j]=np.nan
        if HKEFout[j]>100. and ~np.isnan(HKEFout[j]):
            HKEFout[j]=np.nan
        if VKEFout[j]>100. and ~np.isnan(VKEFout[j]):
            VKEFout[j]=np.nan
            

    return Mu,AA,N0out,Rout,CWout,HKEFout,VKEFout
        

def ratioRo(hei,h_0):#To find the rain rate, is neccessary to find the variation of density in function of height

    hei=hei+h_0#hei is in meters and h_0 in meters
    alfa1=[]
    alfa2=[]
    alfa3=[]
    for i in range(len(hei)):
        he=hei[i]

        ratio=(0.0007*(he/1000.)**2)+0.0431*(he/1000.)+0.9921#extract from data Atlas1973 the ratio is (ro0/ro)^0.4 using polinomical regression  
        alfa1.append(9.65*ratio)#units m s-1
        alfa2.append(10.3*ratio)# units m s-1
        alfa3.append(600)#units m-1
    
    #attention alfa1 and alfa2 have units m s-1 and alfa3 has m-1

    return alfa1,alfa2,alfa3
    



def Deletesporius(vector):
    v=[]

    v.append(vector[0])
    for i in range(len(vector)):

        if i>0 and i<(len(vector)-1):

            a=vector[i-1]
            b=vector[i]
            c=vector[i+1]
            if np.isnan(a) and np.isnan(c):
                v.append(np.nan)
            else:
                v.append(b)
    v.append(vector[-1])

    return v
             
def ZeEval(alt,CTE,ATT,P,elva,snr,raindetc):#nomes ho trobo pel valor de beam 1
##    print('ateniation',ATT)
    snr1=snr[0]
    p1=P[0]
    LimitSNR=10.#dB to a good signal, i think that is 15, but try with 5
    LimitReflec=-10.#less of this value the reflective has not sense
    ##The reflectivity for each beam is calculated
    ze1=[];
    for l in range(len(alt)):# IS POSSIBLE TO DO MEAN FROM DBZ? OR IS BETTER TO DO THE MEAN FROM PR??????????
        #print('valor de rain',raindetc)
        if np.isnan(snr1[l]) or raindetc==0 or snr1[l]<LimitSNR:
            z1=np.nan
            Pp1=np.nan
            #print('dins',snr1[l])
        else:
            valuez=np.nansum(p1[l])
            if valuez<=0:
                z1=np.nan
                Pp1=np.nan
            else:
                Pp1=valuez
                z1=((10*np.log10(valuez)-CTE+ATT[0]+20*np.log10(alt[l],out=np.zeros_like(alt[l]), where=(alt[l]!=0))))

                if z1<LimitReflec:
                    z1=np.nan
        
        ze1.append(z1)
    ze1=Deletesporius(ze1)

    return ze1
    
#def WindEval(v1,v2,v3,v4,v5,alt,elva,snr1,snr2,snr3,snr4,snr5):
def WindEval(V,alt,elva,Nbeam):
    
##    print(v1,len(v1));print(v2,len(v2));print(v3,len(v3));print(v4,len(v4));print(v5,len(v5))
    LimitUV=50.#Differences limit for longitudinal and vertical speed (m/s)
    
    
    
    
    if Nbeam==1:

        v1=np.asarray(V[0],float)
        u=[];v=[];#modeule of vertical and horizontal wind
        for i in range(len(alt)):
            u.append(np.nan)
            v.append(np.nan)
    if Nbeam==2:
        v2=np.asarray(V[1],float)
        Covr2=interpolate(v2,alt,elva)
        u=[];v=[];#modeule of vertical and horizontal wind
        for i in range(len(alt)):
            estv=(Covr2[i])/(-1.*np.cos(np.pi*elva/180))
            if estv>LimitUV:
                estv=np.nan
            
            u.append(np.nan)
            v.append(estv)
    if Nbeam==3:
        v2=np.asarray(V[1],float)
        v3=np.asarray(V[2],float)
        Covr2=interpolate(v2,alt,elva)
        Covr3=interpolate(v3,alt,elva)
        u=[];v=[];#modeule of vertical and horizontal wind
        for i in range(len(alt)):
            estv=(Covr3[i]-Covr2[i])/(2.*np.cos(np.pi*elva/180))
            if estv>LimitUV:
                estv=np.nan
            u.append(np.nan)
            v.append(estv)
       
    if Nbeam==4:
        v2=np.asarray(V[1],float)
        v3=np.asarray(V[2],float)
        v4=np.asarray(V[3],float)
        Covr2=interpolate(v2,alt,elva)
        Covr3=interpolate(v3,alt,elva)
        Covr4=interpolate(v4,alt,elva)

        u=[];v=[];#modeule of vertical and horizontal wind
        for i in range(len(alt)):
            estv=(Covr3[i]-Covr2[i])/(2.*np.cos(np.pi*elva/180))
            estu=(Covr4[i])/(-1.*np.cos(np.pi*elva/180))
            if estv>LimitUV:
                estv=np.nan
            if estu>LimitUV:
                estu=np.nan
            
            u.append(estu)
            v.append(estv)
                
        
    if Nbeam==5:
        v2=np.asarray(V[1],float)
        v3=np.asarray(V[2],float)
        v4=np.asarray(V[3],float)
        v5=np.asarray(V[4],float)
        Covr2=interpolate(v2,alt,elva)
        Covr3=interpolate(v3,alt,elva)
        Covr4=interpolate(v4,alt,elva)
        Covr5=interpolate(v5,alt,elva)

        u=[];v=[];#modeule of vertical and horizontal wind
        for i in range(len(alt)):
            estv=(Covr3[i]-Covr2[i])/(2.*np.cos(np.pi*elva/180))
            estu=(Covr5[i]-Covr4[i])/(2.*np.cos(np.pi*elva/180))
            if estv>LimitUV:
                estv=np.nan
            if estu>LimitUV:
                estu=np.nan
            
            u.append(estu)
            v.append(estv)
                
    #m_U=Deletesporius(u)#positive to east, negative to west
    #m_V=Deletesporius(v)#positive north, negative south
    

    
    
    return u,v


def interpolate(Radialvector,alt,angle):
    vector=np.copy(Radialvector)
    Nalt=alt*np.sin(np.pi*angle/180.)
    if np.isnan(vector).all():
##        print('dins',vector)
        vinterp=np.ones(len(vector))*np.nan
    else:
        #vinterp=np.interp(np.arange(len(vector)),np.arange(len(vector))[np.isnan(vector) == False],vector[np.isnan(vector) == False])
        vinterp=np.interp(alt,Nalt,vector)
    

    return vinterp


def ReadSpectrums(Valors,gates,nfft,eix,eixR):#The signal is quantificate with its noise and reesized in the ref axis
    
    Spectr=[]
    dv=eix[5]-eix[4]
    #first create the matrix height x doppler bins
    for i in range(gates):
        s=Valors[nfft*i:nfft*(i+1)]
        
        Spectr.append(s)

    SpectrF=[];Noi=[];NoidB=[]
    for i in range(len(Spectr)):
        vector=Spectr[i]
                
        vectorI=np.asarray(vector)
##        IS NECESSARY DELETE THE VALUES FOR V=0, index 63-64-65, Not it'is necessary
        #vectorCs=np.copy(vectorI)
        #vectorCs[63:66]=np.nan
        #vectorC=np.interp(np.arange(len(vectorCs)),np.arange(len(vectorCs))[np.isnan(vectorCs) == False],vectorCs[np.isnan(vectorCs) == False])#new array where the values 63,64 and 65 are interpolated
        #vectorP=np.copy(vectorC)
        vectorP=np.copy(vectorI)
        
        LenSignal=len(vector)
        MinSignal=[]

        #divideixo la senyal en 32, aixi que cada mostra tindra 4 valors
        for m in range(int(LenSignal/8)):
            WorkingVector=vectorP[m*16:16*(m+1)]
            MinSignal.append(np.nanmean(WorkingVector))
        NoiseW=np.nanmin(MinSignal)
##        print('soroll',NoiseW,'dv',dv)
        Noi.append(NoiseW)
                   
##        print('nois in dB from w',10*np.log10(NoiseW*dv))
        #valNoidB=10*np.log10(NoiseW*dv)#is multiplicated by 10 because is noise from power
        valNoidB=10*np.log10(NoiseW)#is multiplicated by 10 because is noise from power, he tret el dv
        
##        print('soroll dB',valNoidB)
        Snr_pro=10.*np.log(np.nansum(vectorP)/NoiseW)
        if Snr_pro<=5.:#if the signal has higher noise is deleted, if SNR<5 --> sum(signal)>3*Noise
            SpectrF.append(np.ones(len(vectorP))*np.nan)
            NoidB.append(np.nan)
        else:
            NoidB.append(valNoidB)
            #resize the vector from speed vector to speed reference vector
            Ns=np.interp(eixR,eix,vectorP)
            cond=True
            if eix[0]>eixR[0]:
                for o in range(len(eix)):
                    if eix[0]<eixR[o] and cond==True:
    ##                    print('index minim', o)
                        indMin=o
                        cond=False
                    if eix[-1]<eixR[o]:
    ##                    print('index maxim', o)
                        indMax=o
                        break
                carc=np.concatenate([np.ones(indMin)*np.nan,np.ones(indMax-indMin),np.nan*np.ones(len(eix)-indMax)])
                NvectorCorr=Ns*carc
                        
                SpectrF.append(NvectorCorr)
            else:
                SpectrF.append(Ns)
            

    return SpectrF,NoidB,Noi

def convertHalfFloat(OldArray,Nele):#Nele is the number of elements. This fucntion must be optimized (in 524 files the time is 4 minutes)
    OldArray=bytearray(OldArray)
    #print(OldArray)
    array=[hex(x)[2:].zfill(2) for x in OldArray]#convert the array to hex

    test=np.asarray(array)

    Z=[]
    for i in range(int(Nele/2)):

        v1=test[(i*2)+1]

        v2=test[i*2]

        v12=v1+v2
        v3=bin(int(v12,16))[2:].zfill(16)


##    ##Calculem signe
        Signe=1 if v3[0]=='0' else -1

##    ##Calculem expoennt i mantissa
        exp=int(v3[1:9],2)
    
        Reste=int(v3[9:],2)
        R=(int(v3[9:],2))/128.
 
    ##special conditions

        if exp==255 and R!=0:
            z=np.nan
        if  R==0.:
            if exp==255:
                z=np.nan
            if exp==0:
                z=0.
        
        if  exp<255 and exp>0:
            z=Signe*(2**(exp-127))*(1+R)
        if  exp==0 and Reste!=0.:
            z=Signe*(2**(exp-126))*(0+R)


        
        Z.append(z)#create the new array
    return Z


def date2unix(date):
    return calendar.timegm(date.timetuple())
def unix2date(unix):
    return datetime.datetime.utcfromtimestamp(unix)


def Caract(fileid):
    f=open(fileid,'rb')
    byte=f.read()
    f.close()
    TotalBytes=len(byte)
    LHeader=np.ndarray((1,),dtype='<i2',buffer=byte[0:2])
    Size=int(LHeader[0])
    Header=np.ndarray((Size,),dtype='<i2',buffer=byte[0:Size*2])

    Any=int(Header[15])
    Mes=int(Header[21])
    Dia=int(Header[20])

    if Mes<10:
        mes='0'+str(Mes)
    else:
        mes=str(Mes)
    if Dia<10:
        dia='0'+str(Dia)
    else:
        dia=str(Dia)


    FileNameOut=str(Any)+str(mes)+str(dia)



    
    longitude=[int(Header[34]),int(Header[35]),int(Header[36])]
    latitude=[int(Header[37]),int(Header[38]),int(Header[39])]
    totalbeam=int(Header[1])
    lenHe=int(Header[0])
    nptsp=int(Header[4])
    freq=int(Header[25])/10.

    

    BEAM=np.modf(TotalBytes/float(totalbeam))

    off=False

    if BEAM[0]==0:

        numberBeams=int(BEAM[1])
        Gatesnumber=int(Header[2])
        NMODE=1
        #print('number of de beams',numberBeams)
        

    else:

        count=0
        #print('mida de llsita',len(list1))
        while off==False:
            header=np.ndarray((Size,),dtype='<i2',buffer=byte[totalbeam*count:totalbeam*(count+1)*2])
            if int(header[1])!=totalbeam:
                numberBeams=count
                GatesnumberL=int(header[2])
                GatesnumberH=int(Header[2])
                off=True
                NMODE=2
                #print('number of de beams',numberBeams)
                
            else:
                count=count+1

    

    return GatesnumberH,GatesnumberL,numberBeams,longitude,latitude,NMODE,nptsp,freq,FileNameOut

def Read1File(NameFile,freqR,nciR,beams):#read the file
    
    if beams==5:
        Pbeam1=[];Pbeam2=[];Pbeam3=[];Pbeam4=[];Pbeam5=[]
    if beams==4:
        Pbeam1=[];Pbeam2=[];Pbeam3=[];Pbeam4=[]
    if beams==3:
        Pbeam1=[];Pbeam2=[];Pbeam3=[]
    if beams==2:
        Pbeam1=[];Pbeam2=[]
    if beams==1:
        Pbeam1=[]

    f=open(NameFile,'rb')

   

    byte=f.read() 
    f.close()
    
    LHeader=np.ndarray((1,),dtype='<i2',buffer=byte[0:2])
    Size=int(LHeader[0])

    Header=np.ndarray((Size,),dtype='<i2',buffer=byte[0:Size*2])
    LenHeader=int(Header[0])
    TotalDades=int(Header[1])/2
    NPTSP=int(Header[4])
    
    NHTS=int(Header[2])
    N_Pic=int(Header[45])
    NSPEC=int(Header[5])#number of inchorent integrations

##    print('intgr incohe',int(Header[5]))
    
    
    
    Type='High'
    Contador=0
    NewOffset=0
    TamanyFitxer=0
    Atenu=[]
    
    while True:

                
        LenMoments=((4*N_Pic)+1)*NHTS#fiquem 4 perque tenim vel,sig,noise,snr,skew. he tret la quality
        OffsetValues=(LenMoments+Size)*2
        LonTrama=(Size+LenMoments+(NPTSP*NHTS))*2
        #the moments are not necessary toi read, because then we process the signal
        #Moments=convertHalfFloat(byte[NewOffset+(LonTrama*Contador)+(Size*2):NewOffset+(LonTrama*Contador)+2*(Size+LenMoments)],LenMoments*2)
        
        
        Values=convertHalfFloat(byte[NewOffset+(LonTrama*Contador)+OffsetValues:NewOffset+(LonTrama*Contador)+OffsetValues+(2*NPTSP*NHTS)],2*NPTSP*NHTS)
        
                 
        Contador=Contador+1
        TamanyFitxer=TamanyFitxer+(LonTrama/2)


        Att=int(Header[75])#read the attenuation existence
        #print('Is there attenuation?',Att)
        Atenu.append(Att)
        

##        print(Header)
        Any=int(Header[15])
        Mes=int(Header[21])
        Dia=int(Header[20])
        Hora=int(Header[17])
        Min=int(Header[18])
        Sec=int(Header[19])
##        print(Any,Mes,Dia,Hora,Min,Sec)

        dat=datetime.datetime(year=Any,month=Mes,day=Dia,hour=Hora,minute=Min,second=Sec)
######        print('data',dat)
        dat=int(date2unix(dat))

        TimestampTitle=Type+'-'+str(Any)+str(Mes)+str(Dia)+'--'+str(Hora)+'-'+str(Min)+'-'+str(Sec)
        
        ##CALCULO ALGUNS PAARAMETRES
           
        pp=int(Header[7])*int(Header[32])*(10**-9)
        fri=1/pp
        
        ncivrai=int(Header[6])*int(Header[63])*int(Header[41])#is possible that it is wrong
##        ncivrai=int(Header[5])
        freq=int(Header[25])/10.
        #print('frequencia',freq)
##CALCULATE THE REFERENCE VALUES
        friRef=freqR
        ncivraiRef=nciR
        VambRef=(friRef*c)/(4.*ncivraiRef*freq*(10**6))
        dvref=2*VambRef/int(Header[4])
        eixVRef=np.arange(-VambRef,VambRef,dvref)
##CALCULATE THE REAL VALUES
        Vamb=(fri*c)/(4*ncivrai*freq*(10**6))

        dv=2*Vamb/int(Header[4])
##        print('increment vel',dv,ncivrai,freq,c,fri)
        eixV=np.arange(-Vamb,Vamb,dv)
       
        Tprt1=int(Header[9])*int(Header[32])*(10**-9)-int(Header[29])*(10**-9)

        Sprt1=(c/2)*(Tprt1-int(Header[8])*int(Header[32])*(10**-9)*(int(Header[42])-0.5))
        dh=(c/2)*(int(Header[10])*int(Header[32])*(10**-9))
        hmax=Sprt1+(int(Header[2])-1)*dh
        H=np.arange(0,hmax,dh)
####        READ SOME PARAMETERS
        az=int(Header[24])/10.#azimut from beam
        el=int(Header[28])/10.# la elevacio del haz
        pluie=int(Header[54])#it is the same value for all beams, 1 with rain, 0 wihtout rain
        
        


        if el==90 and az==0:
            Pbeam1,NoidB1,Noi1=ReadSpectrums(Values,NHTS,NPTSP,eixV,eixVRef)
            #M1=Moments
            
            
        if el!=90 and az==0:
            Pbeam2,NoidB2,Noi2=ReadSpectrums(Values,NHTS,NPTSP,eixV,eixVRef)
            #M2=Moments
            elev=el

        if el!=90 and az==180:
            Pbeam3,NoidB3,Noi3=ReadSpectrums(Values,NHTS,NPTSP,eixV,eixVRef)
            #M3=Moments
            
                
        if el!=90 and az==90:
##            print('dins 4')
            Pbeam4,NoidB4,Noi4=ReadSpectrums(Values,NHTS,NPTSP,eixV,eixVRef)
            #M4=Moments
##            print('fora 4')
                
        if el!=90 and az==270:
            Pbeam5,NoidB5,Noi5=ReadSpectrums(Values,NHTS,NPTSP,eixV,eixVRef)
            #M5=Moments
 

        ##tornem a llegir el arxiu
        if TamanyFitxer>=len(byte)/2:

            break

        LHeader=np.ndarray((1,),dtype='<i2',buffer=byte[NewOffset+(LonTrama*Contador):NewOffset+(LonTrama*Contador)+2])
        Size=int(LHeader[0])

        Header2=np.ndarray((Size,),dtype='<i2',buffer=byte[NewOffset+(LonTrama*Contador):NewOffset+(LonTrama*Contador)+(Size*2)])
        
    ######Loop to detect the height change
        if int(Header2[2])==int(Header[2]):
            Header=Header2
        else:
            #start the low mode
            NewOffset=LonTrama*Contador
            Contador=0
            if beams==5:
                P1h=Pbeam1;P2h=Pbeam2;P3h=Pbeam3;P4h=Pbeam4;P5h=Pbeam5;
                NoidB1h=NoidB1;NoidB2h=NoidB2;NoidB3h=NoidB3;NoidB4h=NoidB4;NoidB5h=NoidB5;
                Noi1h=Noi1;Noi2h=Noi2;Noi3h=Noi3;Noi4h=Noi4;Noi5h=Noi5
            if beams==4:
                P1h=Pbeam1;P2h=Pbeam2;P3h=Pbeam3;P4h=Pbeam4
                NoidB1h=NoidB1;NoidB2h=NoidB2;NoidB3h=NoidB3;NoidB4h=NoidB4
                Noi1h=Noi1;Noi2h=Noi2;Noi3h=Noi3;Noi4h=Noi4
            if beams==3:
                P1h=Pbeam1;P2h=Pbeam2;P3h=Pbeam3
                NoidB1h=NoidB1;NoidB2h=NoidB2;NoidB3h=NoidB3
                Noi1h=Noi1;Noi2h=Noi2;Noi3h=Noi3
            if beams==2:
                P1h=Pbeam1;P2h=Pbeam2;NoidB1h=NoidB1;NoidB2h=NoidB2;Noi1h=Noi1;Noi2h=Noi2
            if beams==1:
                P1h=Pbeam1;NoidB1h=NoidB1;Noi1h=Noi1

            Header=Header2
            eixVh=eixVRef
            hig=H
            att_h=Atenu
            Atenu=[]
            
        

         
        LenHeader=int(Header[0])
        NPTSP=int(Header[4])
        NHTS=int(Header[2])
        N_Pic=int(Header[45])
                          
    if beams==5:                      
        NoidBh=np.asarray([NoidB1h,NoidB2h,NoidB3h,NoidB4h,NoidB5h]);NoidBl=np.asarray([NoidB1,NoidB2,NoidB3,NoidB4,NoidB5])
        Noih=np.asarray([Noi1h,Noi2h,Noi3h,Noi4h,Noi5h]);Noil=np.asarray([Noi1,Noi2,Noi3,Noi4,Noi5])
        Ph=np.asarray([P1h,P2h,P3h,P4h,P5h]);Pl=np.asarray([Pbeam1,Pbeam2,Pbeam3,Pbeam4,Pbeam5])
    if beams==4:
        NoidBh=np.asarray([NoidB1h,NoidB2h,NoidB3h,NoidB4h]);NoidBl=np.asarray([NoidB1,NoidB2,NoidB3,NoidB4])
        Noih=np.asarray([Noi1h,Noi2h,Noi3h,Noi4h]);Noil=np.asarray([Noi1,Noi2,Noi3,Noi4])
        Ph=np.asarray([P1h,P2h,P3h,P4h]);Pl=np.asarray([Pbeam1,Pbeam2,Pbeam3,Pbeam4])
    if beams==3:
        NoidBh=np.asarray([NoidB1h,NoidB2h,NoidB3h]);NoidBl=np.asarray([NoidB1,NoidB2,NoidB3])
        Noih=np.asarray([Noi1h,Noi2h,Noi3h]);Noil=np.asarray([Noi1,Noi2,Noi3])
        Ph=np.asarray([P1h,P2h,P3h]);Pl=np.asarray([Pbeam1,Pbeam2,Pbeam3])
    if beams==2:
        NoidBh=np.asarray([NoidB1h,NoidB2h]);NoidBl=np.asarray([NoidB1,NoidB2]);Noih=np.asarray([Noi1h,Noi2h]);Noil=np.asarray([Noi1,Noi2])
        Ph=np.asarray([P1h,P2h]);Pl=np.asarray([Pbeam1,Pbeam2])
    if beams==1:
        NoidBh=np.asarray([NoidB1h]);NoidBl=np.asarray([NoidB1]);Noih=np.asarray([Noi1h]);Noil=np.asarray([Noi1])
        Ph=np.asarray([P1h]);Pl=np.asarray([Pbeam1])
    
##    print('leng',len(P1h),len(P2h),len(P3h),len(P4h),len(P5h),len(Pbeam1),len(Pbeam2),len(Pbeam3),len(Pbeam4),len(Pbeam5))
    return Ph,Pl,eixVRef,hig,elev,H,att_h,Atenu,dat,NSPEC,pluie,NoidBh,NoidBl,Noih,Noil#,pluie1h,pluie2h,pluie3h,pluie4h,pluie5h,pluie1,pluie2,pluie3,pluie4,pluie5#dat is the time mark for the last beam

def peaksfinder2(vector):#funcio per detectar el pic o pics de la senyal
    v1=np.copy(vector)
    std1=np.nanstd(v1)
    v_final=np.ones(len(v1))*np.nan
    #v1[v1<np.nanmean(v1)+0.*std1]=np.nan
    #trobo la posició dels dos màxims
    if ~np.isnan(v1).all():
        Ind1=np.nanargmax(v1)
        val1=np.nanmax(v1)
        for i in range(len(v1)-Ind1):
            if v1[Ind1+i]<=val1:
                val1=v1[Ind1+i]
                v1[Ind1+i]=np.nan
            else:
                Indmax1=Ind1+i
                break
        val1=np.nanmax(v1)
        for i in range(Ind1):
            if v1[Ind1-i]<=val1:
                val1=v1[Ind1-i]
                v1[Ind1-i]=np.nan
            else:
                Indmin1=Ind1-i
                break
    Inter1=np.arange(Indmin1,Indmax1,1)
    for i in range(len(Inter1)):
        v_final[Inter1[i]]=vector[Inter1[i]]
    if ~np.isnan(v1).all():
        Ind2=np.nanargmax(v1)
        val2=np.nanmax(v1)
        for i in range(len(v1)-Ind2):
            if v1[Ind1+i]<=val1:
                val2=v1[Ind1+i]
                v1[Ind1+i]=np.nan
            else:
                Indmax2=Ind1+i
                break
        val1=np.nanmax(v1)
        for i in range(Ind2):
            if v1[Ind1-i]<=val2:
                val2=v1[Ind1-i]
                v1[Ind1-i]=np.nan
            else:
                Indmin2=Ind2-i
                break
    Inter2=np.arange(Indmin2,Indmax2,1)
    
    for i in range(len(Inter2)):
        v_final[Inter2[i]]=vector[Inter2[i]]

    return v_final


def group_consecutives(vals, step=1):
    """Return list of consecutive lists of numbers from vals (number list)."""
    run = []
    result = [run]
    expect = None
    for v in vals:
        if (v == expect) or (expect is None):
            run.append(v)
        else:
            run = [v]
            result.append(run)
        expect = v + step
    return result


def ElimDupli(vector):
    
    
    b=abs(np.round(vector,1))

    for i in range(len(b)):
        
        for j in range(len(b)):
            if i!=j and b[i]==b[j]:
                b[i]=np.nan
                b[j]=np.nan
            
    b=np.asarray(b)
    divb=np.where(b==0,1,b)       
    bnorm=b/divb
    #print(b)
    resul=bnorm*vector
    return resul 

def CleanMatrix(matrix):
    new_matrix=np.copy(matrix)
    for i in range(len(matrix)):
        if i>1 and i<len(matrix)-2:
            for j in range(len(matrix[i])):
                if np.isnan(matrix[i-1][j]) and ~np.isnan(matrix[i][j]) and np.isnan(matrix[i+1][j]):
                    new_matrix[i][j]=np.nan
    return new_matrix



####################CONSTANTS#############
c=299792458.
cte_high=116.5
cte_low=111.3
np.warnings.filterwarnings('ignore')#to avoid the error messages

########INCLUDE THE OPTIONS IN EXECUTATION
if len(sys.argv)==1:
    option=0
c1_op=0;c2_op=0;c3_op=0;h0_opt=0.#c_opt=0;c1=0#;c2=0;c3=0;h0_opt=np.nan
if len(sys.argv)>1:
    for i in sys.argv:
   
        if i=='-Pot':
            print('\nYour chosen option is to save the raw potency from each beam in two mode\n')
            c1_op=1
        if i=='-NameFile':
            c2_op=1
        if i[0:2]=='-h':
            #print('The first height has been changed\n')
            h0_opt=float(i[2:])
            c3_op=1
        if i[0:3]=='-cH':#calibration constant High mode
            cte_high=float(i[3:])
        if i[0:3]=='-cL':#calibration constant Low mode
            cte_low=float(i[3:])
            

Count_t0=dt.datetime.now()
tt00=dt.datetime.now()
##print(dt.datetime.now())

print('Insert the path where the dat files are --for instance d:/UHF/')
Root=input()  #input from the user 

if not os.path.isdir(Root):
    print('\nthe folder is not exist, please check your syntax')
    exit()
os.chdir(Root)
folder=Root    
dircf=glob.glob(Root+'*.dat')
if len(dircf)==0:
    print('In this folder there are not dat files')
    exit()
dircf=np.sort(dircf)
print('\nIn this folder there are '+str(len(dircf))+' dat files')
print('\nThe script generates netcdf file.')


d=glob.glob(Root+'*.dat')
NumberFiles=len(d)


NumberJoin=2#average for three samples the value is 2

gatesH,gatesL,NBeams,lonG,laT,Nmode,nptsp,Freq,PossNameFile=Caract(d[0])
if c2_op==1:
    print('Please write the netcdf name (without extension.nc)\n')
    filenameplot=input()
else:
    filenameplot=PossNameFile

count=0

print('\nUHF LOCATION\nlongitude [degrees,minutes,seconds] ',lonG,'\nlatitude  [degrees,minutes,seconds] ',laT)
if Nmode==2:

    print('\nThere are two modes, high and low')
else:
    print('\nThere is one mode, high or low ')
print('Radar Frequency ',Freq,' MHz and the wavelenght is ',round(c/(1000.*Freq),1), ' mm' )
print('Number of beams ',NBeams)


#READ ALL FILES TO FIND THE LONGEUR IMPULSION CORRESPON AT MINIMUM NCIVRAI

Ncivrai=[];NCItota=[];Fri=[];NciTotal=[]
for name in d:
    

    ncivrai,NCI,ffri=Read1FileNci(name)

    NciTotal=np.concatenate([NciTotal,NCI])
    Fri=np.concatenate([Fri,ffri])

FreqRef=Fri[np.argmin(NciTotal)]
NCIref=min(NciTotal)
##print('frequency repetition impulse',FreqRef,' ncivrai', NCIref)
#print('read all files fi',dt.datetime.now())
##TiF=TimeIntervalFile(d[0],FreqRef,NCIref)



####STARTING THE NETCDF FILE

folder=Root

##################CREATE THE NETDF############
dataset=Dataset(folder+filenameplot+'.nc','w',format='NETCDF4')
dataset.description='Processed data from UHF radar'
dataset.author='Albert Garcia Benad'+u'\xed'
dataset.orcid='0000-0002-5560-4392 '
dataset.acknowledgement='Thanks a lot for the support, patience and goods advises of Dr. Bernard Campistron and the Dr. Philipp Currier.'    
dataset.location_Longitude=str(lonG[0])+'degrees,'+str(lonG[1])+'minutes,'+str(lonG[2])+'seconds'
dataset.location_Latitude=str(laT[0])+'degrees,'+str(laT[1])+'minutes,'+str(laT[2])+'seconds.'
dataset.Number_Of_Beams=str(NBeams)
dataset.Number_Of_modes=str(Nmode)

####creates the menu
if c1_op==1:
    Origin=dataset.createGroup("Signal from file without processing")

High=dataset.createGroup("High mode")
Low=dataset.createGroup("Low mode")


##Create the dimensions
dataset.createDimension('time',None)#time hasn't dimension!!! It is very important
####dataset.createDimension('timeM',None)#time hasn't dimension!!! It is very important
####dataset.createDimension('time_utc_M',None)#to mean
dataset.createDimension('time_utc',None)#instantaneous
dataset.createDimension('Hi_height',gatesH)
dataset.createDimension('Low_height',gatesL)
dataset.createDimension('Number_Bins',nptsp)
##dataset.createDimension('Rain',None)#time hasn't dimension!!! It is very important

##Creates the variables
nc_times=dataset.createVariable('Time','float64',('time',))
####nc_timesMean=dataset.createVariable('TimeM','float64',('timeM',))
nc_ranges_H=dataset.createVariable('Hi_height','f',('Hi_height',))
nc_ranges_L=dataset.createVariable('Low_height','f',('Low_height',))
nc_ranges_Nbins=dataset.createVariable('Number_Bins','f',('Number_Bins',))
nc_Format_times=dataset.createVariable('time_utc', 'float64', ('time_utc',))
####nc_Format_timesM=dataset.createVariable('time_utc_M', 'float64', ('time_utc_M',))

##nc_rain=dataset.createVariable('rain','float64',('Rain',))

##Create the units
nc_times.units = 'UNIX Time Stamp, SECOND SINCE 1970-01-01'
nc_times.description='Time in unix format'
####nc_timesMean.units = 'UNIX Time Stamp, SECOND SINCE 1970-01-01'
####nc_timesMean.description='Time in unix format'

nc_Format_times.units='seconds since 1970-01-01'
nc_Format_times.calendar='standard'
nc_Format_times.decription='time UTC'

####nc_Format_timesM.units='seconds since 1970-01-01'
####nc_Format_timesM.calendar='standard'
####nc_Format_timesM.decription='time UTC'
nc_ranges_H.units = 'm'
nc_ranges_L.units = 'm'
if c3_op==0:
    nc_ranges_H.description = 'Heights in high mode above ground level'
    nc_ranges_L.description = 'Heights in low mode above ground level'
else:
    nc_ranges_H.description = 'Heights in high mode above sea level'
    nc_ranges_L.description = 'Heights in low mode above sea level'


nc_ranges_Nbins.units='none'
nc_ranges_Nbins.description='Number of Doppler bins'


##nc_rain.description='This paremeter is the Pluie from UHF, 0 is not precipiatation and 1 is precipitation'
##nc_rain.units = 'None'



count=0#counter for create dimension in time variable from sinal in netcdf
countM=0#counter for create dimension in time variable from average signal in netcdf
##M1H=[];M2H=[];M3H=[];M4H=[];M5H=[];M1L=[];M2L=[];M3L=[];M4L=[];M5L=[];
PLUIE=[]
##Noi1h=[];Noi2h=[];Noi3h=[];Noi4h=[];Noi5h=[];Noi1l=[];Noi2l=[];Noi3l=[];Noi4l=[];Noi5l=[];
MH=[];ML=[];NoiH=[];NoiL=[];Att_h=[];Att_l=[]
countwork=0;Filesdone=0
##print 'Work in progress(0%%)',
dat_0=int(date2unix(dt.datetime.now()))
for name in d:
    Filesdone+=1
    
    
    

    #print('Read 1 file ini',dt.datetime.now())
    MaxClose=3#valor respecte la std,Bernard aconsella 3#ATEMCIO VALOR MOLT IMPORTANT
    Count_t1=dt.datetime.now()
##    print('arxiu analitzat',name)
    PH,PL,velX,H,eleva,h,att_h,att_l,TimeStamp,NUMINT,Pluie,noidBh,noidBl,noih,noil=Read1File(name,FreqRef,NCIref,NBeams)
##    if c3_op==1:
##        H=H+h0_opt
##        h=h+h0_opt


    #print('parametre pluja',Pluie,unix2date(TimeStamp))
    

    p_h,z_h,c2n_h,vel_h,velM_h,sig_h,sk_h,kur_h,snr_h,mask_h,v_air_h,v_hidro_h=PeaksAll(PH,H,velX,noih,att_h,cte_high,TimeStamp)#llegim i treiem els paràmetres
    

    p_l,z_l,c2n_l,vel_l,velM_l,sig_l,sk_l,kur_l,snr_l,mask_l,v_air_l,v_hidro_l=PeaksAll(PL,h,velX,noil,att_l,cte_low,TimeStamp)#llegim i treiem els paràmetres
    
    U_h,V_h=WindEval(velM_h,H,eleva,NBeams)
    
    alfa1_h,alfa2_h,alfa3_h=ratioRo(H,h0_opt)
    mu_h,A_h,N0_h,R_h,Cw_h,Hkef_h,Vkef_h=rainrate(alfa1_h,alfa2_h,alfa3_h,velM_h[0],z_h[0],U_h,V_h)

    U_l,V_l=WindEval(velM_l,h,eleva,NBeams)
    
    alfa1_l,alfa2_l,alfa3_l=ratioRo(h,h0_opt)
    mu_l,A_l,N0_l,R_l,Cw_l,Hkef_l,Vkef_l=rainrate(alfa1_l,alfa2_l,alfa3_l,velM_l[0],z_l[0],U_l,V_l)

    
    if count==0:
        Z0_h=z_h[0];W_h=velM_h[0];Sig_h=sig_h[0];Sk_h=sk_h[0];Kur_h=kur_h[0];C2N_h=c2n_h[0]
        SnR_h=snr_h[0];Rr_h=R_h
        u_h=U_h;v_h=V_h;Mu_h=mu_h;a_h=A_h;n0_h=N0_h;CW_h=Cw_h
        HKEF_h=Hkef_h;VKEF_h=Vkef_h
        MASK_h=mask_h
##        W_air_h=v_air_h[0]
##        W_hidro_h=v_hidro_h[0]

        Z0_l=z_l[0];W_l=velM_l[0];Sig_l=sig_l[0];Sk_l=sk_l[0];Kur_l=kur_l[0];C2N_l=c2n_l[0]
        SnR_l=snr_l[0];Rr_l=R_l
        u_l=U_l;v_l=V_l;Mu_l=mu_l;a_l=A_l;n0_l=N0_l;CW_l=Cw_l
        HKEF_l=Hkef_l;VKEF_l=Vkef_l
        MASK_l=mask_l
##        W_air_l=v_air_l[0]
##        W_hidro_l=v_hidro_l[0]

    else:
        Z0_h=np.vstack((Z0_h,z_h[0]));W_h=np.vstack((W_h,velM_h[0]));Sig_h=np.vstack((Sig_h,sig_h[0]));Sk_h=np.vstack((Sk_h, sk_h[0]));Kur_h=np.vstack((Kur_h,kur_h[0]));C2N_h=np.vstack((C2N_h,c2n_h[0]))
        SnR_h=np.vstack((SnR_h,snr_h[0]));Rr_h=np.vstack((Rr_h,R_h))
        u_h=np.vstack((u_h,U_h));v_h=np.vstack((v_h,V_h));Mu_h=np.vstack((Mu_h,mu_h));a_h=np.vstack((a_h, A_h));n0_h=np.vstack((n0_h, N0_h));CW_h=np.vstack((CW_h, Cw_h))
        HKEF_h=np.vstack((HKEF_h,Hkef_h));VKEF_h=np.vstack((VKEF_h,Vkef_h))
        MASK_h=np.concatenate((MASK_h,mask_h))

##        W_air_h=np.vstack((W_air_h,v_air_h[0]));W_hidro_h=np.vstack((W_hidro_h,v_hidro_h[0]))

        Z0_l=np.vstack((Z0_l,z_l[0]));W_l=np.vstack((W_l,velM_l[0]));Sig_l=np.vstack((Sig_l,sig_l[0]));Sk_l=np.vstack((Sk_l, sk_l[0]));Kur_l=np.vstack((Kur_l,kur_l[0]));C2N_l=np.vstack((C2N_l,c2n_l[0]))
        SnR_l=np.vstack((SnR_l,snr_l[0]));Rr_l=np.vstack((Rr_l,R_l));
        u_l=np.vstack((u_l,U_l));v_l=np.vstack((v_l,V_l));Mu_l=np.vstack((Mu_l,mu_l));a_l=np.vstack((a_l, A_l));n0_l=np.vstack((n0_l, N0_l));CW_l=np.vstack((CW_l, Cw_l))
        HKEF_l=np.vstack((HKEF_l,Hkef_l));VKEF_l=np.vstack((VKEF_l,Vkef_l))
        MASK_l=np.concatenate((MASK_l,mask_l))
##        W_air_l=np.vstack((W_air_l,v_air_l[0]));W_hidro_l=np.vstack((W_hidro_l,v_hidro_l[0]))
        


    PH_c=ExtractNoise(PH,noih,NBeams)
    PL_c=ExtractNoise(PL,noil,NBeams)


        
    Count_t2=dt.datetime.now()
    
    

    ##adding valu in time axes

    nc_times[count:count+1]=TimeStamp
    nc_Format_times[count:count+1]=date2num(unix2date(TimeStamp),units=nc_Format_times.units,calendar=nc_Format_times.calendar)
    ##nc_rain[count:count+1]=Pluie
    

##nc_timesUV[:]=Tco
    if c3_op==1:
        nc_ranges_H[:]=np.array(H,dtype='f4')+h0_opt
        nc_ranges_L[:]=np.array(h,dtype='f4')+h0_opt
    else:
        nc_ranges_H[:]=np.array(H,dtype='f4')
        nc_ranges_L[:]=np.array(h,dtype='f4')
    
    nc_ranges_Nbins[:]=np.arange(0,nptsp,1)
##    print(NPTSP)
##    print(np.arange(0,128,1),len(np.arange(0,128,1)))

    

    
    ncShape2DH = ('time_utc','Hi_height',)#per muntar les altures
    ncShape2DL = ('time_utc','Low_height',)
####    ncShape2DMeanH = ('time_utc_M','Hi_height',)#per muntar les altures
####    ncShape2DMeanL = ('time_utc_M','Low_height',)
    ncShape3D=('time_utc','Hi_height','Number_Bins',)
    ncShape3Dlow=('time_utc','Low_height','Number_Bins',)
    
    if count==0:
######        We declare the variables
##        in high mode

        


        ###COMBINED SPEED
        nc_V1_h=High.createVariable('W_high','f',ncShape2DH)
        nc_V1_h.description='Vertical air and hidrometeor speed, where negative is upward movement, in high mode'
        nc_V1_h.units='m/s'

        nc_Sig1_h=High.createVariable('Spectral width_high','f',ncShape2DH)
        nc_Sig1_h.description='Spectral width from the Vertical air and hidrometeor speed'
        nc_Sig1_h.units='m/s'

        nc_Sk1_h=High.createVariable('Skewness_high','f',ncShape2DH)
        nc_Sk1_h.description='Skewness from the Vertical air and hidrometeor speed'
        nc_Sk1_h.units='none'

        nc_Kur1_h=High.createVariable('Kurtosis_high','f',ncShape2DH)
        nc_Kur1_h.description='Kurtosis from the Vertical air and hidrometeor speed'
        nc_Kur1_h.units='none'

        ####FROM AIR

##        nc_V1_air_h=High.createVariable('W_air_high','f',ncShape2DH)
##        nc_V1_air_h.description='Vertical air speed, where negative is upward movement, in high mode'
##        nc_V1_air_h.units='m/s'
##
##        ####FROM HIDRO
##
##        nc_V1_hidro_h=High.createVariable('W_hidro_high','f',ncShape2DH)
##        nc_V1_hidro_h.description='Vertical hidrometeor speed, where negative is upward movement, in high mode'
##        nc_V1_hidro_h.units='m/s'

        

        nc_Z_h=High.createVariable('Z_high','f',ncShape2DH)
        nc_Z_h.description='Reflectivity'
        nc_Z_h.units='dBZ'

        nc_U_h=High.createVariable('U_high','f',ncShape2DH)
        nc_U_h.description='Eastward wind or Zonal wind, data from beam 4(to east) and 5 (to west)'
        nc_U_h.units='m/s'

        nc_V_h=High.createVariable('V_high','f',ncShape2DH)
        nc_V_h.description='Northward wind or Meridional wind, data from beam 2(to north) and 3 (to south)'
        nc_V_h.units='m/s'


        nc_Snr1_h=High.createVariable('Snr_high','f',ncShape2DH)
        nc_Snr1_h.description='snr from beam 1 in high mode'
        nc_Snr1_h.units='dB'

        nc_Rr_h=High.createVariable('RR_hm','f4',ncShape2DH)
        nc_Rr_h.description='Rain rate in high mode for each height'
        nc_Rr_h.units='mm h-1'

        nc_EType_h=High.createVariable('Type_hm','f4',ncShape2DH)
        nc_EType_h.description='Estimation from hydrometeor -10 snow, 0 mixed, 10 rain, and 20 unknown in hm'
        nc_EType_h.units='none'

        nc_EType2_h=High.createVariable('Type_2_hm','f4',ncShape2DH)
        #nc_EType2_h.description='EStimation from hydrometeor --10 for snow, 0 stratiform rain, and 10 convective rain, 20 unkown'
        nc_EType2_h.description='Estimation from hydrometeor from Ralph --10 for snow, 0 mixed, and 10 convective/stratiform rain, 20 unkown'
        nc_EType2_h.units='none'

        nc_Hkef_h=High.createVariable('HKEF_hm','f4',ncShape2DH)
        nc_Hkef_h.description='Horizontal kinetic energy flux in high mode for each height'
        nc_Hkef_h.units='g s-3'

        nc_Vkef_h=High.createVariable('VKEF_hm','f4',ncShape2DH)
        nc_Vkef_h.description='Vertical kinetic energy flux in high mode for each height'
        nc_Vkef_h.units='g s-3'

        nc_No_h=High.createVariable('log10(N0)_hm','f4',ncShape2DH)
        nc_No_h.description='log10(N0) in high mode for each height'
        nc_No_h.units='log10 (m-4-'+u'\u03BC'+' )'

        nc_A_h=High.createVariable('slope_parameter_hm','f',ncShape2DH)
        nc_A_h.description='Slope parameter in high mode for each height'
        nc_A_h.units='m-1'

        nc_Mu_h=High.createVariable('shape_parameter_hm','f',ncShape2DH)
        nc_Mu_h.description='shape parameter in high mode for each height'
        nc_Mu_h.units='None'

        nc_Cw_h=High.createVariable('Liquid_water_content_hm','f4',ncShape2DH)
        nc_Cw_h.description='Amount Liquid Water in high mode for each height'
        nc_Cw_h.units='g m-3'

        nc_C2n_h=High.createVariable('C2n_hm','f4',ncShape2DH)
        nc_C2n_h.description=' refractive index for high mode in log'
        nc_C2n_h.units='log(m(-2/3))'

####HO TREC PER ANAR MES RÀPID EN LES PROBES
        if NBeams==1 and c1_op==1:
            nc_P0h=Origin.createVariable('P_0_hm','f',ncShape3D)
            nc_P0h.description='Potency without noise from beam 0, vertical, in high mode'
            nc_P0h.units='watt'

            nc_P0l=Origin.createVariable('P_0_lm','f',ncShape3Dlow)
            nc_P0l.description='Potency without noise from beam 0, vertical, in low mode'
            nc_P0l.units='watt'
        if NBeams==2 and c1_op==1:
            nc_P0h=Origin.createVariable('P_0_hm','f',ncShape3D)
            nc_P0h.description='Potency without noise from beam 0, vertical, in high mode'
            nc_P0h.units='watt'
            nc_P1h=Origin.createVariable('P_1_hm','f',ncShape3D)
            nc_P1h.description='Potency without noise from beam 1 in high mode'
            nc_P1h.units='watt'

            nc_P0l=Origin.createVariable('P_0_lm','f',ncShape3Dlow)
            nc_P0l.description='Potency without noise from beam 0, vertical, in low mode'
            nc_P0l.units='watt'
            nc_P1l=Origin.createVariable('P_1_lm','f',ncShape3Dlow)
            nc_P1l.description='Potency without noise from beam 1 in low mode'
            nc_P1l.units='watt'
        if NBeams==3 and c1_op==1:
            nc_P0h=Origin.createVariable('P_0_hm','f',ncShape3D)
            nc_P0h.description='Potency without noise from beam 0, vertical, in high mode'
            nc_P0h.units='watt'
            nc_P1h=Origin.createVariable('P_1_hm','f',ncShape3D)
            nc_P1h.description='Potency without noise from beam 1 in high mode'
            nc_P1h.units='watt'
            nc_P2h=Origin.createVariable('P_2_hm','f',ncShape3D)
            nc_P2h.description='Potency without noise from beam 2 in high mode'
            nc_P2h.units='watt'

            nc_P0l=Origin.createVariable('P_0_lm','f',ncShape3Dlow)
            nc_P0l.description='Potency without noise from beam 0, vertical, in low mode'
            nc_P0l.units='watt'
            nc_P1l=Origin.createVariable('P_1_lm','f',ncShape3Dlow)
            nc_P1l.description='Potency without noise from beam 1 in low mode'
            nc_P1l.units='watt'
            nc_P2l=Origin.createVariable('P_2_lm','f',ncShape3Dlow)
            nc_P2l.description='Potency without noise from beam 2 in low mode'
            nc_P2l.units='watt'
        if NBeams==4 and c1_op==1:
            nc_P0h=Origin.createVariable('P_0_hm','f',ncShape3D)
            nc_P0h.description='Potency without noise from beam 0, vertical, in high mode'
            nc_P0h.units='watt'
            nc_P1h=Origin.createVariable('P_1_hm','f',ncShape3D)
            nc_P1h.description='Potency without noise from beam 1 in high mode'
            nc_P1h.units='watt'
            nc_P2h=Origin.createVariable('P_2_hm','f',ncShape3D)
            nc_P2h.description='Potency without noise from beam 2 in high mode'
            nc_P2h.units='watt'
            nc_P3h=Origin.createVariable('P_3_hm','f',ncShape3D)
            nc_P3h.description='Potency without noise from beam 3 in high mode'
            nc_P3h.units='watt'

            nc_P0l=Origin.createVariable('P_0_lm','f',ncShape3Dlow)
            nc_P0l.description='Potency without noise from beam 0, vertical, in low mode'
            nc_P0l.units='watt'
            nc_P1l=Origin.createVariable('P_1_lm','f',ncShape3Dlow)
            nc_P1l.description='Potency without noise from beam 1 in low mode'
            nc_P1l.units='watt'
            nc_P2l=Origin.createVariable('P_2_lm','f',ncShape3Dlow)
            nc_P2l.description='Potency without noise from beam 2 in low mode'
            nc_P2l.units='watt'
            nc_P3l=Origin.createVariable('P_3_lm','f',ncShape3Dlow)
            nc_P3l.description='Potency without noise from beam 3 in low mode'
            nc_P3l.units='watt'
        if NBeams==5 and c1_op==1:
            nc_P0h=Origin.createVariable('P_0_hm','f',ncShape3D)
            nc_P0h.description='Potency without noise from beam 0, vertical, in high mode'
            nc_P0h.units='watt'
            nc_P1h=Origin.createVariable('P_1_hm','f',ncShape3D)
            nc_P1h.description='Potency without noise from beam 1 in high mode'
            nc_P1h.units='watt'
            nc_P2h=Origin.createVariable('P_2_hm','f',ncShape3D)
            nc_P2h.description='Potency without noise from beam 2 in high mode'
            nc_P2h.units='watt'
            nc_P3h=Origin.createVariable('P_3_hm','f',ncShape3D)
            nc_P3h.description='Potency without noise from beam 3 in high mode'
            nc_P3h.units='watt'
            nc_P4h=Origin.createVariable('P_4_hm','f',ncShape3D)
            nc_P4h.description='Potency without noise from beam 4 in high mode'
            nc_P4h.units='watt'

            nc_P0l=Origin.createVariable('P_0_lm','f',ncShape3Dlow)
            nc_P0l.description='Potency without noise from beam 0, vertical, in low mode'
            nc_P0l.units='watt'
            nc_P1l=Origin.createVariable('P_1_lm','f',ncShape3Dlow)
            nc_P1l.description='Potency without noise from beam 1 in low mode'
            nc_P1l.units='watt'
            nc_P2l=Origin.createVariable('P_2_lm','f',ncShape3Dlow)
            nc_P2l.description='Potency without noise from beam 2 in low mode'
            nc_P2l.units='watt'
            nc_P3l=Origin.createVariable('P_3_lm','f',ncShape3Dlow)
            nc_P3l.description='Potency without noise from beam 3 in low mode'
            nc_P3l.units='watt'
            nc_P4l=Origin.createVariable('P_4_lm','f',ncShape3Dlow)
            nc_P4l.description='Potency without noise from beam 4 in low mode'
            nc_P4l.units='watt'


        

        
        
        
##        in average mode low
        
##        nc_Mask_l=Low.createVariable('Mask_low','f',ncShape2DMeanL)
##        nc_Mask_l.description='mascara en low'
##        nc_Mask_l.units='none'
        ##SPEED COMBINED FRO AIR AND HIDRO
        nc_V1_l=Low.createVariable('W_low','f',ncShape2DL)
        nc_V1_l.description='Vertical air and hidrometeor speed, where negative is upward movement, in low mode'
        nc_V1_l.units='m/s'
        ##SPEED COMBINED FROM AIR
##        nc_V1_air_l=Low.createVariable('W_air_low','f',ncShape2DL)
##        nc_V1_air_l.description='Vertical air speed, where negative is upward movement, in low mode'
##        nc_V1_air_l.units='m/s'
##        ##SPEED COMBINED FROM HIDRO
##        nc_V1_hidro_l=Low.createVariable('W_hidro_low','f',ncShape2DL)
##        nc_V1_hidro_l.description='Vertical hidrometeor speed, where negative is upward movement, in low mode'
##        nc_V1_hidro_l.units='m/s'

        nc_Sig1_l=Low.createVariable('Spectral width_low','f',ncShape2DL)
        nc_Sig1_l.description='Spectral width from the Vertical air speed'
        nc_Sig1_l.units='m/s'

        nc_Sk1_l=Low.createVariable('Skewness_low','f',ncShape2DL)
        nc_Sk1_l.description='Skewness from the Vertical air speed'
        nc_Sk1_l.units='none'

        nc_Kur1_l=Low.createVariable('Kurtosis_low','f',ncShape2DL)
        nc_Kur1_l.description='Kurtosis from the Vertical air speed'
        nc_Kur1_l.units='none'

        nc_Z_l=Low.createVariable('Z_low','f',ncShape2DL)
        nc_Z_l.description='Reflectivity'
        nc_Z_l.units='dBZ'

        nc_U_l=Low.createVariable('U_low','f',ncShape2DL)
        nc_U_l.description='Eastward wind or Zonal wind, data from beam 4(to east) and 5 (to west)'
        nc_U_l.units='m/s'

        nc_V_l=Low.createVariable('V_low','f',ncShape2DL)
        nc_V_l.description='Northward wind or Meridional wind, data from beam 2(to north) and 3 (to south)'
        nc_V_l.units='m/s'

        nc_Snr1_l=Low.createVariable('Snr_low','f',ncShape2DL)
        nc_Snr1_l.description='snr from beam 1 in low mode'
        nc_Snr1_l.units='dB'
        
        nc_Rr_l=Low.createVariable('RR_lm','f4',ncShape2DL)
        nc_Rr_l.description='Rain rate in low mode for each height'
        nc_Rr_l.units='mm h-1'

        nc_EType_l=Low.createVariable('Type_lm','f4',ncShape2DL)
        nc_EType_l.description='EStimation from hydrometeor -10 snow, 0 mixed, 10 rain, and 20 unknown in lm'
        nc_EType_l.units='none'

        nc_EType2_l=Low.createVariable('Type_2_lm','f4',ncShape2DL)
        #nc_EType2_l.description='Estimation from hydrometeor -10 for snow, 0 stratiform rain, and 10 convective rain, 20 unkown'
        nc_EType2_l.description='Estimation from hydrometeor from Ralph -10 for snow, 0 mixed, and 10 convective/stratiform rain, 20 unkown'
        nc_EType2_l.units='none'

        nc_Hkef_l=Low.createVariable('HKEF_lm','f4',ncShape2DL)
        nc_Hkef_l.description='horizontal kinetic energy flux in low mode for each height'
        nc_Hkef_l.units='g s-3'

        nc_Vkef_l=Low.createVariable('VKEF_lm','f4',ncShape2DL)
        nc_Vkef_l.description='Vertical kinetic energy flux in low mode for each height'
        nc_Vkef_l.units='g s-3'

        nc_No_l=Low.createVariable('log10(N0)_lm','f4',ncShape2DL)
        nc_No_l.description='log10(N0) in low mode for each height'
        nc_No_l.units='log10 (m-4-'+u'\u03BC'+' )'

        nc_A_l=Low.createVariable('slope_parameter_lm','f',ncShape2DL)
        nc_A_l.description='Slope parameter in low mode for each height'
        nc_A_l.units='m-1'

        nc_Mu_l=Low.createVariable('shape_parameter_lm','f',ncShape2DL)
        nc_Mu_l.description='shape parameter in low mode for each height'
        nc_Mu_l.units='None'

        nc_Cw_l=Low.createVariable('Liquid_Water_Content_lm','f4',ncShape2DL)
        nc_Cw_l.description='Amount of Liquid Water in low mode for each height'
        nc_Cw_l.units='g m-3'

        nc_C2n_l=Low.createVariable('C2n_lm','f4',ncShape2DL)
        nc_C2n_l.description='Refractive index for low mode in log'
        nc_C2n_l.units='log(m(-2/3))'


        

########O TREC PER ANAR MES RAPID EN LE SPROBES
##The number of beams is the same in low and high mode
    
    if NBeams==1  and c1_op==1:
        nc_P0h[count,:,:]=np.array(np.ma.masked_invalid(PH_c[0]),dtype='f4')
        nc_P0l[count,:,:]=np.array(np.ma.masked_invalid(PL_c[0]),dtype='f4')
    if NBeams==2 and c1_op==1:
        nc_P0h[count,:,:]=np.array(np.ma.masked_invalid(PH_c[0]),dtype='f4')
        nc_P1h[count,:,:]=np.array(np.ma.masked_invalid(PH_c[1]),dtype='f4')

        nc_P0l[count,:,:]=np.array(np.ma.masked_invalid(PL_c[0]),dtype='f4')
        nc_P1l[count,:,:]=np.array(np.ma.masked_invalid(PL_c[1]),dtype='f4')
    if NBeams==3 and c1_op==1:
        nc_P0h[count,:,:]=np.array(np.ma.masked_invalid(PH_c[0]),dtype='f4')
        nc_P1h[count,:,:]=np.array(np.ma.masked_invalid(PH_c[1]),dtype='f4')
        nc_P2h[count,:,:]=np.array(np.ma.masked_invalid(PH_c[2]),dtype='f4')

        nc_P0l[count,:,:]=np.array(np.ma.masked_invalid(PL_c[0]),dtype='f4')
        nc_P1l[count,:,:]=np.array(np.ma.masked_invalid(PL_c[1]),dtype='f4')
        nc_P2l[count,:,:]=np.array(np.ma.masked_invalid(PL_c[2]),dtype='f4')
    if NBeams==4 and c1_op==1:
        nc_P0h[count,:,:]=np.array(np.ma.masked_invalid(PH_c[0]),dtype='f4')
        nc_P1h[count,:,:]=np.array(np.ma.masked_invalid(PH_c[1]),dtype='f4')
        nc_P2h[count,:,:]=np.array(np.ma.masked_invalid(PH_c[2]),dtype='f4')
        nc_P3h[count,:,:]=np.array(np.ma.masked_invalid(PH_c[3]),dtype='f4')

        nc_P0l[count,:,:]=np.array(np.ma.masked_invalid(PL_c[0]),dtype='f4')
        nc_P1l[count,:,:]=np.array(np.ma.masked_invalid(PL_c[1]),dtype='f4')
        nc_P2l[count,:,:]=np.array(np.ma.masked_invalid(PL_c[2]),dtype='f4')
        nc_P3l[count,:,:]=np.array(np.ma.masked_invalid(PL_c[3]),dtype='f4')

    if NBeams==5 and c1_op==1:
        nc_P0h[count,:,:]=np.array(np.ma.masked_invalid(PH_c[0]),dtype='f4')
        nc_P1h[count,:,:]=np.array(np.ma.masked_invalid(PH_c[1]),dtype='f4')
        nc_P2h[count,:,:]=np.array(np.ma.masked_invalid(PH_c[2]),dtype='f4')
        nc_P3h[count,:,:]=np.array(np.ma.masked_invalid(PH_c[3]),dtype='f4')
        nc_P4h[count,:,:]=np.array(np.ma.masked_invalid(PH_c[4]),dtype='f4')

        nc_P0l[count,:,:]=np.array(np.ma.masked_invalid(PL_c[0]),dtype='f4')
        nc_P1l[count,:,:]=np.array(np.ma.masked_invalid(PL_c[1]),dtype='f4')
        nc_P2l[count,:,:]=np.array(np.ma.masked_invalid(PL_c[2]),dtype='f4')
        nc_P3l[count,:,:]=np.array(np.ma.masked_invalid(PL_c[3]),dtype='f4')
        nc_P4l[count,:,:]=np.array(np.ma.masked_invalid(PL_c[4]),dtype='f4')
    
    count=count+1
    if count==2:
        dat_1=int(date2unix(dt.datetime.now()))
        dif_time=(dat_1-dat_0)/60.#minutes

        print('')
        print('The processing time estimated is around '+str(round(dif_time*(int((NumberFiles))),2))+' minutes\n')
            
    if count>2:
        print ('Work in progress( '+str(int(100*Filesdone/NumberFiles))+' %)',end='\r')
            

#########AVOIDE THE BAD SIGNALS
W_h=AvoidBadSignal(W_h,MASK_h)
##W_air_h=AvoidBadSignal(W_air_h,MASK_h);W_hidro_h=AvoidBadSignal(W_hidro_h,MASK_h)
Sig_h=AvoidBadSignal(Sig_h,MASK_h);Sk_h=AvoidBadSignal(Sk_h,MASK_h);Kur_h=AvoidBadSignal(Kur_h,MASK_h)
Z0_h=AvoidBadSignal(Z0_h,MASK_h)
u_h=AvoidBadSignal(u_h,MASK_h);v_h=AvoidBadSignal(v_h,MASK_h)
SnR_h=AvoidBadSignal(SnR_h,MASK_h)
Mu_h=AvoidBadSignal(Mu_h,MASK_h);a_h=AvoidBadSignal(a_h,MASK_h)
HKEF_h=AvoidBadSignal(HKEF_h,MASK_h);VKEF_h=AvoidBadSignal(VKEF_h,MASK_h)
Rr_h=AvoidBadSignal(Rr_h,MASK_h);n0_h=AvoidBadSignal(n0_h,MASK_h)
CW_h=AvoidBadSignal(CW_h,MASK_h);C2N_h=AvoidBadSignal(C2N_h,MASK_h)


##W_h=SmoothMatrix(W_h);W_air_h=SmoothMatrix(W_air_h)

Type_h=Type(Z0_h,W_h,Sig_h,velX,H+h0_opt)

Type2_h=Type2(Z0_h,W_h,Sig_h,velX,H+h0_opt,Sk_h)#s'ha de ficar el sig del vel hidro!!!!!!!!!!!!!!!!!!!!

#nc_V1_h[:,:]=np.array(np.ma.masked_invalid(np.multiply(V1_h,Mask_h)),dtype='f4')
nc_V1_h[:,:]=np.array(np.ma.masked_invalid(W_h),dtype='f4')
nc_Sig1_h[:,:]=np.array(np.ma.masked_invalid(Sig_h),dtype='f4')
nc_Sk1_h[:,:]=np.array(np.ma.masked_invalid(Sk_h),dtype='f4')
nc_Kur1_h[:,:]=np.array(np.ma.masked_invalid(Kur_h),dtype='f4')
nc_Z_h[:,:]=np.array(np.ma.masked_invalid(Z0_h),dtype='f4')
nc_U_h[:,:]=np.array(np.ma.masked_invalid(u_h),dtype='f4')
nc_V_h[:,:]=np.array(np.ma.masked_invalid(v_h),dtype='f4')
nc_Snr1_h[:,:]=np.array(np.ma.masked_invalid(SnR_h),dtype='f4')
nc_Mu_h[:,:]=np.array(np.ma.masked_invalid(Mu_h),dtype='f4')
nc_A_h[:,:]=np.array(np.ma.masked_invalid(a_h),dtype='f4')
nc_Hkef_h[:,:]=np.array(np.ma.masked_invalid(HKEF_h),dtype='f4')
nc_Vkef_h[:,:]=np.array(np.ma.masked_invalid(VKEF_h),dtype='f4')
nc_Rr_h[:,:]=np.array(np.ma.masked_invalid(Rr_h),dtype='f4')
nc_No_h[:,:]=np.array(np.ma.masked_invalid(n0_h),dtype='f4')
nc_Cw_h[:,:]=np.array(np.ma.masked_invalid(CW_h),dtype='f4')
nc_C2n_h[:,:]=np.array(np.ma.masked_invalid(C2N_h),dtype='f4')

##nc_V1_air_h[:,:]=np.array(np.ma.masked_invalid(W_air_h),dtype='f4')
##nc_V1_hidro_h[:,:]=np.array(np.ma.masked_invalid(W_hidro_h),dtype='f4')



nc_EType_h[:,:]=np.array(np.ma.masked_invalid(Type_h),dtype='f4')
nc_EType2_h[:,:]=np.array(np.ma.masked_invalid(Type2_h),dtype='f4')

#########AVOIDE THE BAD SIGNALS
W_l=AvoidBadSignal(W_l,MASK_l)
##W_air_l=AvoidBadSignal(W_air_l,MASK_l);W_hidro_l=AvoidBadSignal(W_hidro_l,MASK_l)
Sig_l=AvoidBadSignal(Sig_l,MASK_l);Sk_l=AvoidBadSignal(Sk_l,MASK_l);Kur_l=AvoidBadSignal(Kur_l,MASK_l)
Z0_l=AvoidBadSignal(Z0_l,MASK_l)
u_l=AvoidBadSignal(u_l,MASK_l);v_l=AvoidBadSignal(v_l,MASK_l)
SnR_l=AvoidBadSignal(SnR_l,MASK_l)
Mu_l=AvoidBadSignal(Mu_l,MASK_l);a_l=AvoidBadSignal(a_l,MASK_l)
HKEF_l=AvoidBadSignal(HKEF_l,MASK_l);VKEF_l=AvoidBadSignal(VKEF_l,MASK_l)
Rr_l=AvoidBadSignal(Rr_l,MASK_l);n0_l=AvoidBadSignal(n0_l,MASK_l)
CW_l=AvoidBadSignal(CW_l,MASK_l);C2N_l=AvoidBadSignal(C2N_l,MASK_l)

##W_l=SmoothMatrix(W_l);W_air_l=SmoothMatrix(W_air_l)

Type_l=Type(Z0_l,W_l,Sig_l,velX,h+h0_opt)

Type2_l=Type2(Z0_l,W_l,Sig_l,velX,h+h0_opt,Sk_l)
#print('rseultat',Type_l)


nc_V1_l[:,:]=np.array(np.ma.masked_invalid(CleanMatrix(W_l)),dtype='f4')
#nc_V1_l[:,:]=np.array(np.ma.masked_invalid(np.multiply(V1_l,Mask_l)),dtype='f4')
nc_Sig1_l[:,:]=np.array(np.ma.masked_invalid(CleanMatrix(Sig_l)),dtype='f4')
nc_Sk1_l[:,:]=np.array(np.ma.masked_invalid(CleanMatrix(Sk_l)),dtype='f4')
nc_Kur1_l[:,:]=np.array(np.ma.masked_invalid(CleanMatrix(Kur_l)),dtype='f4')
nc_Z_l[:,:]=np.array(np.ma.masked_invalid(CleanMatrix(Z0_l)),dtype='f4')
nc_U_l[:,:]=np.array(np.ma.masked_invalid(CleanMatrix(u_l)),dtype='f4')
nc_V_l[:,:]=np.array(np.ma.masked_invalid(CleanMatrix(v_l)),dtype='f4')
nc_Snr1_l[:,:]=np.array(np.ma.masked_invalid(SnR_l),dtype='f4')
nc_Mu_l[:,:]=np.array(np.ma.masked_invalid(Mu_l),dtype='f4')
nc_A_l[:,:]=np.array(np.ma.masked_invalid(a_l),dtype='f4')
nc_Hkef_l[:,:]=np.array(np.ma.masked_invalid(HKEF_l),dtype='f4')
nc_Vkef_l[:,:]=np.array(np.ma.masked_invalid(VKEF_l),dtype='f4')
nc_Rr_l[:,:]=np.array(np.ma.masked_invalid(Rr_l),dtype='f4')
nc_No_l[:,:]=np.array(np.ma.masked_invalid(n0_l),dtype='f4')
nc_Cw_l[:,:]=np.array(np.ma.masked_invalid(CW_l),dtype='f4')
nc_C2n_l[:,:]=np.array(np.ma.masked_invalid(C2N_l),dtype='f4')
##nc_Mask_l[:,:]=np.array(np.ma.masked_invalid(Mask_l),dtype='f4')

##nc_V1_air_l[:,:]=np.array(np.ma.masked_invalid(W_air_l),dtype='f4')
##nc_V1_hidro_l[:,:]=np.array(np.ma.masked_invalid(W_hidro_l),dtype='f4')

nc_EType_l[:,:]=np.array(np.ma.masked_invalid(Type_l),dtype='f4')
nc_EType2_l[:,:]=np.array(np.ma.masked_invalid(Type2_l),dtype='f4')


dataset.close()

print('starts at ',tt00,'end at ',dt.datetime.now(),end='\r')
print('\nfiles number processed',count)
