#### Manual extraction pipeline. Extracts arc, identifies lines, extracts spectrum for multiple file and wavelength calibrates. 
#### Flux calibrates to sensitivity function, then corrects to spectrostandard. 
#### Automatically removes previous extraction files when extract_spectra = 'y'

from pyraf import iraf
import os
import glob
import numpy as np
from astropy.io import fits

############################
obj = 'AT2019'+''
z = 0.0
extract_arc = 'y'
extract_spectra = 'y'
flux_calibrate  = 'y'
apply_flux_correction = 'y'
plotgal = 'y'
plotvel = 'n'
velocity = 15500
E = 1e-5
############################

####### Define function
def dop(lamr,vel):
    vel=vel*-1
    a = (1.+vel/3.e5)**0.5
    b=(1.-vel/3.e5)**0.5
    return lamr*a/b

def dered(wav,flux,Rv,Ebv):
   lam=wav*0.0001
   Av=Ebv*Rv
   x=1/lam
   y=x-1.82
   a=1+(0.17699*y)-(0.50477*y**2)-(0.02427*y**3)+(0.72085*y**4)+(0.01979*y**5)-(0.77530*y**6)+(0.32999*y**7)
   b=(1.41338*y)+(2.28305*y**2)+(1.07233*y**3)-(5.38434*y**4)-(0.62251*y**5)+(5.30260*y**6)-(2.09002*y**7)
   AlAv=a+b/Rv
   Al=AlAv*Av
   #print Al
   F=10**(Al/2.5)
   delF= flux*F
   return delF
   
def isarc(fits_file):
    hdul = fits.open(fits_file)
    hdr = hdul[0].header
    if 'OBSTYPE' in hdr:
        if hdr['OBSTYPE']=='ARC':
            return True
    else:
        return False 

### misc definitions
spectrafilelist='@speclist'
sensfile='newsen'
hthresh2=90000


### make a list of the files
## create the speclist but importing all fits files in folder but rejecting the sensfunc file, arc, and previous extraction files
file_list = glob.glob('./*.fits');file_list.sort()

### find arcfile
for f in file_list:
    if isarc(f) == True:
        arcfile = f
        
### find object file        
final=[]
for j in range(len(file_list)):
    if sensfile not in file_list[j]:
        if arcfile not in file_list[j]:
            if '.ms.fits' not in file_list[j]:
                if '.w.fits' not in file_list[j]: 
                    final.append(file_list[j])
            
        
### Save the list of spectra        
np.savetxt('./speclist',final,fmt="%s")

### apall on arc
iraf.twodspec()
iraf.apextract(dispaxis=1)

if extract_arc !='n':
    print('\n'+'Extracting arc file')
    iraf.apall(input=arcfile, output='arc.ms',apertur=1,interac='yes',find='yes',nfind=1,recente='yes',resize='yes',background='none')

    ## identify lines
    print('\n'+'Identify arc lines')
    iraf.identify(images='arc.ms',coordlist='Xe.txt', function='spline3',order=3,maxfeatures=10)

if extract_spectra !='n':
  
    ### apall on spectrum
    for n in final:
        
        ### Remove previous extractions
        if os.path.exists(n+'.ms.fits')==True:
            print('Removing',n+'.ms.fits')
            os.remove(n+'.ms.fits')
        if os.path.exists(n+'.w.fits')==True:
            print('Removing',n+'.w.fits')
            os.remove(n+'.w.fits')
                    
        ### Commence extraction
        print('\n'+'Spectrum file = '+n+'[0]')
        print('\n'+'Extracting Spectra...')
        iraf.apall(input=n+'[0][0]', output=n+'.ms',apertur=1,interac='yes',find='yes',nfind=1,recente='yes',resize='yes',background='fit')
    
        print('\n'+'Editing the header to point to the wavelength calibration...')
        iraf.hedit(images=n+'.ms', addonly='yes',fields='REFSPEC1', value='arc.ms')
    
        print('\n'+'Applying the dispersion correction...')
        iraf.dispcor(input=n+'.ms',output=n+'.w')
    
        print('\n'+'Wavelength calibrated spectrum output as '+ n+'.w.fits')
    

print('\nBeginning flux calibration...')
    
## create the speclist but importing all fits files in folder but rejecting the sensfunc file
file_list = glob.glob('./*.fits');file_list.sort()
final=[]
for j in range(len(file_list)):
    if 'fits.w.fits' in file_list[j]: 
        final.append(file_list[j])
        
np.savetxt('./speclist',final,fmt="%s")

## extract the date. file_list[3] is used because arc.fits, arc.ms.fits, and the sensfile take up the first three positions
date=''
it=0
for j in range(len(file_list[3])):
    if file_list[3][j]=='_':
        it+=1
    if it == 2 and file_list[3][j]!='_':
        date = date+file_list[3][j]
print('Date = '+ date)

if flux_calibrate != 'n':    
    ## begin IRAFing
    
    print('\n'+'Setting airmass...')
    iraf.setairmass('@speclist', observ='lapalma', ra='cat-ra', dec='cat-dec', equi='cat-equi', st='lst', ut='utstart') 
    
    print('\n'+'Combining spectra...')        
    iraf.scombine(input='@speclist',output='allspec',group='all',combine='median', gain = 2.45, reject='crreject',lthresh=1e-30, hthresh=hthresh2)
    
    print('\n'+'Applying flux calibration to sensitivity function...')   
    iraf.calibrate('@speclist', 'allspeccal', obs='lapalma', sens=sensfile,extinct='no', ignoreaps='yes') ## If extinct='yes', try extinction ='onedstds$ctioextinct.dat'
    
    print('\n'+'Extracting the 1D spectra from the calibrated multispec file...') 
    iraf.scopy('allspeccal[0]', 'allspeccal2', format='onedspec', w1=4000,w2=8000, rebin='no')
    
    
    print('\n'+'1D spectrum output...')
    iraf.wspectext('allspeccal2.0001', obj+'_'+date+'.txt',header='no')
    
    print('\n'+'Tidying up...')
    a = glob.glob('./*')
    for f in a:
        if 'allspec' in f:
            os.remove(f)
    
if apply_flux_correction != 'n':
        
    s = np.loadtxt('./ManualExtractionFluxCorrections.txt',unpack=True, usecols=(0,1))
    
    file_to_use = obj+'_'+date+'.txt'
    print('\n'+'Using '+ file_to_use+' as reference...')
    a=np.loadtxt('./'+file_to_use,unpack=True)
    
    for q in range(len(a[0])):
        find = np.argmin(abs(s[0]-a[0][q]))
        a[1][q] = a[1][q]*s[1][find]

         
    print('\n'+file_to_use + ' corrected to ManualExtractionFluxCorrections.txt')
    master = zip(a[0],a[1])

    ### Save the corrected spectrum
    np.savetxt('./'+file_to_use[:-4]+'.w.txt',master,fmt="%s")
        
    ### Load the corrected spectrum
    n = np.loadtxt('./'+file_to_use[:-4]+'.w.txt',dtype='float',unpack=True)
    
    ### Plot the corrected and dereddened spectrum
    import matplotlib.pyplot as plt
    m = np.max(a[1])
    a[1] = dered(a[0],a[1],3.1,E)
    plt.plot(a[0]/(1.+z),a[1]/max(a[1]),color='k',linewidth=1,label=obj+'\n$z='+str(z)+'$')
    
    ### if necessary plot the galaxy lines
    if plotgal != 'n':
        l = [6548,6583, 4959,5007,5890,6717,6731]
        for line in l:
            ys = np.linspace(0,1.)
            xs = [line for op in ys]
            plt.plot(xs,ys,color='grey',linestyle='dotted',zorder=0,linewidth=0.7)
        l = [6563,4861,4341,4100]
        for line in l:
            ys = np.linspace(0,1.)
            xs = [line for op in ys]
            plt.plot(xs,ys,color='red',linestyle='dashed',zorder=0,linewidth=0.7)
    
    if plotvel !='n':
        l = [6355]
        for line in l:
            ys = np.linspace(0,1.)
            xs = [dop(line,velocity) for op in ys]
            plt.plot(xs,ys,color='green',linestyle='dashed',zorder=0,linewidth=0.7)
                    
    
    plt.legend()
    plt.ylim([0,1.1])
    plt.xlabel('Rest-frame wavelength [$\AA$]')
    plt.ylabel('Scaled flux')
    plt.savefig('./'+obj+'.pdf',bbox_inches='tight')
    plt.close()
    
    
    ### Compare the two spectra
    a = np.loadtxt(file_to_use,unpack=True)
    plt.plot(a[0]/(1.+z),a[1]/max(a[1]),color='k',alpha=1,linewidth=1,label='Original')
    a = np.loadtxt(file_to_use[:-4]+'.w.txt',unpack=True)
    plt.plot(a[0]/(1.+z),a[1]/max(a[1]),color='red',alpha=0.7,linewidth=1,label='Corrected')
    plt.legend()
    plt.xlabel('Rest-frame wavelength [$\AA$]')
    plt.ylabel('Scaled flux')
    plt.savefig('./'+obj+'_compare.pdf',bbox_inches='tight')
    plt.close()