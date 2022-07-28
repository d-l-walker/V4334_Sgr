##############################################################
# Based on http://almanas.jb.man.ac.uk/amsr/Atomium/SelfCal.py
# see http://atomium.pbworks.com/w/page/131757642/Assignments
# amsr July 2022
#
# This version is for 2021.1.00917.S which is in TDM (i.e. 128 x 15.625 MHz channels per spw)
# with two sets of spectral scans in each of Bands 6 and 7
# The starting point is the target split out from the pipeline data for each set of scans
# The data for each band are self-calibrated first, and then the bands are combined for final rounds.
# I suggest that you make and work in an imaging directory, e.g.
# /raid/scratch/aaz/2021.1.00917.S/Imaging  (or you could further subdivide for bands etc.)
# Further rounds of self-cal can be added if needed
#
# The final round of applycal makes a copy of the data before calibrating weights.
# Failed solutions are not flagged as all data have already had phase ref solutions and
# some spw just have less continuum so lower S/N
#
# 1) You might want solnorm=True in the 'a' gaincal as the final flux density
#    and S/N is very slightly less than 'p' only
#    - however, the spectral slope of the continuum is less wiggly.
#
# 2) The different SG had somewhat different max baselines, hence different resolutions.
#    At present the final cube has restoringbeam = 'common' which forces it to the largest size present,
#    which at least means all the B6 have the same resolution - but you could force it to be a bit
#    smaller if you wanted e.g. for consistently with B7
#
# 3) The automasking is quite cautious, you might want to reduce the thresholds a bit or otherwise make
#    it more assertive .
#    Also, if you are doing it interactively, once the automasking residual is <~1 mJy, I suggest
#    masking the central ~1" for all planes manually, to pick up weak central emission.
#    If you restart tclean, set calcres=False, calcpsf=False.
#
##############################################################
# Check for best CASA version, currently /usr/local/casa-6.3.0-48/bin/casa

# Please enter user-defined parameters (some already done):
pth = '/raid/scratch/aaz/2021.1.00917.S/'
target='V4334_Sgr'

# Band 6

msin=[pth+'Band_6/MSes/band_6_lower_half.s',pth+'Band_6/MSes/band_6_upper_half.ms']
mvis=target+'_B6.ms'

# Enter a list of refants, e.g. the first few in the pipeline priority lists for each dataset (weblog task hif_refant).
antrefs='DA63, DA57, DV08, DV24'

Nspw=40
rms_1_1 = 0.25  # mJy rms noise predicted mid B6 242 GHz for 1 GHz, 1 min. Dec -17 from sensitivity calculator
# Time for sensitivity calculation (should be able to get from data but utility not compatible at present)
ToS= 100./10.    #100 min over 10 tunings
# Solution intervals - refine after step 4
# Ideally there should be approx. integer solints per scan and solint should be ~integer multiple of integration time
# The actual minimum solint may be >> theoretical as some spw have fewer line-free channels.
solintp0='66s'     # increase if less than minsolint
solintp1='21s'     # no less than minsolint
solinta1='66s'     # usually longer than solintp1

contchans='17:0~88;101~112,19:0~123, 21:0~9;34~55;59~72;85~89;97~123, 23:0~17;38~90,\
29:0~90;93~123,31:0~8;33~40;74~123, 33:25~123,35:0~87;95~106;114~120,\
41:0~41;51~63;68~85;116~123, 43:0~48;90~101,45:0~88;108~123,47:0~33;71~90;107~123,\
53:0~123,55:0~120,57:12~24;54~95,59:15~44;65~80;93~104;116~123,\
65:46~58;75~80;86~89;97~123,67:0~54;71~78;122~123,69:0~27;77~93;108~123,71:0~123,\
133:0~34;52~76;95~123,135:15~123,137:0~104,139:37~80;103~111,\
145:0~9;18~26;34~42;64~78;84~89;99~111,147:0~39;56~77;92~96,149:0~33;72~83;118~123,151:14~29;38~45;50~86;115~123,\
157:0~123,159:26~123,161:0~32;41~58;69~101,163:45~101;108~115,\
169:0~110,171:0~8;32~79;91~123,173:0~5;37~123,175:0~7;24~64;73~123,\
181:0~51;69~76;91~117,183:0~123,185:0~58;115~123,187:0~20;52~123'

#########################################

# repeat for B7
#msin=[]
#mvis=target+'_B7.ms'
#antrefs=''
#contchans=''
#rms_1_1 = 0.4  # mJy rms noise mid B7 for 1 GHz, 1 min. Dec -17
#ToS= 150./10.   # 150 min over 10 tunings

##########################################
# Parameters for all

cell='0.05arcsec'
R=0.75
imsz=540
contchwd=0.015626283 # GHz (in fact all channels, here)
Ncorr = 2.           # 2 polarizations

##################################
# Cube parameters

spwsel=''    # to do all - if so, see parameters in step 13 to avoid excessively noisy or blank channels
                              # which bias automasking and beamsize **change for B7?
# spssel='21'  # comment out to do all or select test channel
     # One first to test, or '' to do all
interact=False   # to check progress/tweak masks manually
masktype='auto-multithresh' # 'automultithresh' for automasking or 'user' for mask named in mask parameter or default
thismask=''     # default name will be used but you can specify an existing mask
nthresh=3.5     # rms multiplier for initial automasking
mbf=0.3         # min beam fraction for automasking

# open file to send continuum stats to.
statfile=open(mvis+'_imstats1.txt', 'a')

#***
# After each continuum imaging step (in the next or the same step) the script
# writes the peak, rms, S/N to screen and to a text file target*_imstats.txt
# This is cumulative including repeated steps so you can check for improvement.
#
# set thesteps here (e.g. to only do band 6, do [0], 1
# or in terminal mysteps=[0] etc. (will overrule thesteps)
thesteps = [1]
step_title = {0:'Concatenate and list',
              1:'Print time on source, predicted sensitivity etc.',
              3:'First continuum image',
              4:'Get image statistics and predict gaincal solint',
              5:'First phase-only self-calibration, plotms solutions',
              6:'Apply calibration and re-image',
              7:'Next phase-only self-calibration, plotms solutions',
              8:'Apply calibration and re-image',
              9:'First amp self-calibration, plotms solutions',
              10:'Apply calibration and re-image',
              12:'Subtract continuum and check line visibilities',
              13:'Make image cube',
              15: 'Export images to FITS format'}


import analysisUtils as aU
from os import write

# Estimate continuum sensitivity
def totcont(contchans):
#    for s in sps:
#        if s in contchans:
#            contchans=contchans.replace(s, '(-')
    contchans=re.sub(r'.?.?:','(-', contchans)                    # replace spw numbers and :
    contchans=contchans.replace('1(','(')                         # remove 1 from spw numbers >100
    contchans=contchans.replace('~','+').replace(';',')+(-').replace(',',')+')
    contchans=contchans+(')')
    return contchans

if len(contchans) > 1:
    ncontchan=eval(totcont(contchans))
    contGHz=contchwd*ncontchan               # total line-free continuum

# sensitivity as ideal * sqrt(ncontchan/maxcontchan)
if len(contchans) > 1:
  if os.path.exists(mvis):
      try: Nants
      except NameError:
          Nants=len(aU.getUnflaggedAntennas(mvis)[0])
#      try: ToS    # utility does not work at present, ToS hard-wired
#      except NameError:
#          ToS=aU.timeOnSource(msin[0],field=target,scienceSpwsOnly=True)[0]['minutes_on_source_with_latency'] # average time on target
      predicted_rms_cont=rms_1_1*(ToS*contGHz)**-0.5
      predicted_rms_line = predicted_rms_cont*(float(ncontchan))**0.5 # potential sensitvity per chan

# Set thresholds for imaging
if 'predicted_rms_cont' in locals():
    thresh= '%2.3f mJy' %(predicted_rms_cont)

if 'predicted_rms_line' in locals():
    linethresh= '%2.3f mJy' %(predicted_rms_line)

###############
try:
  print( 'List of steps to be executed ...', mysteps)
  thesteps = mysteps
except:
  print( 'global variable mysteps not set.')
if (thesteps==[]):
  thesteps = range(0,len(step_title))
  print( 'Executing all steps: ', thesteps)

##########

# concatenate
mystep = 0
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print( 'Step ', mystep, step_title[mystep])

  os.system('rm -rf '+mvis)
  concat(vis=msin,
         concatvis=mvis,
         freqtol='2MHz',
         dirtol='1arcsec',
         copypointing=False)


  listobs(vis=mvis,
          verbose=True,
          overwrite=True,
          listfile=mvis+'.listobs')

# Print useful information
mystep = 1
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print( 'Step ', mystep, step_title[mystep])

  print('concatenated MS '+mvis)
  print( 'Total line-free continuum %5.2f GHz' %(contGHz))
  print( '\nTime on source (avg. per tuning)  %5.1f min' %(ToS))
  print( str(Nants)+'  antennas present (not necessarily in all spw).')
  print( 'Predicted continuum rms  %5.3f mJy' %(predicted_rms_cont))
  print( 'This will be set as "thresh", the continuum cleaning threshold.  \n**Change predicted_rms_cont as needed: higher in early stages or lower at the end for good conditions.\n\n')
  print( 'Predicted line rms  %5.3f mJy' %(predicted_rms_line))
  print( 'This will be set as "linethresh" for cube cleaning and automasking, tweak as needed.\n')

  statfile.write('concatenated MS '+mvis+'\n')
  statfile.write( 'Total line-free continuum %5.2f GHz\n' %(contGHz))
  statfile.write( '\nTime on source (avg. per tuning)  %5.1f min\n' %(ToS))
  statfile.write( str(Nants)+'  antennas present (not necessarily in all spw).\n')
  statfile.write( 'Predicted continuum rms  %5.3f mJy\n' %(predicted_rms_cont))
  statfile.write( 'Predicted line rms  %5.3f mJy\n' %(predicted_rms_line))

# # Plot to check continuum selection
# mystep = 2
# if(mystep in thesteps):
#   casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
#   print( 'Step ', mystep, step_title[mystep])
#
#   chs=''
#   if len(contchans) > 0:
#       chs=contchans
#
#   plotms(vis=mvis,
#          xaxis='channel', yaxis='amp',
#          spw=chs,
#          avgtime='999999',
#          avgscan=True,
#          avgbaseline=True,
#          coloraxis='observation',#showtsky=True,
#          iteraxis='spw')


# first image continuum
mystep = 3
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print( 'Step ', mystep, step_title[mystep])

  print('Stop cleaning when residual is all noise, maybe higher than threshold\n')

  clearcal(vis=mvis,
           addmodel=True)                             # **

  os.system('rm -rf '+mvis+'_cont.clean*')
  tclean(vis=mvis,
         imagename=mvis+'_cont.clean',
         spw=contchans,
         imsize=imsz,
         cell=cell,
         weighting = 'briggs',
         robust=0.5,
         interactive=False,
         perchanweightdensity=False,
         threshold=thresh,
         niter=150,
         savemodel='modelcolumn',
         parallel=True,
         usemask=masktype,
         sidelobethreshold=2.5,
         noisethreshold=nthresh,
         minbeamfrac=mbf,
         lownoisethreshold=1.2,
         growiterations=75)

  print('Check in logger or by plotting that the model is saved.\n\n')

# Stats and prediction
mystep = 4
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print( 'Step ', mystep, step_title[mystep])

  rms=imstat(imagename=mvis+'_cont.clean.image',
             box='10,10,'+str(imsz-10)+','+str(0.2*imsz))['rms'][0]
  peak=imstat(imagename=mvis+'_cont.clean.image',
              box=str(0.25*imsz)+','+str(0.25*imsz)+','+str(0.75*imsz)+','+str(0.75*imsz))['max'][0]
  print( mvis+'_cont.clean.image')
  print( 'Predicted continuum rms  %5.3f mJy; achieved %5.3f mJy' %(predicted_rms_cont, rms*1000.))
  print( 'Continuum peak  %5.3f mJy; S/N %5i' %(peak*1000., peak/rms))
  print( 'You might want to reduce threshold if the rms is lower.')

  statfile.write('\n'+mvis+'_cont.clean.image\n')
  statfile.write('Predicted continuum rms  %5.3f mJy; achieved %5.3f mJy\n' %(predicted_rms_cont, rms*1000.))
  statfile.write('Continuum peak  %5.3f mJy; S/N %5i\n' %(peak*1000., peak/rms))

  N=float(Nants)
  minsol=60.*(3.*predicted_rms_cont/(peak*1000.))**2 *ToS*(N*(N-1)/(2*(N-3)))*Ncorr*float(Nspw)
  print( 'Minimum solint for S/N 3 per antenna, per spw, per polarization %5.1f sec' %(minsol) )
  statfile.write( 'Minimum solint for S/N 3 per antenna, per spw, per polarization %5.1f sec\n' %(minsol) )

  print( 'Note these values (also and enter solintp0, solintp1, solinta1 values in variables at start\n')
  print( 'If min. solint is very short, start with a longer solint e.g. 60s, to refine model first\n')

  image1=mvis+'_cont.clean.image'
  fit=imfit(imagename=image1)['results']['component0']['shape']
  fitcenRA=fit['direction']['m0']['value']
  if fitcenRA < 0:
      fitcenRA=fitcenRA+2.*pi
  fitcenDec=fit['direction']['m1']['value']
  fitcen='ICRS '+str(fitcenRA)+'rad  '+str(fitcenDec)+'rad'
  fitcenradec=aU.rad2radec(fitcenRA, fitcenDec)
  imcen=aU.rad2radec(imhead(image1,mode='get',hdkey='crval1')['value'], imhead(image1,mode='get',hdkey='crval2')['value'])
#  fitcendegra= aU.rad2deg(fit['direction']['m0']['value'])
#  if fitcendegra <0:
#      fitcendegra=360.+fitcendegra

# raoff=3600000.*degrees(fitcenRA-imhead(image1,mode='get',hdkey='crval1')['value'])*cos((fitcenDec+imhead(image1,mode='get',hdkey='crval2')['value'])/2)
#  decoff=3600000.*degrees(fitcenDec-imhead(image1,mode='get',hdkey='crval2')['value'])

  statfile.write( '\npeakpos = "'+fitcen+'"\n')
  statfile.write( '='+imcen+'\n')
#  statfile.write( 'Peak position - pointing centre = (%4.1f, %4.1f) mas' % (raoff, decoff))

  print( 'peakpos = "'+fitcen+'"' )
  print( '='+imcen+'\n')
#  print( 'Peak position - pointing centre = (%4.1f, %4.1f) mas' % (raoff, decoff))


# first phase self-cal
mystep = 5
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])

#  ft(vis=mvis,
#    model=target+'_'+config+'_cont.clean.model',
#    usescratch=True)

  os.system('rm -rf '+mvis+'.p0')
  gaincal(vis=mvis,
          field=target,
          caltable=mvis+'.p0',
          spw=contchans,
          solint=solintp0,
          refant=antrefs,
          refantmode='flex',
          calmode='p')

  # plotms(vis=mvis+'.p0',
  #        xaxis='time',yaxis='phase',
  #        xsharedaxis=True,xselfscale=True,
  #        ysharedaxis=True,yselfscale=True,
  #        gridrows=5,gridcols=2,
  #        iteraxis='antenna')

# applycal and re-image
mystep = 6
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])

# flag data seen to have bad solutions
  flagdata(vis=mvis,
           mode='manual',
           antenna='DA50',
           scan='48',
           spw='17,19,21,23')

# Do not flag nor weight, to avoid a imperfect model biasing results
  applycal(vis=mvis,
           gaintable=mvis+'.p0',
           applymode='calonly',
           calwt=False,
           flagbackup=False)

  os.system('rm -rf '+mvis+'_contp0.clean*')
  tclean(vis=mvis,
         imagename=mvis+'_contp0.clean',
         spw=contchans,
         imsize=imsz,
         cell=cell,
         weighting = 'briggs',
         robust=0.5,
         interactive=False,
         perchanweightdensity=False,
         threshold=thresh,
         niter=200,
         savemodel='modelcolumn',
         parallel=True,
         usemask=masktype,
         sidelobethreshold=2.5,
         noisethreshold=nthresh,
         minbeamfrac=mbf,
         lownoisethreshold=1.2,
         growiterations=75)


  rms=imstat(imagename=mvis+'_contp0.clean.image',
             box='10,10,'+str(imsz-10)+','+str(0.2*imsz))['rms'][0]
  peak=imstat(imagename=mvis+'_contp0.clean.image',
              box=str(0.25*imsz)+','+str(0.25*imsz)+','+str(0.75*imsz)+','+str(0.75*imsz))['max'][0]

  print(mvis+'_contp0.clean.image\n')
  print( 'Predicted continuum rms  %5.3f mJy; achieved %5.3f mJy' %(predicted_rms_cont, rms*1000.))
  print( 'Continuum peak  %5.3f mJy; S/N %5i' %(peak*1000., peak/rms))
  print( 'You might want to reduce threshold if the rms is lower.')

  statfile.write('\n'+mvis+'_contp0.clean.image\n')
  statfile.write('Predicted continuum rms  %5.3f mJy; achieved %5.3f mJy\n' %(predicted_rms_cont, rms*1000.))
  statfile.write('Continuum peak  %5.3f mJy; S/N %5i\n' %(peak*1000., peak/rms))

######################################
# second phase self-cal
mystep = 7
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])

  os.system('rm -rf '+mvis+'.p1')
  gaincal(vis=mvis,
          field=target,
          caltable=mvis+'.p1',
          gaintable=mvis+'.p0',
          gaintype='T',    # average polarisations to allow shorter solint
          spw=contchans,
          solint=solintp1,
          refant=antrefs,
          refantmode='flex',
          calmode='p')

  # plotms(vis=mvis+'.p1',
  #        xaxis='time',yaxis='phase',
  #        xsharedaxis=True,xselfscale=True,
  #        ysharedaxis=True,yselfscale=True,
  #        gridrows=5,gridcols=2,
  #        iteraxis='antenna')

# applycal and re-image
mystep = 8
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])

# Do not flag nor weight, to avoid a imperfect model biasing results
  applycal(vis=mvis,
           gaintable=[mvis+'.p0',mvis+'.p1'],
           applymode='calonly',
           calwt=False,
           flagbackup=False)

  os.system('rm -rf '+mvis+'_contp1.clean*')
  tclean(vis=mvis,
         imagename=mvis+'_contp1.clean',
         spw=contchans,
         imsize=imsz,
         cell=cell,
         deconvolver='mtmfs', nterms=2,
         weighting = 'briggs',
         robust=0.5,
         interactive=False,
         perchanweightdensity=False,
         threshold=thresh,
         niter=300,
         savemodel='modelcolumn',
         parallel=True,
         usemask=masktype,
         sidelobethreshold=2.5,
         noisethreshold=nthresh,
         minbeamfrac=mbf,
         lownoisethreshold=1.2,
         growiterations=75)

  rms=imstat(imagename=mvis+'_contp1.clean.image.tt0',
             box='10,10,'+str(imsz-10)+','+str(0.2*imsz))['rms'][0]
  peak=imstat(imagename=mvis+'_contp1.clean.image.tt0',
              box=str(0.25*imsz)+','+str(0.25*imsz)+','+str(0.75*imsz)+','+str(0.75*imsz))['max'][0]

  print(mvis+'_contp1.clean.image.tt0\n')
  print( 'Predicted continuum rms  %5.3f mJy; achieved %5.3f mJy' %(predicted_rms_cont, rms*1000.))
  print( 'Continuum peak  %5.3f mJy; S/N %5i' %(peak*1000., peak/rms))

  statfile.write('\n'+mvis+'_contp1.clean.image.tt0\n')
  statfile.write('Predicted continuum rms  %5.3f mJy; achieved %5.3f mJy\n' %(predicted_rms_cont, rms*1000.))
  statfile.write('Continuum peak  %5.3f mJy; S/N %5i\n' %(peak*1000., peak/rms))

######################
# amp self-cal
mystep = 9
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])

  os.system('rm -rf '+mvis+'.a1')
  gaincal(vis=mvis,
          field=target,
          caltable=mvis+'.a1',
          gaintable=[mvis+'.p0',mvis+'.p1'],
          spw=contchans,
          solint=solinta1,
          refant=antrefs,
          refantmode='flex',
          calmode='a')

  # plotms(vis=mvis+'.a1',
  #        xaxis='time',yaxis='amp',
  #        xsharedaxis=True,xselfscale=True,
  #        ysharedaxis=True,yselfscale=True,
  #        gridrows=5,gridcols=2,
  #        iteraxis='antenna')

# applycal and re-image
mystep = 10
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])

  if not os.path.exists(mvis+'_prefinalapplycalwt'):
      os.system('cp -r '+mvis+' '+mvis+'_prefinalapplycalwt')
  applycal(vis=mvis,
           gaintable=[mvis+'.p0',mvis+'.p1',mvis+'.a1'],
           applymode='calonly')  # see notes at start
#           calwt=False,          # see notes at start
#           flagbackup=False)

  os.system('rm -rf '+mvis+'_contpa1.clean*')
  tclean(vis=mvis,
         imagename=mvis+'_contpa1.clean',
         spw=contchans,
         imsize=imsz,
         cell=cell,
         deconvolver='mtmfs', nterms=2,
         weighting = 'briggs',
         robust=0.5,
         interactive=False,
         threshold=thresh,
         niter=300,
         parallel=True,
         usemask=masktype,
         sidelobethreshold=2.5,
         noisethreshold=nthresh,
         minbeamfrac=mbf,
         lownoisethreshold=1.2,
         growiterations=75)

  rms=imstat(imagename=mvis+'_contpa1.clean.image.tt0',
             box='10,10,'+str(imsz-10)+','+str(0.2*imsz))['rms'][0]
  peak=imstat(imagename=mvis+'_contpa1.clean.image.tt0',
              box=str(0.25*imsz)+','+str(0.25*imsz)+','+str(0.75*imsz)+','+str(0.75*imsz))['max'][0]

  print(mvis+'_contpa1.clean.image.tt0\n')
  print( 'Predicted continuum rms  %5.3f mJy; achieved %5.3f mJy' %(predicted_rms_cont, rms*1000.))
  print( 'Continuum peak  %5.3f mJy; S/N %5i' %(peak*1000., peak/rms))

  statfile.write('\n'+mvis+'_contpa1.clean.image.tt0\n')
  statfile.write('Predicted continuum rms  %5.3f mJy; achieved %5.3f mJy\n' %(predicted_rms_cont, rms*1000.))
  statfile.write('Continuum peak  %5.3f mJy; S/N %5i\n' %(peak*1000., peak/rms))

# # Check continuum selection,  no spw entirely flagged, corrected continuum slope more regular
# mystep = 11
# if(mystep in thesteps):
#   casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
#   print( 'Step ', mystep, step_title[mystep])
#
#   chs=''
#   if len(contchans) > 0:
#       chs=contchans
#
#   plotms(vis=mvis,
#          xaxis='freq', yaxis='amp',
#          ydatacolumn='corrected',
#          spw=chs,
#          avgtime='999999',
#          avgscan=True,
#          avgbaseline=True,
#          coloraxis='spw')

# uvcontsub
mystep = 12
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])

  os.system('rm -rf '+mvis+'.contsub')
  os.system('rm -rf '+mvis+'.contsub.flagversions')
  mstransform(vis= mvis,
            douvcontsub=True,
            reindex=False,    #*** CHECK WORKS should keep original spw numbering
            outputvis=mvis+'.contsub',
            fitspw=contchans,
            fitorder=1)

  # plotms(vis=mvis+'.contsub',
  #        xaxis='freq',
  #        yaxis='amp',
  #        ydatacolumn='corrected',
  #        avgtime='999999',
  #        avgscan=True,
  #        avgbaseline=True,
  #        coloraxis='spw')

# image cubes
# might need to experiment with automasking and threshold...
mystep = 13
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])

  st=0
  nch=-1
  if spwsel=='' or spwsel=='17':
      st=5
  if spwsel=='':
      nch=4040

  os.system('rm -rf '+mvis+'_spw'+spwsel+'.clean*')
  tclean(vis=mvis+'.contsub',
         imagename=mvis+'_spw'+spwsel+'.clean',
         spw=spwsel,      # NB if blank=all, set start & nchan as per preamble to avoid bad/blank end chans
         start=st, nchan=nch,
         specmode='cube',
         outframe='lsrk',
         perchanweightdensity=False,
         imsize=imsz,
         cell=cell,
         mask=thismask,
         weighting = 'briggs',
         robust=0.5,
         pbcor=True,
         interactive=False,
         restoringbeam='common',    # otherwise differences between tunings.  Or can set explicitly
         cycleniter=200,
         usemask=masktype,
         sidelobethreshold=2.5,     # increase for much bright emission
         noisethreshold=nthresh,    #
         minbeamfrac=mbf,           #
         lownoisethreshold=1.2,     # increase for much bright emission
         growiterations=75,         # for speed
         threshold=linethresh,
         niter=50000,
         parallel=True)               # increase if needed



statfile.close()
