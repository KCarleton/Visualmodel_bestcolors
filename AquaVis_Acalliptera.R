######################################
#
#     AquaVis_Massoko R set of functions for visual modeling
#
#####################################
# Karen Carleton May 2024          
#
# This gets rid of wavelength in VP file
#
# This library contains a set of functions that are useful for modeling vision
# in an aquatic environment: AquaVis.  It calculates quantum catch
# for colors illuminated by light taking into account a visual system with 
# Short, Medium, and Long wavelength sensitive cones and a transparent lens.  
# It then calculates the quantum catches for the three receptors and uses 
# the Receptor Noise Limited model to determine the
# Just Noticeable Differences (JNDs) between those colors to determine what is 
# discriminable. 
# Refs: Vorobyev and Osorio 1998 for RNL  
# Champ et al 2016 and Olsson et al 2015 for the absolute calculations
# The functions contained in this version include: 
#   This includes 
#     Opsin_gov to calculate visual pigment templates 
#     Absorptance_calc to determine absorptance
#     VPAbsolute_calc to determine absolute absorptance (incl photoreceptor diameter and length)
#     Mono_dataset to calculate a set of monochromatic colors for comparison in SML space
#     Peak_data or Step_data to calculate a single peak or single step color
#     Qcatch_calc to calculate quantum catch of S, M and L cones for a set of colors
#     XY_plotS (or L) to plot SML data in an SML space (S or L cone plotted on top vertex)
#     JND_calc to calculate JNDs for a set of colors
#
#================================================================================
#   FUNCTION: Opsin_gov 
#      Calculate visual pigments absorption using Govardovskii et al 2000 templates
#           lbeg, lend = wavelength range
#           lmax = peak wavelength
#           A1_chrom = % of A1 chromophore
#--------------------------------------------------------------------------------
Opsin_gov<-function(lbeg, lend, lmax, A1_chrom)  {
  # fit A1 peak
  a1<-0.8795+0.0459*exp(-(lmax-300)^2/11940)
  lmb1<-189+0.315*lmax
  p1<--40.5+0.195*lmax  #Govardovskii calls this b but this confuses it with the alpha peak equation
  i<-lbeg:lend
  x<-lmax/i
  Salpha_A1<-1/((exp(69.7*(a1-x))+exp(28*(0.922-x))+exp(-14.9*(1.104-x))+0.674))
  Sbeta_A1<-0.26*exp(-((i-lmb1)/p1)^2)
  Speak_A1<-Salpha_A1+Sbeta_A1
  # fit A2 peak
  a2<-0.875+0.0268*(exp((lmax-665)/40.7))  
  A2<- 62.7 + 1.834*exp((lmax-625)/54.2)
  lmb2<-216.7+0.287*lmax
  p2<-317-1.149*lmax+0.00124*lmax*lmax  #Govardovskii calls this b but this confuses it with the alpha peak equation
  Salpha_A2<-1/((exp(A2*(a2-x))+exp(20.85*(0.9101-x))+exp(-10.37*(1.1123-x))+0.5343))
  Sbeta_A2<-0.26*exp(-((i-lmb1)/p1)^2)
  Speak_A2<-Salpha_A2+Sbeta_A2
  # weight by chromophore, sum and normalize
  Speaktot<-A1_chrom*Speak_A1+(100-A1_chrom)*Speak_A2
  Snorm<-Speaktot/max(Speaktot)
  Gov_results<-data.frame(Snorm)
  return(Gov_results)}
#================================================================================
#   FUNCTION: Absorptance_calc 
#      Calculate visual pigments absorptance which is more exact and doesn't 
#         small absorption assumption
#           Vispig = matrix of visual pigment absorbances (mu)
#           Lengths = photoreceptor lengths (l)
#         Absorptance = 1 - exp(-normalized absorption (waveln) * absorption peak * photoreceptor length)
#           see Land and Nilsson Animal Eyes box 3.2
#--------------------------------------------------------------------------------
Absorptance_calc<-function(Vispig, Abs_peak, Length) {
  AbsorptanceS<-1-exp(-Vispig$VPs*Abs_peak[1]*Length[1])
  AbsorptanceM<-1-exp(-Vispig$VPm*Abs_peak[2]*Length[2])
  AbsorptanceL<-1-exp(-Vispig$VPl*Abs_peak[3]*Length[3])
  Absorptance_results<-data.frame(Vispig$Wavelength, AbsorptanceS, AbsorptanceM, AbsorptanceL)
  colnames(Absorptance_results)<-c("Wavelength", "VPs", "VPm", "VPl")
  return(Absorptance_results) }
#
#==============================================================================
#  FUNCTION: mono_dataset
#     Calculate monochromatic color spectra (1 nm wide peaks across the spectrum)
#     The monochromatic colors essentially bound the color space of "normal" colors
#          lbeg, lend = wavelength range
#          monostart = first monochromatic color (e.g. 350 nm)
#          monostep = spacing between mono colors (e.g 5nm)
#          monoend = last monochromatic color (e.g. 700 nm)
# ------------------------------------------------------------------------------
Mono_dataset<-function(lbeg, lend, monostart, monostep, monoend) {
  mono1<-(monoend-monostart)/(monostep)+1
  wavetotal<-lend-lbeg+1
  waves<-matrix(lbeg:lend)
  peak<-matrix(nrow=mono1)
  mono_data<-data.frame(matrix(ncol=mono1, nrow=wavetotal))
  mono_data[is.na(mono_data)]<-0
  for (j in 1:mono1) {
    peak[j]<-monostart+(j-1)*monostep
    for (i in 1:wavetotal) {
      if (waves[i]==peak[j]) {
        mono_data[i,j]<-1
      }  
    } # i loop
  }   # j loop
  colnames(mono_data)<-peak
  rownames(mono_data)<-waves
  return(mono_data)   }
#
#==============================================================================
#  FUNCTION: Peak_data
#     Calculate ONE Gaussian peaked color spectra (variable peak and width)
#          lbeg, lend = wavelength range
#          peak = peaked color (e.g. 400 nm)
#          peaksd = standard deviation of peak
#          peakoffset = offset from baseline
# ------------------------------------------------------------------------------
Peak_data<-function(lbeg, lend, peak, peaksd, peakoffset) {
  wavetotal<-lend-lbeg+1
  waves<-matrix(lbeg:lend)
  peak_data<-matrix(nrow=wavetotal)
  for (i in 1:wavetotal) {
    coeff<-sqrt(2*pi)*peaksd*(1-peakoffset)+peakoffset
    peak_data[i]<-coeff*dnorm(waves[i], peak, peaksd)+peakoffset
  }   # i loop
  rownames(peak_data)<-waves
  return(peak_data)   }
#==============================================================================
#  FUNCTION: Step_data
#     Calculate ONE Gaussian stepped color (variable peak and rise rate)
#          lbeg, lend = wavelength range
#          peak = first stepped color (e.g. 400 nm)
#          peaksd = standard deviation of peak
#          peakoffset = offset from baseline
# ------------------------------------------------------------------------------
Step_data<-function(lbeg, lend, peak, peaksd, peakoffset) {
  wavetotal<-lend-lbeg+1
  waves<-matrix(lbeg:lend)
  peak_data<-matrix(nrow=wavetotal)
  for (i in 1:wavetotal) {
    coeff<-sqrt(2*pi)*peaksd*(1-peakoffset)+peakoffset
    peak_data[i]<-coeff*dnorm(waves[i], peak, peaksd)+peakoffset
    if (waves[i]>peak)  {
      peak_data[i]<-1
    }
  }   # i loop
  rownames(peak_data)<-waves
  return(peak_data)   }
#  
#==============================================================================
#   FUNCTION:  Qcatch_calc  
#     Calculate the stimulation of the S, M and L cones in terms of quantum catch
#         light = the illuminant
#         fishlens = the lens transmission
#         Vispig = 3 visual pigments which can be coexpresing and involve A1 or A2
#         targets = the set of color targets
#     This version works on target file where there is no wavelength column
# ----------------------------------------------------------------------
Qcatch_calc<-function(light, fishlens, Vispig, targets) {
  target_labels<-colnames(targets)
  number_colors<-length(target_labels)   # number of colors in set
  #
  # Calculate the von Kries correction for accommodation to the illuminant
  vonK_illum_QCs=0
  vonK_illum_QCm=0
  vonK_illum_QCl=0
  vonK_illum_QCs<-sum(light$IRRAD*fishlens$Transmission*Vispig$VPs)
  vonK_illum_QCm<-sum(light$IRRAD*fishlens$Transmission*Vispig$VPm)
  vonK_illum_QCl<-sum(light$IRRAD*fishlens$Transmission*Vispig$VPl)
  qcatchs=NULL
  qcatchm=NULL
  qcatchl=NULL
  for (j in 1:number_colors) {
    qcatchs[j]<-sum(light$IRRAD*targets[,j]*fishlens$Transmission*Vispig$VPs)/vonK_illum_QCs
    qcatchm[j]<-sum(light$IRRAD*targets[,j]*fishlens$Transmission*Vispig$VPm)/vonK_illum_QCm
    qcatchl[j]<-sum(light$IRRAD*targets[,j]*fishlens$Transmission*Vispig$VPl)/vonK_illum_QCl
  }
  Qcatchdata_results<-data.frame(qcatchs,qcatchm,qcatchl)
  return(Qcatchdata_results)  }
#==============================================================================
#  FUNCTION : XYplotS        XY plot of quantum catch data
#       This will convert quantum catch data into normalized S, M, L data for plotting
#       in an SML triangular color space.  To do this, each datapoint is converted into X and Y for plotting.
#
#       This version has S on top, M lower left and L lower right
#       White point is at (0,0)
# ----------------------------------------------------------------------
XY_plotS<-function(QCdataset) {   
  S=NULL
  M=NULL
  L=NULL
  Xcolor=NULL
  Ycolor=NULL
  number_colors<-nrow(QCdataset)
  for (j in 1:number_colors) {
    S[j]<-QCdataset$qcatchs[j]/(QCdataset$qcatchs[j]+QCdataset$qcatchm[j]+QCdataset$qcatchl[j])
    M[j]<-QCdataset$qcatchm[j]/(QCdataset$qcatchs[j]+QCdataset$qcatchm[j]+QCdataset$qcatchl[j])
    L[j]<-QCdataset$qcatchl[j]/(QCdataset$qcatchs[j]+QCdataset$qcatchm[j]+QCdataset$qcatchl[j])
    Xcolor[j]=(L[j]-M[j])/sqrt(3)
    Ycolor[j]=S[j]-1/3
  }
  XYplotS_results<-data.frame(Xcolor,Ycolor)
  return(XYplotS_results)   }
#
#==============================================================================
#  FUNCTION : JND_calc   
#     Calculate JND values between two sets of targets 
#       based on Vorobyev and Osario 1998
#
# ----------------------------------------------------------------------
JND_calc<-function(QCset1, colorlabel1, QCset2, colorlabel2, Webers)  {
  number_colors1<-length(colorlabel1)
  number_colors2<-length(colorlabel2)
  JND_results=data.frame(matrix(nrow=number_colors1, ncol=number_colors2))
  ws<-Webers[1]
  wm<-Webers[2]
  wl<-Webers[3]
  denom<-(ws*wm)^2+(ws*wl)^2+(wm*wl)^2
  for (i in 1:number_colors1) {
    for (j in 1:number_colors2) {
      Dfs<-log(QCset1[i,1]/QCset2[j,1])
      Dfm<-log(QCset1[i,2]/QCset2[j,2])
      Dfl<-log(QCset1[i,3]/QCset2[j,3])
      JND_results[i,j]<-sqrt((ws^2*(Dfm-Dfl)^2+wm^2*(Dfs-Dfl)^2+wl^2*(Dfs-Dfm)^2)/denom)
    } } 
  rownames(JND_results)<-colorlabel1
  colnames(JND_results)<-colorlabel2
  return(JND_results)  }
#
#==============================================================================