######################################
#
#     VisualModel_Acalliptera.R    Trichromat
#       
#####################################
# Karen Carleton March 2022
#
# This program will calculate JNDs between a set of Gaussian colors (including black and white)
# for a trichromatic visual system containing S, M and L cones.  
# It compares two different Massoko illuminants (shallow and deep).  
# It then calculates the quantum catches for the three receptors and uses
# Receptor Noise Limited model to determine the Just Noticeable Differences (JNDs)
# between those colors to determine which color combos are most discriminable
#  see Vorobyev and Osorio 1998 for RNL and Champ et al 2016
# Functions are contained in the separate file AquaVis_Acalliptera.R and should be loaded first
#   
# Inputs:
#   Light spectra for Massoko shallow and deep, binned to 1 nm resolution with Binning_spectra.R
#   lmax determined from other cichlid species, including A burtoni (see Opsin_align.xlsx)
#   opsin gene expression from Massoko_opsin_expression_data averaging shallow and deep individuals
#      For coexpression, we assume single cones coexpress all SWS genes - weighted sum
#                double cone#1 coexpresses all RH2 genes - weighted sum
#                double cone#2 expresses LWS - this is a bit diff from Malawi (Dalton et al 2014)
#   Gaussian colors are calculated as peaks (425-560) and steps (560-625) plus broad black and white "peak"
#    
# ===================================================================================
# Start of program - INSERT DATA HERE
# ===================================================================================
#
# Load functions first: 
data_folder<-"~/Programs/Rwork/Color/Acalliptera_opsin/Final4github/"
Aquavis_functions<-paste(data_folder,"AquaVis_Acalliptera.R",sep="")
source(Aquavis_functions) # load functions
#
# 1. Choose data files to use: illuminants
illum_file1<-paste(data_folder,"Sidewelling_Sept2022_shallow_NormAvg2-10m_irradiance_binned_350.txt",sep="")
illum_file2<-paste(data_folder,"Sidewelling_Sept2022_deep_NormAvg20-24m_irradiance_binned_350.txt",sep="")
illum_name1<-"Massoko_Shallow"
illum_name2<-"Massoko_Deep"

# Note illumination files will set the wavelength range for the visual pigments
# These are sidewelling irradiances which illuminate fish colors
  illum1<-read.table(illum_file1, header=TRUE)
  colnames(illum1)<-c("Wavelength","IRRAD")
  illum2<-read.table(illum_file2, header=TRUE)
  colnames(illum2)<-c("Wavelength","IRRAD")
  number_wave<-length(illum1$Wavelength)
  lambda_begin<-illum1$Wavelength[1] # wavelength set here
  lambda_end<-illum1$Wavelength[number_wave]
  waves<-matrix(lambda_begin:lambda_end)
  #
# Choose colors and legend labels for the plots of colors in SML triangle
  legend1<-illum_name1
  legend2<-illum_name2
  legend_mono<-"Mono"
  illum_name<-"Massoko"
  color1<-"blue"
  color2<-"green"
  color_mono<-"gray"
#
# 2. Choose a set of evenly spaced colors to compare
  stdev_peak<-55     # 55 for peaks    
  stdev_step<-30     # 30 for steps    
  baseln_peak<-0.1   # 0.1 for peaks
  baseln_step<-0.1  # 0.1 for steps 
  peak_start<-375    # peak colors are shorter wavelength (violet to green)
  peak_end<-549
  step_start<-550    # step colors are longer wavelength (yellow to red)
  step_end<-625
  color_step<-1
  peak_BW<-500      # black and white are broad
  stdev_BW<-350
  coeffB<-0.1
  plot_color<-c("purple","blueviolet","blue","aquamarine","chartreuse","green3","greenyellow","yellow","orange","red","grey70","black")
  color_name<-"Gaussian colors"
#  
# 3. Choose the visual system
# 3a. Set visual pigments and chromophore - taken from A burtoni and other cichlid expression values based on A calliptera sequences
  lmax1 <- 378   # matches P acei
  lmax2 <- 420   # average other cichlids
  lmax3 <- 455   # matches A burtoni
  lmax4 <- 481   # average other cichlids
  lmax5 <- 512   # average other cichlids
  lmax6 <- 523   # matches A burtoni
  lmax7 <- 562   # matches A burtoni
  A1 <- 100  #  % A1 chromophore
  VP_color<-c("blue","green","red")
#  
# 3b. Set opsin %coexpression. Assume SWS genes coexpressed in single cones
#     RH2 genes coexpressed in one double cone and LWS in other double cone
#     Use violet lens from M auratus
# Select vision
  visual_palette<-"long" # choose from "shallow", "deep", or "long"
# If shallow vision use this coexpression
  if (visual_palette=="shallow") {
    SWS1exp<-7.98; SWS2Bexp<-72.56; SWS2Aexp<-19.46
    RH2Bexp<-3.56; RH2Abexp<-7.70; RH2Aaexp<-40.15
    LWSexp<-48.59
  }
# If deep vision use this coexpression  
  if (visual_palette=="deep") {
    SWS1exp<-2.1; SWS2Bexp<-32.7; SWS2Aexp<-65.19
    RH2Bexp<-3.12; RH2Abexp<-11.80; RH2Aaexp<-48.04
    LWSexp<-37.04
  }
# If pure visual pigments with long palette, i.e. no coexpression
  if (visual_palette=="long") {
    SWS1exp<-0; SWS2Bexp<-0; SWS2Aexp<-100
    RH2Bexp<-0; RH2Abexp<-0; RH2Aaexp<-100
    LWSexp<-100
  }
#
# 4. Sensitivities as set by weber fraction and cone ratios
# For JNDs need the Weber fraction or the ratio of intensities that is just 
# discriminable at threshold. We'll use Weber = 0.1 for birds but could increase based on Escobar-Camacho 2019 for cichlids
# Actual value doesn't matter as we are making relative comparisons  
  WeberL<-0.1
# The Weber fraction for each cone type is determined by the absolution Weber fraction WeberL
# and the number ratio of S, M and L cone types in the retinal mosaic.  
# Cichlids have S : M : L of 1: 2 : 2  where wi = wL sqrt(nL/ni) where L is long cone 
# 
  ws<-WeberL*sqrt(2)
  wm<-WeberL
  wl<-WeberL
  W<-c(ws, wm, wl)

# 5. Other parameters
  peak_abss<-0.009  # the peak absorbance determined by MSP in units of per um Carleton 2000
  peak_absm<-0.015
  peak_absl<-0.015
  peak_abs<-c(peak_abss,peak_absm,peak_absl)
  Ls<-5.5   # single cone length um estimated from MSP photographs
  Lm<-26    # double cone length um
  Ll<-Lm
  Length<-c(Ls, Lm, Ll)
  method<-"Abs"     #Choose either "Gov" or "Abs" for Govardovskii or absorptance. Abs is most accuriate but not that different from Gov
  lens_filename<-"Malawi_violetlens.txt"   # Mauratus_lens1048
#  
# 6. Choose the output files for JND
  directory<-"~/Programs/Rwork/Color/Acalliptera_opsin/Final4github/" # can be different from data folder but for now they are the same
  text1<-"Acalliptera_VisModel_"    # prefix for output file names
  text2<-paste(visual_palette,"VP",sep="")
  out_file1<-paste(directory,text1,text2,legend1,".txt",sep="") # add visual system and illuminant to file name
  out_file2<-paste(directory,text1,text2,legend2,".txt", sep="")
  plot_file<-paste(directory,text1,text2,"Massoko",sep="")
  
  save2file<-"Y"  # whether or not to save JND values to file Y or N
  
# 7. Stuff for plotting
  Rymax<-1
  Rxmin<-350
  Rxmax<-750
  Lymax<-1
  color4<-"black"

# For "Best" plot, need to tell it which column of reflect has best colors, e.g. peak 402 (reflect [,28]) and step 572 (reflect [,198])
    bestpeakNo<-28
    beststepNo<-198

#  Choose type of plot
  which_plot<-"Best"   
 
# Plot options 
#  N = none
#  R = reflectances Note: if 1nm spacing can't really see them; increase color_step to perhaps 25nm
#  V = visual pigments, could add lens transmission to this plot as do in Best plot
#  T = SML triangle with colors plotted in this visual space
#  RV = Reflectance and VP   again don't use too small of a color_step 
#  VL = Vispig + illum  
#  Best = vision & best colours  For Best plot, need to tell it which columns of reflect matrix are the best peak and step
#  plots visual pigments * lens and then illuminant and then illuminant * best peak and step
#
# For SML triangle plot (T) need to set end and spacing of monochromatic colors that bound space. Monostart is set by visual-palette.
  monostep<-5
  monoend<-700
#
# -----------------Done with data entry ---------------------------------------
# Let the calculations begin
# -------------------------------------------------------------------------
# Calculate visual pigments
# -------------------------------------------------------------------------
# Calculate the absorbances for pure pigments
  SVP1<-Opsin_gov(lambda_begin, lambda_end, lmax1,A1)
  SVP2<-Opsin_gov(lambda_begin, lambda_end, lmax2,A1)
  SVP3<-Opsin_gov(lambda_begin, lambda_end, lmax3,A1)
  DVP4<-Opsin_gov(lambda_begin, lambda_end, lmax4,A1)
  DVP5<-Opsin_gov(lambda_begin, lambda_end, lmax5,A1)
  DVP6<-Opsin_gov(lambda_begin, lambda_end, lmax6,A1)
  DVP7<-Opsin_gov(lambda_begin, lambda_end, lmax7,A1)
#
# Calculate the coexpressed pigments as weighted combinations of pure pigments
# store the s, m and l visual pigments in data.frame VP
    VPs<-(SVP1*SWS1exp+SVP2*SWS2Bexp+SVP3*SWS2Aexp)/100
    VPm<-(DVP4*RH2Bexp+DVP5*RH2Abexp+DVP6*RH2Aaexp)/(RH2Bexp+RH2Abexp+RH2Aaexp)
    VPl<-DVP7
    monostart<-400   # 400 for medium or long palettes
  VP<-data.frame(waves,VPs$Snorm,VPm$Snorm,VPl$Snorm)
  colnames(VP)<-c("Wavelength","VPs","VPm","VPl")  #set the column names for calculations below
#
# If using Absorptance, calculate it here based on peak absorption and photoreceptor lengths
  if (method=="Abs") {
     Absorptance<-Absorptance_calc(VP, peak_abs, Length)}
#
# Read in the lens transmission and set column names.
#
  lens_file=paste(data_folder,lens_filename,sep="")
  lens_trans<-read.table(lens_file, header=TRUE)
  colnames(lens_trans)<-c("Wavelength","Transmission")
# -------------------------------------------------------------------------
# Calculate colors - Gaussian peaks and steps
# -------------------------------------------------------------------------
  peak_number<-(peak_end-peak_start)/color_step+1
  step_number<-(step_end-step_start)/color_step+1
  color_number<-peak_number+step_number+2
  waveln<-lambda_end-lambda_begin+1
  reflect_name<-matrix(nrow=color_number)
  reflect<-matrix(nrow=waveln, ncol=color_number)
  for (i in 1:peak_number) {
    peak1<-peak_start+(i-1)*color_step
    reflect[,i]<-as.matrix(Peak_data(lambda_begin, lambda_end, peak1, stdev_peak, baseln_peak))
    reflect_name[i]<-paste("P",peak1,sep="")
  }
  for (i in 1:step_number) {
    step1<-step_start+(i-1)*color_step
    reflect[,peak_number+i]<-Step_data(lambda_begin, lambda_end, step1, stdev_step, baseln_step)
    reflect_name[peak_number+i]<-paste("S",step1,sep="")
  }
  # add white and black
  reflect[,color_number-1]<-Peak_data(lambda_begin, lambda_end, peak_BW, stdev_BW, baseln_peak)
  reflect_name[color_number-1]<-"White"
  reflect[,color_number]<-coeffB*Peak_data(lambda_begin, lambda_end, peak_BW, stdev_BW, 0)
  reflect_name[color_number]<-"Black"
  colnames(reflect)<-reflect_name
#
# Calculate quantum catch for colors with two illuminants, shallow and deep
#
  if (method=="Gov") {
    Qcatchdata1<-Qcatch_calc(illum1, lens_trans, VP, reflect)
    Qcatchdata2<-Qcatch_calc(illum2, lens_trans, VP, reflect)
    }
    
  if (method=="Abs") {
    Qcatchdata1<-Qcatch_calc(illum1, lens_trans, Absorptance, reflect)
    Qcatchdata2<-Qcatch_calc(illum2, lens_trans, Absorptance, reflect)
    }

# For SML trichromatic visual space, calculate color (X,Y) positions    
  Plotdata1<-XY_plotS(Qcatchdata1)
  number_colors1<-nrow(Qcatchdata1)         # number of colors being analyzed
#
  Plotdata2<-XY_plotS(Qcatchdata2)
  number_colors2<-nrow(Qcatchdata2)         # number of colors being analyzed
#
# Calculate monochromatic colors and their (X,Y) positions (can plot as a line)
#
  reflect_mono<-Mono_dataset(lambda_begin, lambda_end, monostart, monostep, monoend)
  if (method=="Gov") {
    Qcatch_mono<-Qcatch_calc(illum1, lens_trans, VP, reflect_mono)}
  if (method=="Abs") {
    Qcatch_mono<-Qcatch_calc(illum1, lens_trans, Absorptance, reflect_mono)}
  Plotmono<-XY_plotS(Qcatch_mono)
#
# -------------------------------------------------------------------------
#  Calculate JNDs Here we calculate JND between each pair of colors
#
# -------------------------------------------------------------------------
# Do this for both illuminants - are doing all possible pairwise colors 
  JND1<-JND_calc(Qcatchdata1, reflect_name, Qcatchdata1, reflect_name, W)
  JND2<-JND_calc(Qcatchdata2, reflect_name, Qcatchdata2, reflect_name, W)
#
# Output JND results for each illuminant to a file
  if (save2file=="Y") {
     write.table(JND1,file=out_file1, sep="\t", quote= FALSE)
     write.table(JND2,file=out_file2, sep="\t", quote= FALSE)
  }
  
# Make whichever plot you selected above with which_plot 
  if (which_plot=="R") {
    par(xpd=TRUE)
    par(mar=c(4,4,2,4))  # set margins to reduce white space
    plot(waves, reflect[,1], type="n", axes=TRUE, xlab="Wavelength", ylab="Reflectance", xlim=c(Rxmin, Rxmax), ylim=c(0,Rymax), cex=0.7)
    for (k in 1:color_number) {
      lines(waves, reflect[,k], type="l", col=plot_color[k])
#      lines(waves, reflect[,k], type="l", col=color4)
    }
  } # of which_plot for reflectances

  # visual pigment plot
  if (which_plot=="V") {
    par(xpd=TRUE)
    par(mar=c(4,4,2,4))  # set margins to reduce white space
    plot(VP[,1], VP[,2], type="n", axes=TRUE, xlab="Wavelength", ylab="Absorption", xlim=c(Rxmin, Rxmax), ylim=c(0,1), cex=0.7)
    for (k in 2:4) {
      lines(VP[,1], VP[,k], type="l", col=VP_color[k-1])
    }
  }  # of which_plot for visual pigments
  
  if (which_plot=="RV") { # Reflective colors and Visual pigments
    #  par(xpd=TRUE)
    par(mar=c(4,4,2,4))  # set margins to reduce white space
    plot(VP[,1], VP[,2], type="n", axes=TRUE, xlim=c(Rxmin, Rxmax), ylim=c(0,1), xlab="", ylab="", cex=0.7)
    axis(2, ylim=c(0,1))
    mtext(2, text="Absorption or Reflectance", col=color1)
    for (k in 2:4) {
      lines(VP[,1], VP[,k], type="l", lty=2, col=VP_color[k-1]) # plot the visual pigments
      axis(side=2, ylim=c(0,1))
    }
    par(new=T)
    plot(waves, reflect[,1], type="n", axes=FALSE, xlab="Wavelength", xlim=c(Rxmin, Rxmax), ylab="", ylim=c(0,Rymax), cex=0.7)
    for (m in 1:color_number) {
        lines(waves, reflect[,m], type="l", col=plot_color[m])
    }
    axis(side=4,  ylim=c(0,Rymax))
    mtext(side=4, text="Reflection", col=color1)
  }  # of which_plot RV
  
  if (which_plot=="VL") { # Visual pigments and Light / illuminants
    par(mar=c(4,4,2,4))  # set margins to reduce white space
    plot(VP[,1], VP[,2], type="n", axes=TRUE, xlim=c(Rxmin, Rxmax), ylim=c(0,1), xlab="", ylab="", cex=0.7)
    axis(2, ylim=c(0,1))
    mtext(2, text="Absorption", col=color1)
    for (k in 2:4) {        # plot visual pigments
      lines(VP[,1], VP[,k], type="l", col=VP_color[k-1])
      axis(side=2, ylim=c(0,1))
    }
    par(new=T)
    if (legend1=="Massoko_Shallow") {
      lines(illum1[,1], illum1[,2], type="l", col=color1, lty=3)
    }
    if (legend2=="Massoko_Deep") {
      lines(illum1[,1], illum2[,2], type="l", col=color2, lty=3)
    }
    axis(side=4, ylim=c(0,1))
    mtext(side=4, text="Side IRRAD", col=color4)
  }  # of which_plot VL
  
  if (which_plot=="T") {# Plot the data in SML visual space, including the triangular axes. Here the S cone is on top
    maintitle1<-paste(color_name,"comparing", illum_name, "illuminants",sep=" ")
    maintitle2<-paste(visual_palette,"vision", A1,"% A1", sep=" ")
    main_text<-paste(maintitle1,"\n",maintitle2)
    legend_text<-c(legend1, legend2, legend_mono)
    legend_color<-c(color1,color2, color_mono)
    triangle<-matrix(c(-1/sqrt(3),-1/3,1/sqrt(3),-1/3,0,2/3,-1/sqrt(3),-1/3), nrow=4, ncol=2, byrow="TRUE")
    if (save2file=="Y") {
      plot_name<-paste(plot_file,".tiff",sep="")
      png(filename=plot_name, width=1024, height=686)}
  # otherwise plot on screen
    par(mar=c(0,0,2,0))  # set margins to reduce white space
    plot(triangle, type="n", asp=1, xlim=c(-0.6,0.6), ylim=c(-0.4,0.7), axes= FALSE)
    points(triangle, type="l")
    points(Plotdata1, col=color1, pch=20, cex=1)
    points(Plotdata2, col=color2, pch=20, cex=1)
    points(Plotmono, col=color_mono, type="l", lty=2)  # plot mono as dashed line
    legend("topright",pch=c(20,20), col=legend_color, legend = legend_text, cex = 0.7)
    title(main = main_text, cex.main=1)
  }  # of which_plot T 
 
  if (which_plot=="Best") { # vision and best colors - NEED TO specify colors by HAND HERE
    par(mar=c(4,4,2,4))  # set margins to reduce white space
    plot(VP[,1], VP[,2], type="n", axes=TRUE, xlim=c(Rxmin, Rxmax), ylim=c(0,1), xlab="", ylab="Absorption or Reflectance", cex=0.7)
    axis(2, ylim=c(0,1))
#    mtext(2, text="Absorption or Reflectance", col=color1)
    for (k in 2:4) {
      lines(VP[,1], VP[,k]*lens_trans[,2], type="l", col=VP_color[k-1])
      axis(side=2, ylim=c(0,1))
    }
    par(new=T)
    plot(waves, reflect[,1], type="n", axes=FALSE, xlab="Wavelength", xlim=c(Rxmin, Rxmax), ylab="", ylim=c(0,1), cex=0.7)
    lines(waves, reflect[,bestpeakNo]*illum1[,2], type="l", col="cyan", lty=2, lwd=2)
    lines(waves, reflect[,beststepNo]*illum1[,2], type="l", col="orange", lty=2,lwd=2)
    lines(waves, illum1[,2], type="l", col="grey50", lty=3, lwd=2)
#    axis(side=4,  ylim=c(0,Rymax))
#    mtext(side=4, text="Reflection", col=color1)
  }  # of which_plot Best
    
  
