#Home Range

library(tidyverse)
library(ctmm)
data("buffalo") #these data are a `ctmm` telemetry object
cilla_c <- buffalo$Cilla
cilla_df <- as(cilla_c, "data.frame")
head(cilla_df)
plot(cilla_c)

#Transform into different objects
library(move)
library(sf)
library(adeHabitatHR)
cilla_m <- move(x = cilla_df$x,
                y = cilla_df$y,
                time = cilla_df$timestamp,
                proj = CRS(cilla_c@info$projection), #I took the CRS argument from cilla since it was a bit odd
                animal="cilla", 
                sensor="GPS", 
                data = cilla_df)
cilla_spdf <- as(cilla_m, "SpatialPointsDataFrame") #To a SpatialPointsDataFrame. The other way to do this is to:
cilla_spdf2 <- cilla_df #duplicate the data frame
coordinates(cilla_spdf2) <- ~x+y
proj4string(cilla_spdf2) <- CRS(cilla_c@info$projection)
cilla_sf <- st_as_sf(cilla_df, 
                     coords = c("x", "y"),
                     crs = CRS(cilla_c@info$projection))


# Minimum Convex Polygon
# The simplest delineation of a home range is an MCP – creating the polygon of minimum area around a certain percentage of relocation points. The MCP is simple and used widely in ecology.
# 
# If you are curious to estimate the overall area of your animal’s home range the move package includes a function to bootstrap the mcp estimation:
  
hrBootstrap(x=cilla_m, rep=25, unin='km', unout='km2')

library(adehabitatHR)
library(sf)
mcp <- mcp(cilla_spdf2, percent=90)

ggplot() + 
  geom_sf(data = st_as_sf(mcp)) + 
  geom_sf(data=st_as_sf(cilla_spdf2))+
  theme_bw()

# We’ll get to better metrics shortly but if you want to compare the area of your mcp across percentages, the mcp.area function works well:
  
mcp.area(cilla_spdf2, percent = seq(20,100, by = 5),
           unin = c("m", "km"),
           unout = c("ha", "km2", "m2"), plotit = TRUE)

# If you are curious to see what’s going on under the hood of the adehabitatHR  mcp functions, I recommend checking out this blog post on the subject by Mitchell Gritts.
# 
# https://www.r-bloggers.com/home-range-estimation-mcp



# Kernel Density Estimation
# Worton Kernel UD
# The “classic” utilization distribution: Worton (1995)
# 
# > The Utilization Distribution (UD) is the bivariate function giving the probability density that an animal is found at a point according to its geographical coordinates. Using this model, one can define the home range as the minimum area in which an animal has some specified probability of being located.
# 
# Let’s go back to our girl Cilla to work some of this out.

kud <- kernelUD(cilla_spdf2)  # h = href is the default - ad hoc method for determining h
## Warning in kernelUD(cilla_spdf2): xy should contain only one column (the id of the animals)
## id ignored
image(kud) + title("Cilla the Buffalo's UD")


#We can get the area of this
jj <- kernel.area(kud) ## home range size from the kernel Utilization Distribution
jj
plot(jj, xaxt = "n", xlab = "Kernel UD", ylab = "Area (ha)")
axis(1, at=1:16, labels=seq(20,95, by = 5))## Plots home range size

#Often we also want to plot the contours of different percentages of the home range along with the UD

ver95 <- getverticeshr(kud, percent = 95) ## get the 95% home-range contours
ver80  <- getverticeshr(kud, percent = 80) #get the 80th percentile

sp::plot(ver95)+ #plot the full 95% Contour
  sp::plot(ver80, add=TRUE, col="green")+ ## Plots contours
  points(cilla_spdf2)

# Autocorrelated Gaussian reference function Kernel Density Estimate (aKDE)
#https://cran.r-project.org/web/packages/ctmm/index.html

#The previous home range estimates are pretty straight forward & the methods to get numbers out are canned. The aKDE however needs a little bit of extra attention to make sure that the correct model is chosen to get adequate representation of the autocorrelation structure in the data. 

# As we read, autocorrelation can make you really underestimate the true home range of the animals you are studying. One method around this is the aKDE that separates out the movement trajectory from the home range calculation in a way that gives stable and robust estimates of true home ranges (it seems).
# 
# We’ll calculate these in ctmm. Fitting an aKDE is a little more involved, but is worth the effort to get better estimates. If you type in browseVignettes(package = "ctmm") in your console, you’ll see a few different tutorials that help you with this process. We’ll go through some of it now, and start with the move object since the ctmm authors suggest taking your data directly from Movebank and importing with with as.telemetry(). Our object cilla is already a telemetry object (the kind required by ctmm), so we’ll look at that and move on.

str(cilla_c)
plot(cilla_c)

# Variograms
# https://cran.r-project.org/web/packages/ctmm/vignettes/variogram.html

# Variograms are an unbiased way to visualize autocorrelation structure when migration, range shifting, drift, or other translations of the mean location are not happening. When drift occurs in the data, then the variogram represents a mixture of both the drift and the autocorrelation structure, each of which contains distinct movement behaviors.

# A variogram can give us an indiciation of how long it takes for our data to start reflecting true home range. This will be based on the range crossing time of the animal.

#We can get a guess at a variogram by fitting one to start with & then modifying it.
SVF <- variogram(cilla_c)
level <- c(0.5,0.95) # 50% and 95% CIs
xlim <- c(0,12 %#% "hour") # 0-12 hour window

par(mfrow=c(1,2)) #this is a graphics command that says I want 1 row, 2 columns
plot(SVF,xlim=xlim,level=level)
title("zoomed in")
plot(SVF,fraction=0.65,level=level)
title("zoomed out")

# There are many caveats to fitting a variogram and a few different methods to do so based on the distibution of the data (IID (independent, identical distribution), OUF Ornstein-Uhlenbeck foraging movement). It’s important to have a good variogram estimated for the data since this will be used later in calculating the range. Check out the variogram vignette when you type in browseVignettes(package = "ctmm"). We’re just going to jump to fitting a variogram the easy way. This will let us adjust the parameters to get an optimal fit to the data and then save those into an object called GUESS

variogram.fit(SVF)

#You can use the sliders to optimize the fit of the variogram model, then save it as GUESS

#We can plot how well that GUESS fits vs the preliminary variogram

plot(SVF, CTMM=GUESS, level=level, col.CTMM="purple", xlim=xlim) #let's see how well GUESS fits the variogram


#How well do that work? And what model do we choose?
# If you have a complex, hypothetical model in mind, say the OUF m.ouf as you would get from variogram fitting the easy way, then you can perform model selection more conveniently with the ctmm.select function. ctmm.select considers the initial guess (hypothesis) and then iterates this model to select the best model based upon an information criteria.

FITZ <- ctmm.select(cilla_c, GUESS, verbose=TRUE, cores=2) #This will take a little time
summary(FITZ)

# The isotropic and anisotropic (isotropic=FALSE) flags correspond to circular and elliptical covariances respectively—an option we did not consider above. The OUf model is a special case of the OUF model where the two autocorrelation timescales, τ position and τ velocity, cannot be distinguished. This model is usually only relevant for short tracks of data. The IID model was never considered here by ctmm.select because it first requires selecting OU over OUF in the nested model hierarchy. See help("ctmm") for more options.

#This tells us that the OUF anisotropic model is best.



# Variograms as a Diagnostic for Maximum Likelihood
# Now its time to make sure that our selected model is explaining the most significant features of the animal’s movement. Let us plot the variogram again with our fit models

par(mfrow=c(1,2))
plot(SVF, CTMM=FITZ, col.CTMM=c("red", "purple", "blue", "green"), fraction=0.65, level=0.5)
title("zoomed out")
plot(SVF, CTMM=FITZ, col.CTMM=c("red", "purple", "blue", "green"), xlim=xlim, level=0.5)
title("zoomed in")

#Notice that the blue OU model is significantly biased downward and is underestimating diffusion. This is because the continuous-velocity behavior at short time lags, which the OU model does not account for, is throwing off the estimate. The green OUf model where timelags are idential (removing autocorrleation) underestimates diffusion altogether.

# While the red OUF is the selected model among all candidates, it looks slightly biased upwards in comparison to the variogram fit m.ouf. Some of this is due to sampling variability, which be (partially) remedied in the variogram method by increasing the res argument. Differences here can also arise from not accounting for telemetry error, which is unfortunately not annotated in this data.

# Fitting an aKDE
# Now that we have an estimated variogram for our buffalo, we can actually fit an estimate of her home range. The ctmm.fit functions can take a lot of time.

M.IID <- ctmm.fit(cilla_c) # no autocorrelation timescales. By not specifying the variogram guess, ctmm.fit creates an IID Kernel Density Estimate corrected for small sample sizes (KDEc)
m.ouf <- ctmm.guess(cilla_c, interactive=FALSE) # automated model guess
M.OUF <- ctmm.fit(cilla_c, m.ouf) #our best model fit.

# M.IID is the inappropriate, IID model, which will result in a conventional kernel-density estimate, while M.OUF is the better, continuous-velocity OUF model. Note that you want the best model for each individual, even if that differs by individual. Different movement behaviors and sampling schedules will reveal different autocorrelation structures in the data.
# 
# Now we can calculate an akde object for each model. These also can take some time to run.

UD0 <- akde(cilla_c, M.IID)
UD2 <- akde(cilla_c, M.OUF)
UD2w <- akde(cilla_c, M.OUF, weights=TRUE)

#Finally we calculate UDs with and without accounting for autocorrelation (M.OUF versus M.IID), with and without optimal weighting of the data (weights=TRUE). Now let us plot the results.

# calculate one extent for all UDs to make plotting look consistent
EXT <- extent(list(UD0,UD2,UD2w),level=0.95)

par(mfrow=c(2,2))
plot(cilla_c, UD=UD0, xlim=EXT$x, ylim=EXT$y)
title(expression("IID KDE"["C"]))
plot(cilla_c, UD=UD2, xlim=EXT$x, ylim=EXT$y)
title(expression("OUF AKDE"["C"]))
plot(cilla_c, UD=UD2w, xlim=EXT$x, ylim=EXT$y)
title(expression("weighted OUF AKDE"["C"]))


# By default both the density function and its 95% contours are plotted along with the location data. The middle contour represent the maximum likelihood area where the animal spends 95% of its time. This percentage can be changed with the level.UD option (see help(plot.telemetry)). The inner and outer contours correspond to confidence intervals on the magnitude of the area, which can be adjusted with the level option.
# 
# The optimal bandwidth determines the “resolution” of the kernel density estimate. By default we plot grid lines with dimensions matching the standard deviations of the individual kernels. This gives a rough guideline as to what spatial details are and are not important in the density estimate. One can see that the IID KDEC estimate fits tightly to the data and reports many significant details in the buffalo’s home range. The autocorrelated estimates predict future space use more accurately, based on the diffusion model, and yield a more honest account of their uncertainties.
# 
# Finally, we can compare the area estimates and effective sample sizes.
#DOF = Degrees of Freedom.

summary(UD0, units = TRUE) #IID
summary(UD2, units = TRUE) #M.OUF
summary(UD2w, units = TRUE) #M.OUF weighted AKDE



