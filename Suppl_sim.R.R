#####################################################################
# In this program we generate bioclim environmental variables to 
# simulate the species distributions.
#####################################################################

# Load required Packages
library(rgdal)
library(raster)
library(gstat)

setwd("...")  # set working directory
########################################################################################
# 	Function to simulate a spatially correlated variable using a variogram.
#	x.range	=	the range of long values to simulate the varible
#	y.range =	the range of lat values to simulate the varible
#	res		=	the size of the cell, can be a single value or a vector defining the 
#				width and height of the cell
#	ext		=	extent of the raster. This can be used instead of x.range and y.range
#	psill	=	Limit of the variogram
#	model	=	model to use to simulate the variable. see vgm() for models. 
#	range	=	Distance in which the difference of the variogram from the sill becomes negligible.
#				(definition from wikipedia)
#	nugget	=	the height of the jump of the semivariogram at the discontinuity at the origin.
#				(definition from wikipedia)
#	nsim	=	number of varibles to simulate using the same parameters
#	beta	= 	the expected value of the variable
#	nmax	=	the number of nearest observations that should be used for a simulation. Nearest is defined in terms
#				of the spatial locations. See function gstat().

rspatvar	<-	function(formula,locations,x.range,y.range,ext=NULL,res,psill=NA,model="Sph",range=NA,nugget=NULL,nsim=1,beta,nmax
                     ,proj4string=CRS(as.character(NA))){
  
  if(!is.null(ext)){x.range<-ext[1:2];y.range<-ext[3:4]}
  
  res.len			<-	length(res)
  if(res.len==1){res	<-	c(res,res)}
  x.seq	<-	seq(from=x.range[1],to=x.range[2],by=res[1])
  y.seq	<-	seq(from=y.range[1],to=y.range[2],by=res[2])	
  
  xy		<-	expand.grid(x=x.seq,y=y.seq)
  
  if(!is.null(nugget)){vgm.mod		<-	vgm(psill=psill,model=model,range=range,nugget=nugget)}else{
    
    vgm.mod		<-	vgm(psill=psill,model=model,range=range)				
    
  }
  
  g.dummy	<-	gstat(formula=formula,locations=locations,dummy=TRUE,beta=beta,model=vgm.mod,nmax=nmax)
  
  pred.vals	<-	predict(g.dummy,newdata=xy,nsim=nsim)
  
  if(nsim==1){
    mat	<-	matrix(pred.vals[,3],ncol=length(x.seq),nrow=length(y.seq),byrow=TRUE)
    if(is.null(ext)){ext	<-	extent(c(x.range,y.range))}		
    ras	<-		raster(mat,crs=proj4string,xmn=x.range[1],ymn=y.range[1],xmx=x.range[2],ymx=y.range[2])
    
  }else{
    pos			<-	3:(nsim+2)
    ras.list	<-	list()
    for(i in 1:nsim){
      
      ith.mat	<-	matrix(pred.vals[,pos[i]],ncol=length(x.seq),nrow=length(y.seq),byrow=TRUE)
      ras.list[[i]]<-	raster(ith.mat,crs=proj4string
                             ,xmn=x.range[1],ymn=y.range[1]
                             ,xmx=x.range[2],ymx=y.range[2])
      
    }
    
    ras	<-	stack(ras.list)
  }
  return(ras)
  
}


########################################################################################
# Function to simulate the distribution of species from a probability matrix
# p = probability matrix giving the presence probability of the species in a location
# p can be a raster file or a matrix representing the extent of the species output.
# n = the number of realized distributions of the species. 

spdist.rnd	<-	function(p,n){
  
  if(class(p)=="RasterLayer"){p.orig	<-	p; p	<-	as.matrix(p);}
  
  ras.dim	<-	dim(p)
  
  if(n==1){
    sp.dist	<-	matrix(0,nrow=ras.dim[1],ncol=ras.dim[2])
    for(i in 1:ras.dim[1]){
      for(j in 1:ras.dim[2]){
        
        sp.dist[i,j]	<-	rbinom(1,1,p[i,j])
        
      }	
    }
  }else{
    
    sp.dist	<-list()
    
    for(i in 1:n){
      
      ith.sp.dist	<- matrix(0,nrow=ras.dim[1],ncol=ras.dim[2])
      for(j in 1:ras.dim[1]){
        for(l in 1:ras.dim[2]){
          
          sp.dist[j,l]	<-	rbinom(1,1,p[j,l])
          
          
        }
      }
      
      sp.dist[[i]]	<-	ith.sp.dist
      
    }
    
  }
  if(is.null(p.orig)){out	<-	raster(sp.dist)}else{
    
    if(n==1){
      out	<-	raster(sp.dist,xmn=xmin(p.orig),xmx=xmax(p.orig),ymn=ymin(p.orig),ymx=ymax(p.orig)
                    ,crs=CRS(proj4string(p.orig)))
      
    }else{
      out	<-	list()
      for(i in 1:n){
        out[[i]]	<-	raster(sp.dist[[i]],xmn=xmin(p.orig),xmx=xmax(p.orig),ymn=ymin(p.orig),ymx=ymax(p.orig)
                           ,crs=CRS(proj4string(p.orig)))
        
      }
    }
  }
  
  return(out)
  
}



# Create random spatial variables. Two different sets, one wiht an Exponential model and
# the second one with a Spherical model.

r.extent	<-	extent(c(1,11.5,1,11.5))
var.set1		<-	rspatvar(formula=z~1,locations=~x+y,ext=r.extent,res=1/100,psill=1,range=10
                      ,model="Exp",nsim=5,beta=0,nmax=10)
var.set2		<-	rspatvar(formula=z~1,locations=~x+y,ext=r.extent,res=1/100,psill=1,range=6
                      ,model="Sph",nsim=5,beta=0,nmax=10)
var.set		<-	stack(var.set1,var.set2)

#save the environmental layers as tiff
for (i in 1:10){
  writeRaster(var.set,paste0('rbio',i),format='GTiff')
}

# we set the number of variables as 3
nvars			<-	3

##########################################################################################
# Simulation procedure:
# 1. Generate a coefficient for each independent variables using rnorm(),change the mean for weak-coefficient phenomenon 
# 2. Use Exp((X*betas)^2) transform the linear combination of variables and coefficients to give
#    a raster of presence probability.
# 3. Simulate the species distribution using the spdist.rnd() function
#    source. The function generates a distribution based on a bernoulli process for each
#    quadrant in the environmental space and the probability of presence of that quadrant.
# 4. Extract all of the quadrants in which the ith species is present.
# 5. Sample at random 50 quadrants in which the species is present.

# REPEAT the steps for strong-coefficient phenomenon
##########################################################################################
file.nams		<-	list.files("...",full.names=TRUE)
nspp			<-	100	# Number of Species to simulate
n.sample		<-	50
var.set			<-	stack(file.nams)
sp.dists		<-	stack()
sp.probs		<-	stack()
betas			<-	matrix(0,nrow=nspp,ncol=4,dimnames=list(paste("Sp",1:nspp,sep="_"),c("Intercept","beta1","beta2","beta3")))
var.id			<-	matrix(NA,nrow=nspp,ncol=3,dimnames=list(paste("Sp",1:nspp,sep="_"),c("Var1","Var2","Var3")))
spp.samp		<-	data.frame(SPP=rep(paste("Spp",1:nspp,sep="_"),each=n.sample),x=rep(0,nspp*n.sample),y=rep(0,nspp*n.sample))
start.ind		<-	seq(1,(nspp*n.sample),by=n.sample)
end.ind			<-	seq(n.sample,nspp*n.sample,by=n.sample)

for(i in 1:nspp){
  r.sample	<-	sample(1:10,3)
  ith.betas	<-	rnorm(3,mean=1,sd=0.5)
  intercept	<-	0
  tmp.stack	<-	var.set[[r.sample]]
  tmp.list	<-	list()
  for(j in 1:length(r.sample)){
    tmp.list[[j]]	<-	tmp.stack[[j]]*ith.betas[j]
  }
  
  betas[i,]	<-	c(intercept,ith.betas)
  var.id[i,]	<-	names(var.set)[r.sample]
  new.stack	<-	stack(tmp.list)
  lin.comb	<-	calc(new.stack,fun=sum)
  lin.comb	<-	intercept + lin.comb
  ith.p		<-	calc(lin.comb,function(x)exp(-(x)^2))
  ith.spp		<-	spdist.rnd(ith.p,1)
  sp.dists	<-	stack(sp.dists,ith.spp)
  sp.probs	<-	stack(sp.probs,ith.p)
  pts.val		<-	extract(ith.spp,coordinates(ith.spp))
  pres.coords	<-	coordinates(ith.spp)[pts.val==1,]
  pres.sampl	<-	sample(1:nrow(pres.coords),n.sample)
  spp.samp[start.ind[i]:end.ind[i],2:3]	<-	pres.coords[pres.sampl,]
  
}

names(sp.dists)	<-	names(sp.probs)	<-	paste("Sp",1:nspp,sep="_")
##########################################################################################
# Save the simulation in Rdata, export the coordinates of the random points and export
# the rasters of the simulated environmental data
save(list=c("sp.dists","sp.probs","betas","var.id")
     ,file="species_simulations.RData")
write.table(spp.samp
            ,"Presence_Points.txt"
            ,sep="\t",row.names=FALSE,quote=FALSE)


