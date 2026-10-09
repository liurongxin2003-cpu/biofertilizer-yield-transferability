# Pure graphical area calibration; not a statistical or predictive model.
refine_voronoi<-function(tm,tolerance=1e-5){
 sites<-do.call(rbind,lapply(tm@cells,function(z)c(z$site,z$weight)))
 target<-vapply(tm@cells,function(z)z$target,numeric(1))*4000000
 box<-sf::st_sfc(sf::st_polygon(list(rbind(c(0,0),c(2000,0),c(2000,2000),c(0,2000),c(0,0)))))
 diagram<-function(w){
  s<-sites;s[,3]<-w
  v<-WeightedTreemaps::cropped_voronoi(s)
  g<-sf::st_sfc(lapply(v,function(z){b<-as.matrix(z$border);if(any(b[1,]!=b[nrow(b),]))b<-rbind(b,b[1,]);sf::st_polygon(list(b))}))
  g<-suppressWarnings(sf::st_intersection(g,box))
  list(polygons=g,areas=as.numeric(sf::st_area(g)))
 }
 w<-sites[,3];current<-diagram(w);n<-length(w)
 history<-data.frame(iteration=0,max_relative_error=max(abs(current$areas/target-1)))
 for(iter in 1:8){
  err<-max(abs(current$areas/target-1));message('Voronoi graphical area refinement ',iter,': max relative error ',signif(err,5))
  if(err<tolerance)break
  jac<-matrix(0,n-1,n-1);step<-.02
  for(j in 1:(n-1)){
   trial<-w;trial[j]<-trial[j]+step
   ar<-diagram(trial)$areas
   stopifnot(length(ar)==n)
   jac[,j]<-(ar[1:(n-1)]-current$areas[1:(n-1)])/step
  }
  dw<-c(solve(jac,target[1:(n-1)]-current$areas[1:(n-1)]),0)
  accepted<-FALSE
  for(alpha in c(1,.5,.25,.125,.0625)){
   trial<-diagram(w+alpha*dw)
   if(length(trial$areas)==n&&max(abs(trial$areas/target-1))<err){w<-w+alpha*dw;current<-trial;accepted<-TRUE;break}
  }
  if(!accepted)stop('Weighted Voronoi area refinement failed to improve. No substitute chart generated.')
  history<-rbind(history,data.frame(iteration=iter,max_relative_error=max(abs(current$areas/target-1))))
 }
 stopifnot(length(current$areas)==n,max(abs(current$areas/target-1))<tolerance)
 list(polygons=current$polygons,areas=current$areas,targets=target,sites=sites,weights=w,
  study_id=as.integer(vapply(tm@cells,function(z)z$name,character(1))),history=history)
}
