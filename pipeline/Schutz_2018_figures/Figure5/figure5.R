#!/usr/bin/env Rscript
# Descriptive evidence-domain structure only. No outcomes or prediction results used.
invisible(Sys.setlocale('LC_ALL','Chinese (Simplified)_China.utf8'))
sa<-grep('^--file=',commandArgs(FALSE),value=TRUE)
out<-dirname(normalizePath(sub('^--file=','',sa[1]),winslash='/'))
base<-normalizePath(file.path(out,'../..'),winslash='/');src<-file.path(out,'figure5_plot_data','v2')
.libPaths(c(normalizePath(file.path(out,'../R_lib'),winslash='/'),.libPaths()))
pk<-c('ggplot2','patchwork','ragg','sf','maps','ComplexUpset','circlize','WeightedTreemaps','gridGraphics','ggrepel','digest','jsonlite')
missing<-pk[!vapply(pk,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))stop('Missing R dependencies: ',paste(missing,collapse=', '))
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(grid)})
dir.create(src,showWarnings=FALSE)
hash<-function(p)digest::digest(file=p,algo='sha256',serialize=FALSE)
read<-function(p)read.delim(p,check.names=FALSE,stringsAsFactors=FALSE,na.strings='NA')
wt<-function(d,f)write.table(d,file.path(src,f),sep='\t',quote=TRUE,row.names=FALSE,na='NA',fileEncoding='UTF-8')
files<-c(matrix='Schutz_2018_phase0/analysis_matrix.tsv',dictionary='Schutz_2018_phase0/data_dictionary.tsv',
 metadata='Schutz_2018_phase_minus1/phase_minus1_cohort_manifest.tsv',study_metadata='Schutz_2018_phase_minus1/phase_minus1_study_audit.tsv',
 crop_order='Schutz_2018_figures/Figure1/figure1_plot_data/crop_annotations_and_order.tsv',
 country_counts='Schutz_2018_figures/Figure1/figure1_plot_data/country_support.tsv',
 F1='Schutz_2018_figures/Figure1/Figure1_v2.pdf',F2='Schutz_2018_figures/Figure2/Figure2_final.pdf',
 F3='Schutz_2018_figures/Figure3/Figure3_final.pdf',F4='Schutz_2018_figures/Figure4/Figure4_final.pdf',
 F4_png='Schutz_2018_figures/Figure4/Figure4_final.png',F4_script='Schutz_2018_figures/Figure4/figure4.R',
 F4_state='Schutz_2018_figures/Figure4/FIGURE4_FINAL_STATE.md')
before<-setNames(vapply(file.path(base,files),hash,character(1)),names(files))
stopifnot(before[['matrix']]=='66fd6e7ddd65426225167bfa3016f518cad26ff836e0748decdd3afa39209094')
rd<-function(k)read(file.path(base,files[[k]]))
dictionary<-rd('dictionary');predictors<-dictionary$column[dictionary$role=='main_predictor']
stopifnot(length(predictors)==12)
wt(dictionary[dictionary$role=='main_predictor',],'frozen_predictor_dictionary.tsv')
dat<-rd('matrix')[,c('observation_id','study_id','pair_id',predictors)]
stopifnot(nrow(dat)==1723,length(unique(dat$study_id))==169,length(unique(dat$crop))==40,!anyDuplicated(dat$observation_id))
meta<-rd('metadata');key<-paste(dat$study_id,dat$pair_id,sep=':');mkey<-paste(meta$study_id,meta$pair_id,sep=':')
stopifnot(!anyDuplicated(mkey),setequal(key,mkey))
dat$country<-meta$country[match(key,mkey)];stopifnot(!anyNA(dat$country),length(unique(dat$country))==31)
wt(dat[,c('observation_id','study_id','pair_id','country','crop')],'observation_country_crop_membership.tsv')
FONT<-'Arial';INK<-'#242A30';MID<-'#66737D';RULE<-'#CBD2D7';BLUE<-'#59798F';GREEN<-'#4F7772'
theme_set(theme_classic(base_size=7,base_family=FONT)+theme(axis.line=element_line(linewidth=.25,colour=INK),
 axis.ticks=element_line(linewidth=.25,colour=INK),axis.text=element_text(size=5.7,colour=INK),axis.title=element_text(size=6.2),
 plot.margin=margin(2,3,2,2),legend.text=element_text(size=5.6),legend.title=element_text(size=6.2),panel.grid=element_blank()))

# A: geographic counts at representative country points, not trial coordinates.
country<-aggregate(rep(1,nrow(dat)),list(country=dat$country),sum);names(country)[2]<-'comparisons'
country$studies<-vapply(country$country,function(c)length(unique(dat$study_id[dat$country==c])),integer(1))
country<-country[order(-country$comparisons,country$country),];country$rank<-1:nrow(country)
prev<-rd('country_counts');stopifnot(all(country$comparisons==prev$comparisons[match(country$country,prev$country)]),
 all(country$studies==prev$studies[match(country$country,prev$country)]))
world<-sf::st_as_sf(maps::map('world',plot=FALSE,fill=TRUE))
crs_robin<-'+proj=robin +lon_0=0 +datum=WGS84 +units=m +no_defs'
world<-sf::st_transform(world,crs_robin);world<-world[world$ID!='Antarctica',]
world<-sf::st_simplify(world,dTolerance=5000)
idx<-match(country$country,world$ID);stopifnot(!anyNA(idx))
countrypts<-lapply(idx,function(i){g<-suppressWarnings(sf::st_cast(sf::st_geometry(world[i,]),'POLYGON'));g<-g[which.max(sf::st_area(g))];suppressWarnings(sf::st_point_on_surface(g))[[1]]})
points<-sf::st_sfc(countrypts,crs=sf::st_crs(world));xy<-sf::st_coordinates(points)
country$x<-xy[,1];country$y<-xy[,2];country$label<-ifelse(country$rank<=8,paste0(country$country,'\n',country$studies,' studies'),'')
nx<-c(India=2000000,Iran=1000000,Turkey=-1500000,Pakistan=2100000,Egypt=-1500000,Brazil=1300000,Canada=-1700000,Argentina=-1300000)
ny<-c(India=-1600000,Iran=2000000,Turkey=1400000,Pakistan=1600000,Egypt=-600000,Brazil=700000,Canada=-700000,Argentina=-600000)
wt(country,'panel_A_country_evidence.tsv')
sf::st_write(world,file.path(src,'panel_A_world_basemap.geojson'),delete_dsn=TRUE,quiet=TRUE)
saveRDS(world,file.path(src,'panel_A_world_basemap.rds'))
pA<-ggplot()+geom_sf(data=world,fill='#EAECED',colour='white',linewidth=.10)+
 geom_point(data=country,aes(x,y,size=comparisons),shape=21,fill=BLUE,colour='#4A6577',stroke=.20,alpha=.72)+
 ggrepel::geom_text_repel(data=country[country$rank<=8,],aes(x,y,label=label),family=FONT,size=2.2,lineheight=1.02,
  nudge_x=nx[country$country[country$rank<=8]],nudge_y=ny[country$country[country$rank<=8]],
  seed=5001,min.segment.length=0,segment.colour='#7E8B94',segment.size=.20,box.padding=.6,point.padding=.6,max.overlaps=Inf,force=4)+
 scale_size_area(max_size=10,breaks=c(50,200,600),name='Comparisons')+
 coord_sf(crs=crs_robin,datum=NA,expand=FALSE)+theme_void(base_family=FONT)+
 theme(legend.position='bottom',legend.direction='horizontal',legend.box='horizontal',legend.key.height=unit(2.6,'mm'),
  legend.key.width=unit(4,'mm'),legend.title=element_text(size=6),legend.text=element_text(size=5.8),plot.margin=margin(0,0,0,0))

# B: dictionary NA tokens only; zero and valid nominal 'unclear' are observed.
miss<-is.na(dat[predictors]);colnames(miss)<-predictors
avail<-rowSums(!miss)
missingness<-data.frame(predictor=predictors,missing_n=colSums(miss),observed_n=nrow(dat)-colSums(miss))
missingness$missing_fraction<-missingness$missing_n/nrow(dat)
stopifnot(identical(as.integer(missingness$missing_n),as.integer(dictionary$missing_n[dictionary$role=='main_predictor'])))
labels<-c(crop='Crop',climate_zone='Climate',soil_texture='Soil texture',soil_organic_matter_pct='Soil OM',soil_ph_water='Soil pH',
 available_p_kg_ha='Available P',inoculant_function='Inoculant function',n_applied_kg_ha='N applied',k_applied_kg_ha='K applied',
 p_applied_kg_ha='P applied',irrigation_status='Irrigation',fertilization_type='Fertilization type')
missingness$display_label<-labels[predictors]
wt(missingness,'panel_B_predictor_missingness.tsv')
wt(cbind(dat[,c('observation_id','study_id')],as.data.frame(miss),available_predictors=avail),'panel_B_observation_missingness.tsv')
# Stable dictionary-order binary patterns: 1 missing, 0 observed.
binary<-apply(miss,1,function(z)paste(as.integer(z),collapse=''))
patterns<-as.data.frame(table(binary));names(patterns)<-c('pattern','comparisons');patterns$pattern<-as.character(patterns$pattern)
patterns<-patterns[order(-patterns$comparisons,patterns$pattern),]
patterns$missing_predictors<-vapply(patterns$pattern,function(p)paste(predictors[strsplit(p,'',fixed=TRUE)[[1]]=='1'],collapse=';'),character(1))
patterns$degree<-vapply(patterns$pattern,function(p)sum(strsplit(p,'',fixed=TRUE)[[1]]=='1'),integer(1))
wt(patterns,'panel_B_all_exclusive_missingness_patterns.tsv')
partial<-missingness[missingness$missing_n>0,];partial<-partial[order(-partial$missing_n,partial$predictor),]
ud<-as.data.frame(miss[,partial$predictor]);names(ud)<-partial$display_label
top<-head(patterns[patterns$degree>0,],12)
intlist<-lapply(top$missing_predictors,function(s)unname(labels[strsplit(s,';',fixed=TRUE)[[1]]]))
wt(top,'panel_B_displayed_exclusive_patterns.tsv')
# Lollipop frequency annotations retain the UpSet combination matrix; no plain bar panel.
frequency<-ggplot(mapping=aes(x=intersection,y=in_exclusive_intersection))+
 stat_summary(fun=sum,fun.min=function(z)0,fun.max=sum,geom='linerange',colour=BLUE,linewidth=.35,na.rm=TRUE)+
 stat_summary(fun=sum,geom='point',colour=BLUE,size=1.4,na.rm=TRUE)+
 stat_summary(fun=sum,geom='text',aes(label=after_stat(y)),vjust=-.5,size=1.9,family=FONT,colour=INK,na.rm=TRUE)
frequency<-frequency+labs(y='Comparisons')+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
 axis.text.y=element_text(size=5.7),axis.title.y=element_text(size=6.2),plot.margin=margin(5,2,0,2))+
 scale_y_continuous(expand=expansion(mult=c(0,.16)))
single<-ComplexUpset::upset_set_size(geom=geom_point(stat='count',size=1.6,colour=BLUE))+
 labs(y='Missing comparisons')+theme(axis.text=element_text(size=5.7),axis.title=element_text(size=6.2),plot.margin=margin(0,2,0,2))
pB<-ComplexUpset::upset(ud,intersect=names(ud),n_intersections=12,sort_sets=FALSE,
 name='Missing predictors',min_degree=1,mode='distinct',encode_sets=FALSE,keep_empty_groups=TRUE,width_ratio=.20,height_ratio=.65,
 base_annotations=list('Comparisons'=frequency),set_sizes=single,
 matrix=ComplexUpset::intersection_matrix(geom=geom_point(size=1.5),segment=geom_segment(linewidth=.35),
  outline_color=list(active=INK,inactive='#DEE3E5')),
 themes=ComplexUpset::upset_modify_themes(list(default=theme(text=element_text(family=FONT,size=6),axis.text=element_text(size=6),
  panel.grid=element_blank(),plot.background=element_rect(fill='white',colour=NA),legend.position='none'),
  intersections_matrix=theme(plot.margin=margin(0,2,0,2),axis.text.y=element_text(size=6,family=FONT),axis.title.x=element_blank(),axis.text.x=element_blank()),
  overall_sizes=theme(axis.text=element_text(size=5.7,family=FONT),axis.title=element_text(size=6.2,family=FONT)))))
# The empty upper-left spacer otherwise imposes a margin at the row boundary.
pB$patches$plots[[1]]<-pB$patches$plots[[1]]+theme(plot.margin=margin(0,0,0,0))
upset_structure<-ComplexUpset::upset_data(ud,names(ud),n_intersections=12,min_degree=1,
 keep_empty_groups=TRUE,sort_sets=FALSE,encode_sets=FALSE)
wt(upset_structure$matrix_frame,'panel_B_displayed_matrix_cells.tsv')
display_names<-rev(upset_structure$sorted$intersections)
displayed<-data.frame(display_rank=seq_along(display_names),intersection=display_names,
 comparisons=unname(upset_structure$sizes$exclusive_intersection[display_names]))
wt(displayed,'panel_B_displayed_intersection_counts.tsv')
stopifnot(setequal(displayed$comparisons,top$comparisons),length(upset_structure$plot_sets_subset)==9)
# Audit the built annotation coordinates as well as their underlying tables.
frequency_built<-ggplot_build(pB$patches$plots[[2]])
freq_points<-frequency_built$data[[2]][,c('x','y')]
freq_points<-freq_points[order(freq_points$x),]
stopifnot(nrow(freq_points)==12,identical(as.integer(freq_points$y),as.integer(displayed$comparisons)))
wt(data.frame(display_rank=as.integer(freq_points$x),rendered_comparisons=freq_points$y),'panel_B_rendered_intersection_count_QC.tsv')
side_built<-ggplot_build(pB$patches$plots[[3]])
side_points<-side_built$data[[2]][,c('x','y')]
stopifnot(nrow(side_points)==9,identical(sort(as.integer(abs(side_points$y))),sort(as.integer(partial$missing_n))))
wt(data.frame(group_rank=as.integer(side_points$x),rendered_missing_n=abs(side_points$y)),
 'panel_B_rendered_marginal_count_QC.tsv')

# C: study-average completeness; retain alternative summaries in source audit only.
study<-aggregate(list(comparisons=rep(1,nrow(dat)),available_predictors=avail),list(study_id=dat$study_id),sum)
study$available_predictors<-study$available_predictors/study$comparisons;study$completeness<-study$available_predictors/12
study$predictors_observed_in_all_rows<-vapply(study$study_id,function(i)sum(colSums(miss[dat$study_id==i,,drop=FALSE])==0),integer(1))
study$predictors_observed_in_any_row<-vapply(study$study_id,function(i)sum(colSums(!miss[dat$study_id==i,,drop=FALSE])>0),integer(1))
study<-study[order(-study$comparisons,study$study_id),];study$rank<-1:169
cache<-file.path(src,'panel_C_refined_voronoi.rds')
if(!file.exists(cache)){
 source(file.path(out,'refine_weighted_voronoi.R'))
 input<-data.frame(study=as.character(study$study_id),size=study$comparisons)
 tm<-WeightedTreemaps::voronoiTreemap(input,levels='study',cell_size='size',filter=0,seed=5003,positioning='random',
  shape='rectangle',maxIteration=500,error_tol=.0001,verbose=FALSE)
 result<-refine_voronoi(tm);saveRDS(result,cache)
}
vo<-readRDS(cache);target<-study$comparisons[match(vo$study_id,study$study_id)]/nrow(dat)*4000000
stopifnot(length(vo$polygons)==169,setequal(vo$study_id,study$study_id),max(abs(vo$areas/target-1))<1e-5)
geometry<-sf::st_sf(study_id=vo$study_id,geometry=vo$polygons)
geometry<-merge(geometry,study,by='study_id',sort=FALSE)
geometry$area_fraction<-as.numeric(sf::st_area(geometry))/sum(as.numeric(sf::st_area(geometry)))
geometry$area_relative_error<-geometry$area_fraction/(geometry$comparisons/nrow(dat))-1
stopifnot(all(sf::st_is_valid(geometry)),max(abs(geometry$area_relative_error))<1e-5)
vertices<-sf::st_coordinates(geometry);vertices<-data.frame(study_id=geometry$study_id[vertices[,'L2']],x=vertices[,'X'],y=vertices[,'Y'])
vertices$completeness<-geometry$completeness[match(vertices$study_id,geometry$study_id)]
wt(sf::st_drop_geometry(geometry),'panel_C_study_size_completeness.tsv');wt(vertices,'panel_C_polygon_vertices.tsv');wt(vo$history,'panel_C_graphical_area_calibration.tsv')
centroids<-sf::st_coordinates(sf::st_centroid(sf::st_geometry(geometry)))
cent<-cbind(sf::st_drop_geometry(geometry),x=centroids[,1],y=centroids[,2]);cent$label<-sprintf('S%03d',cent$study_id)
cent<-cent[cent$rank<=4,];wt(cent,'panel_C_largest_study_labels.tsv')
# Display bins only: unchanged comparison-mean predictor availability, no rounding.
bin_labels<-c('12/12','9–11/12','6–8/12','≤5/12')
classify<-function(v)factor(ifelse(v==12,'12/12',ifelse(v>=9,'9–11/12',ifelse(v>=6,'6–8/12','≤5/12'))),levels=bin_labels)
study$display_completeness_bin<-classify(study$available_predictors)
vertices$display_completeness_bin<-classify(geometry$available_predictors[match(vertices$study_id,geometry$study_id)])
wt(study[,c('study_id','available_predictors','display_completeness_bin')],'panel_C_display_bins.tsv')
bin_summary<-data.frame(display_bin=bin_labels,studies=as.integer(table(study$display_completeness_bin)),comparisons=vapply(bin_labels,function(b)sum(study$comparisons[study$display_completeness_bin==b]),numeric(1)))
wt(bin_summary,'panel_C_display_bin_summary.tsv')
stopifnot(identical(bin_summary$studies,c(66L,85L,16L,2L)),sum(bin_summary$comparisons)==1723)
pC<-ggplot(vertices,aes(x,y,group=study_id,fill=display_completeness_bin))+geom_polygon(colour='white',linewidth=.25)+
 geom_text(data=cent,aes(x,y,label=label),inherit.aes=FALSE,family=FONT,size=2.1,colour=INK)+
 scale_fill_manual(values=setNames(c(GREEN,'#8DA6A0','#C0CFCA','#E5EAE8'),bin_labels),drop=FALSE,
  name='Available predictors / 12',guide=guide_legend(title.position='top',nrow=1,
    keywidth=unit(3,'mm'),keyheight=unit(2.5,'mm')))+
 coord_equal(expand=FALSE)+theme_void(base_family=FONT)+theme(legend.position='bottom',legend.title=element_text(size=6.2),
  legend.text=element_text(size=5.7),legend.margin=margin(2,0,5,0),plot.margin=margin(0,0,0,0))

# D: top domains by structural support only; all remaining rows retained in Other.
crop_order<-rd('crop_order');crop_order<-crop_order[order(crop_order$rank),]
top_country<-head(country$country,8);top_crop<-head(crop_order$crop,8)
dc<-data.frame(observation_id=dat$observation_id,study_id=dat$study_id,
 country=ifelse(dat$country%in%top_country,dat$country,'Other countries'),
 crop=ifelse(dat$crop%in%top_crop,crop_order$display_label[match(dat$crop,crop_order$crop)],'Other crops'))
wt(dc,'panel_D_observation_grouping.tsv')
cc<-aggregate(rep(1,nrow(dc)),list(country=dc$country,crop=dc$crop),sum);names(cc)[3]<-'comparisons'
wt(cc,'panel_D_country_crop_comparisons.tsv')
country_nodes<-c(top_country,'Other countries');crop_nodes<-c(crop_order$display_label[match(top_crop,crop_order$crop)],'Other crops')
country_col<-setNames(colorRampPalette(c('#54758C','#BCCBD3'))(9),country_nodes)
crop_col<-setNames(colorRampPalette(c('#9A846F','#D6CCBC'))(9),crop_nodes)
sector_col<-c(country_col,crop_col);sectors<-c(country_nodes,crop_nodes)
wt(data.frame(sector=sectors,domain=rep(c('Country','Crop'),each=9),colour=unname(sector_col[sectors])),'panel_D_sector_palette_order.tsv')
stopifnot(sum(cc$comparisons)==1723,length(unique(dc$observation_id))==1723)
chord_fun<-function(){
 par(mar=c(0,0,0,0),family=FONT,cex=.75)
 circlize::circos.clear()
 circlize::circos.par(start.degree=270,gap.after=c(rep(2,8),12,rep(2,8),12),
  track.margin=c(.001,.001),cell.padding=c(0,0,0,0),canvas.xlim=c(-1.28,1.28),canvas.ylim=c(-1.22,1.22))
 links<-circlize::chordDiagram(cc,order=sectors,grid.col=sector_col,col=country_col[cc$country],transparency=.57,
  annotationTrack='grid',annotationTrackHeight=.038,preAllocateTracks=list(track.height=.12),
  directional=0,link.border=NA,link.sort=TRUE,link.decreasing=TRUE,reduce=0)
 circlize::circos.trackPlotRegion(track.index=1,panel.fun=function(x,y){
  sector<-circlize::get.cell.meta.data('sector.index');lim<-circlize::get.cell.meta.data('xlim')
  circlize::circos.text(mean(lim),.35,sector,facing='clockwise',niceFacing=TRUE,adj=c(0,.5),cex=.60,col=INK)
 },bg.border=NA)
 circlize::circos.clear()
 invisible(links)
}

txt<-function(label,x,y,size=7,col=INK,face='plain',just='left')grid.text(label,x=unit(x,'npc'),y=unit(y,'npc'),just=just,gp=gpar(fontfamily=FONT,fontsize=size,col=col,fontface=face))
header<-function(tag,title){yy<-1-convertHeight(unit(2.4,'mm'),'npc',valueOnly=TRUE);xx<-.02+convertWidth(unit(4.4,'mm'),'npc',valueOnly=TRUE);txt(tag,.02,yy,9.5,face='bold');txt(title,xx,yy,8.5,face='bold')}
embed<-function(p,x=.5,y=.45,w=.98,h=.76){pushViewport(viewport(x=x,y=y,width=w,height=h));grid.draw(if(inherits(p,'patchwork'))patchwork::patchworkGrob(p)else ggplotGrob(p));popViewport()}
draw_A<-function(){header('A','Global evidence footprint');txt(sprintf('%d studies / %d countries / %d comparisons',length(unique(dat$study_id)),nrow(country),nrow(dat)),.98,.969,6.0,col=MID,just='right');embed(pA,y=.443,h=.866)}
draw_B<-function(){header('B','Predictor missingness');embed(pB,y=.477,w=.978,h=.826)
 txt(paste0('No missing values: ',paste(labels[missingness$predictor[missingness$missing_n==0]],collapse=', ')),.025,.050,5.8,col=MID)
 txt(sprintf('Complete comparisons: %d/%d · Top %d nonempty patterns shown',sum(avail==12),nrow(dat),nrow(top)),.025,.020,5.8,col=MID)}
draw_C<-function(){header('C','Study metadata landscape');embed(pC,y=.439,w=.95,h=.816)
 txt('Polygon area: comparisons per study',.025,.027,5.8,col=MID)}
draw_D<-function(){
 echo_file<-tempfile(fileext='.pdf');on.exit(unlink(echo_file))
 pushViewport(viewport(x=.5,y=.452,width=.97,height=.816))
 gridGraphics::grid.echo(chord_fun,newpage=FALSE,prefix='chord',
  device=function(width,height){grDevices::cairo_pdf(echo_file,width=width,height=height,family=FONT);grDevices::dev.control(displaylist='enable')})
 popViewport()
 header('D','Country–crop evidence coupling')
 txt('Countries',.08,.89,6.2,col=country_col[1]);txt('Crops',.82,.89,6.2,col=crop_col[1])
 txt('Ribbon width: comparisons',.025,.027,5.8,col=MID)}
capture_panel<-function(fun,w,h){tmp<-tempfile(fileext='.pdf');on.exit(unlink(tmp));grid.grabExpr(fun(),width=w,height=h,wrap=TRUE,wrap.grobs=TRUE,device=function(width,height)grDevices::cairo_pdf(tmp,width=width,height=height,family=FONT))}
width_mm<-183;height_mm<-247;hA<-74;hB<-85;hBottom<-88
gA<-capture_panel(draw_A,width_mm/25.4,hA/25.4);gB<-capture_panel(draw_B,width_mm/25.4,hB/25.4)
gC<-capture_panel(draw_C,90/25.4,hBottom/25.4);gD<-capture_panel(draw_D,93/25.4,hBottom/25.4)
wrap<-function(g)wrap_elements(full=g,clip=FALSE)
bottom<-(wrap(gC)|wrap_elements(full=gD,clip=TRUE))+plot_layout(widths=c(90,93))
fig<-(wrap(gA)/wrap(gB)/bottom)+plot_layout(heights=c(hA,hB,hBottom))
fig<-fig & theme(plot.background=element_rect(fill='white',colour=NA),plot.margin=margin(0,0,0,0))
pdf_file<-file.path(out,'Figure5_v2.pdf');png_file<-file.path(out,'Figure5_v2.png')
grDevices::cairo_pdf(pdf_file,width=width_mm/25.4,height=height_mm/25.4,family=FONT);print(fig);invisible(dev.off())
ragg::agg_png(png_file,width=width_mm,height=height_mm,units='mm',res=600,background='white');print(fig);invisible(dev.off())
for(it in list(list(name='preview',plot=fig,height=height_mm,res=180),list(name='A',plot=wrap(gA),height=hA,res=300),list(name='B',plot=wrap(gB),height=hB,res=300),list(name='C_D',plot=bottom,height=hBottom,res=300))){
 ragg::agg_png(file.path(src,paste0('Figure5_QA_',it$name,'.png')),width=width_mm,height=it$height,units='mm',res=it$res,background='white');print(it$plot);invisible(dev.off())}
after<-setNames(vapply(file.path(base,files),hash,character(1)),names(files));stopifnot(identical(before,after))
con<-file(png_file,'rb');invisible(readBin(con,'raw',8));dims<-ppm<-NULL
repeat{len<-readBin(con,'integer',1,size=4,endian='big');if(!length(len))break;typ<-rawToChar(readBin(con,'raw',4));raw<-readBin(con,'raw',len);readBin(con,'raw',4);if(typ=='IHDR')dims<-readBin(raw[1:8],'integer',2,size=4,endian='big');if(typ=='pHYs')ppm<-readBin(raw[1:8],'integer',2,size=4,endian='big');if(typ=='IEND')break};close(con)
stopifnot(all(abs(dims-c(width_mm,height_mm)/25.4*600)<2),all(abs(ppm*.0254-600)<.02))
caption<-c('Figure 5 | Metadata completeness and evidence-domain boundaries.',
 'All panels describe the frozen Schütz main cohort (1,723 comparisons, 169 studies). No outcome, prediction, model error, feature-importance analysis or new statistical test is used. These structures are not interpreted as causal explanations of prediction failure.',
 'A, Country evidence footprint. Bubble area is proportional to comparison count; labels for the eight countries with most comparisons show country names and numbers of studies. Counts are regenerated from the frozen comparison matrix connected to existing country metadata by unique study/pair IDs. Bubbles are located at representative points of each country’s largest mapped land polygon and are not trial locations. The light-grey Robinson-projection basemap uses the locally installed maps world database, omitting Antarctica from display. Study counts across countries are not additive because one study can contribute evidence from several countries.',
 sprintf('B, Missingness of exactly the 12 frozen main predictors, using only the matrix/dictionary NA definition; zero values and valid nominal categories remain observed. Dots at left show total missing comparisons for each of the %d predictors with missing data. The UpSet matrix shows the %d most frequent exclusive, nonempty missingness combinations, and lollipop annotations show their comparison counts. A filled dot marks a missing predictor; a pale dot marks an observed predictor. Predictors without missing values are listed separately. %d comparisons have all 12 predictors observed. Marginal counts use the full cohort, not only displayed combinations; omitted rare patterns and complete cases are retained in source data. No missingness mechanism is inferred.',nrow(partial),nrow(top),sum(avail==12)),
 sprintf('C, Area-weighted study tessellation. Each of all 169 polygons represents one study; target area is comparison count divided by 1,723. An additively weighted Voronoi layout is numerically calibrated to these targets (maximum relative area error %.3g). Layout positions/neighbours do not encode geographic location, study similarity or time. Color is comparison-mean predictor availability within each study: mean number of its 12 observed predictors per comparison, divided by 12. A predictor contributes one unit when nonmissing on a comparison; no imputation is performed. Four display categories preserve this unrounded study mean: 12/12 denotes exactly 12; 9–11/12 denotes [9,12); 6–8/12 denotes [6,9); ≤5/12 denotes <6 available predictors. The interval definitions account for fractional study means; categories contain 66, 85, 16 and 2 studies, respectively. No completeness values or polygon areas are changed. The four largest studies are labelled with frozen study IDs. Alternative any-row/all-row availability summaries are provided in source data for definition transparency, not as additional plotted analyses.',max(abs(geometry$area_relative_error))),
 'D, Country–crop bipartite chord. The eight countries with most comparisons are retained individually; 23 others form Other countries. The eight crops with greatest study support are retained individually, using the unchanged frozen crop order/tie rules; 32 others form Other crops. Ribbon width represents the number of comparisons connecting these groups, summing to 1,723. Countries and crops use separate restrained blue and earth-tone palettes. Grouping is based on structural support only, and all observations are retained. This panel displays evidence-domain coupling and is not intended for precise pairwise ranking. Ribbons do not imply independent experiments, direction of causation or any effect-size association.',
 'Source data include comparison-level country/crop membership and missingness, complete marginal/intersection counts, country representatives/basemap, every study’s completeness/area/vertices, and all aggregated and observation-level chord memberships. All graphics, previews and exports use R. Frozen source files and Figures 1–4 remain unchanged.')
writeLines(caption,file.path(src,'Figure5_v2_caption.txt'),useBytes=TRUE)
writeLines(c('# Figure 5 plot data','',
 'A: country comparison/study counts and representative map points; local maps world database in Robinson projection. No field-coordinate claim.',
 'B: missingness computed from all 12 dictionary main predictors. All exclusive patterns and observation memberships retained; top nonempty patterns alone displayed.',
 'C: 169 area-weighted polygons. Display color = four bins of unchanged mean observed predictors per comparison within study / 12. Numerical area calibration changes drawing geometry only; no models or statistical tests.',
 'D: top-eight countries by comparisons and top-eight crops by study support, plus complete Other groups. All 1723 memberships and aggregate counts saved.',
 'Reproduce with figure5.R. The refined Voronoi RDS is a geometry cache; remove it only to regenerate the layout using refine_weighted_voronoi.R. It contains no outcome or prediction data.',
 'QA previews are R-generated; caption, session info and QC accompany source tables.'),file.path(src,'README.md'),useBytes=TRUE)
writeLines(c('https://krassowski.github.io/complex-upset/articles/Examples_R.html','https://stat.ethz.ch/CRAN/web/packages/WeightedTreemaps/refman/WeightedTreemaps.html','https://jokergoo.github.io/circlize_book/book/the-chorddiagram-function.html'),file.path(src,'plotting_documentation.txt'))
capture.output(sessionInfo(),file=file.path(src,'R_sessionInfo.txt'))
qc<-list(status='RENDERED_PENDING_QC',source_sha256=as.list(before),sources_unchanged=TRUE,
 comparisons=nrow(dat),studies=length(unique(dat$study_id)),countries=nrow(country),predictors=length(predictors),
 no_missing_predictors=missingness$predictor[missingness$missing_n==0],complete_comparisons=sum(avail==12),
 exclusive_missingness_patterns=nrow(patterns),displayed_patterns=nrow(top),study_polygons=169,
 maximum_relative_polygon_area_error=max(abs(geometry$area_relative_error)),
 chord_comparisons=sum(cc$comparisons),outcome_values_used=0,predictions_read=0,models_fitted=0,new_statistical_tests=0,
 font=FONT,width_mm=width_mm,height_mm=height_mm,PNG_pixels=as.list(dims),PNG_dpi=as.list(ppm*.0254),
 output_sha256=list(pdf=hash(pdf_file),png=hash(png_file)),script_sha256=hash(file.path(out,'figure5.R')))
jsonlite::write_json(qc,file.path(src,'figure5_QC.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
message('Figure 5 rendered: ',paste(dims,collapse=' x '),' pixels, 600 dpi.')
capture.output(warnings(),file=file.path(src,'R_render_warnings.txt'))
