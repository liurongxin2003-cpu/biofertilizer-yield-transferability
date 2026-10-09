#!/usr/bin/env Rscript
# Frozen Task A visualization only. No training, resampling, density estimation or new inference.
invisible(Sys.setlocale('LC_ALL','Chinese (Simplified)_China.utf8'))
sa<-grep('^--file=',commandArgs(trailingOnly=FALSE),value=TRUE)
out<-if(length(sa))dirname(normalizePath(sub('^--file=','',sa[1]),winslash='/'))else normalizePath('Schutz_2018_figures/Figure3',winslash='/')
base<-normalizePath(file.path(out,'../..'),winslash='/')
src<-file.path(out,'figure3_plot_data','v4');oldsrc<-file.path(out,'figure3_plot_data','v3')
.libPaths(c(file.path(base,'Schutz_2018_figures','R_lib'),.libPaths()))
required<-c('ggplot2','ggalluvial','patchwork','ragg','digest','jsonlite')
missing<-required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))stop('Missing R packages: ',paste(missing,collapse=', '))
suppressPackageStartupMessages({library(ggplot2);library(ggalluvial);library(patchwork);library(grid)})
dir.create(src,recursive=TRUE,showWarnings=FALSE)
hash<-function(p)digest::digest(file=p,algo='sha256',serialize=FALSE)
read<-function(p)read.delim(p,check.names=FALSE,stringsAsFactors=FALSE)
write_tsv<-function(x,f)write.table(x,file.path(src,f),sep='\t',quote=FALSE,row.names=FALSE,na='NA',fileEncoding='UTF-8')
files<-c(predictions='Schutz_2018_phase1/phase1_predictions.tsv',metrics='Schutz_2018_phase1/phase1_metrics.tsv',
 bootstrap='Schutz_2018_phase1/phase1_bootstrap.tsv',archive='Schutz_2018_phase1/run/bootstrap_draws.npz',
 protocol='Schutz_2018_phase1/evidence/protocol.json',classification='Schutz_2018_phase1/evidence/classification_evidence.json',
 Figure1_pdf='Schutz_2018_figures/Figure1/Figure1_v2.pdf',Figure1_png='Schutz_2018_figures/Figure1/Figure1_v2.png',
 Figure2_pdf='Schutz_2018_figures/Figure2/Figure2_final.pdf',Figure2_png='Schutz_2018_figures/Figure2/Figure2_final.png',
 Figure3_v3_pdf='Schutz_2018_figures/Figure3/Figure3_v3.pdf',Figure3_v3_png='Schutz_2018_figures/Figure3/Figure3_v3.png')
before<-vapply(file.path(base,files),hash,character(1));names(before)<-names(files)
stopifnot(before[['predictions']]=='6cf426153d4b129aee6cb83f7842f4af09b7de1bce8479e2667b662494c28f77',
 before[['metrics']]=='d8cb6df171f5f3733987b22c81c91d882fbf0753a3799fa836d657d8cdd90122',
 before[['bootstrap']]=='56843ff7e2e9661de9a3ec05a3e120e8bc4185b7a135a4a0a508d63819f5ad36',
 before[['protocol']]=='1420e0b80930256655aa60e145bc8b8ae93823a080ffa4c4a40fe0e64dfcbb0d')
inputs<-c(study='study_order_and_observed_mean.tsv',means='panel_A_B_matched_study_means.tsv',
 assignments='panel_B_study_response_state_assignments.tsv',flows='panel_B_alluvial_study_counts.tsv',
 draws='panel_C_existing_paired_bootstrap_draws.tsv',summary='panel_C_frozen_point_and_CI.tsv',
 gains='panel_D_all_study_gain_loss.tsv',metrics='caption_frozen_study_equal_metrics.tsv',conversion='bootstrap_conversion_audit.json')
input_hash<-vapply(file.path(oldsrc,inputs),hash,character(1));names(input_hash)<-names(inputs)
study<-read(file.path(oldsrc,inputs[['study']]))
means<-read(file.path(oldsrc,inputs[['means']]))
assignments<-read(file.path(oldsrc,inputs[['assignments']]))
flows<-read(file.path(oldsrc,inputs[['flows']]))
draws<-read(file.path(oldsrc,inputs[['draws']]))
cs<-read(file.path(oldsrc,inputs[['summary']]))
d<-read(file.path(oldsrc,inputs[['gains']]))
metrics<-read(file.path(oldsrc,inputs[['metrics']]))
conversion<-jsonlite::read_json(file.path(oldsrc,inputs[['conversion']]),simplifyVector=TRUE)
stopifnot(before[['archive']]==conversion$original_archive_sha256,input_hash[['draws']]==conversion$export_sha256,
 conversion$new_bootstrap_draws==0,conversion$paired_draw_order_preserved,
 nrow(study)==169,sum(study$comparisons)==1723,identical(study$rank,1:169),
 identical(order(study$observed_mean,study$study_id),1:169),nrow(means)==507,nrow(draws)==12000,nrow(d)==507)
models<-c('Ridge','RandomForest','XGBoost');LABELS<-c(Ridge='Ridge',RandomForest='RF',XGBoost='XGB')
MODEL_COLORS<-c(Ridge='#2C607C',RandomForest='#548578',XGBoost='#967B64')
states<-c('Negative','Weak positive','Stronger positive');STATE_COLORS<-setNames(c('#748E9D','#BEC4C8','#B39C7D'),states)
FONT<-'Arial';INK<-'#242A30';MID<-'#66737D';RULE<-'#CBD2D7';OBS<-'#AAB3BA'
for(est in models){
 q<-means[means$estimator==est,]
 stopifnot(nrow(q)==169,identical(q$study_id,study$study_id),identical(q$rank,study$rank),
  max(abs(q$observed_mean-study$observed_mean))<1e-12,all(is.finite(q$predicted_mean)))
}
stopifnot(all(tapply(flows$studies,flows$estimator,sum)==169),nrow(cs)==6,all(table(d$estimator)==169))
write_tsv(study,'shared_A_D_study_order.tsv');write_tsv(means,'panel_A_ranked_response_envelopes.tsv')
write_tsv(assignments,'panel_B_study_state_assignments.tsv');write_tsv(flows,'panel_B_alluvial_study_counts.tsv')
stopifnot(file.copy(file.path(oldsrc,inputs[['draws']]),file.path(src,'panel_C_existing_paired_bootstrap_draws.tsv'),overwrite=TRUE))
write_tsv(metrics,'caption_frozen_study_equal_metrics.tsv')
write_tsv(data.frame(estimator=models,label=LABELS[models],colour=MODEL_COLORS[models]),'model_palette.tsv')
write_tsv(data.frame(state=states,colour=STATE_COLORS[states]),'response_state_palette.tsv')

# C: central quantiles of the existing paired draws, with frozen point and 95% bounds.
intervals<-do.call(rbind,lapply(seq_len(nrow(cs)),function(i){
 q<-draws[draws$estimator==cs$estimator[i]&draws$baseline==cs$baseline[i],]
 stopifnot(nrow(q)==2000,identical(sort(q$draw_id),1:2000))
 qq<-quantile(q$delta_RMSE,c(.025,.10,.25,.75,.90,.975),type=7,names=FALSE)
 stopifnot(max(abs(qq[c(1,6)]-unlist(cs[i,c('ci95_lower','ci95_upper')])))<1e-12)
 data.frame(estimator=cs$estimator[i],baseline=cs$baseline[i],row=cs$row[i],
  coverage=c(.95,.80,.50),lower=c(cs$ci95_lower[i],qq[2],qq[3]),upper=c(cs$ci95_upper[i],qq[5],qq[4]),
  point=cs$point[i],bootstrap_n=2000L,source= c('frozen 95% endpoints','existing-draw quantiles','existing-draw quantiles'))
}))
stopifnot(nrow(intervals)==18,all(is.finite(intervals$lower)),all(is.finite(intervals$upper)))
write_tsv(intervals,'panel_C_nested_interval_blocks.tsv');write_tsv(cs,'panel_C_frozen_point_and_CI.tsv')
bf<-flows[flows$studies>0,];bf$flow_id<-seq_len(nrow(bf))
blodes<-rbind(transform(bf,axis='Observed',state=observed_state),transform(bf,axis='Predicted',state=predicted_state))
blodes$axis<-factor(blodes$axis,levels=c('Observed','Predicted'))
blodes$state<-factor(blodes$state,levels=states)
blodes$estimator<-factor(blodes$estimator,levels=models)
write_tsv(blodes,'panel_B_alluvial_lodes.tsv')

# D: a common signed radial displacement from each model's own zero circle.
d$rank<-study$rank[match(d$study_id,study$study_id)]
d$observed_mean<-study$observed_mean[match(d$study_id,study$study_id)]
wheel_base<-c(Ridge=2.0,RandomForest=3.4,XGBoost=4.8)
radial_scale<-4.0;angular_total<-189;angular_offset<-10;radial_max<-6.3
d$angle_position<-d$rank+angular_offset
d$zero_radius<-wheel_base[d$estimator]
d$tip_radius<-d$zero_radius+radial_scale*d$delta_RMSE
d$inner_radius<-pmin(d$zero_radius,d$tip_radius);d$outer_radius<-pmax(d$zero_radius,d$tip_radius)
d$xmin<-d$angle_position-.40;d$xmax<-d$angle_position+.40
stopifnot(all(d$inner_radius>1.38),all(d$outer_radius<radial_max),
 max(abs((d$tip_radius-d$zero_radius)/radial_scale-d$delta_RMSE))<1e-12)
for(est in models)stopifnot(identical(sort(d$rank[d$estimator==est]),1:169))
# Track ranges must not overlap: use exact extrema, with no clipping or per-model scaling.
for(j in 1:2)stopifnot(max(d$outer_radius[d$estimator==models[j]])<min(d$inner_radius[d$estimator==models[j+1]]))
ring<-study;ring$state<-ifelse(ring$observed_mean<0,states[1],ifelse(ring$observed_mean<log1p(.10),states[2],states[3]))
ring$xmin<-ring$rank+angular_offset-.46;ring$xmax<-ring$rank+angular_offset+.46
ring$inner_radius<-1.10;ring$outer_radius<-1.38
zero_circles<-do.call(rbind,lapply(models,function(est)data.frame(estimator=est,x=seq(10.5,179.5,length.out=600),radius=unname(wheel_base[est]))))
track_labels<-data.frame(estimator=models,x=0,radius=wheel_base[models],label=LABELS[models],colour=MODEL_COLORS[models])
rank_labels<-data.frame(x=c(11,179),radius=6.05,label=c('1','169'))
write_tsv(d[order(match(d$estimator,models),d$rank),],'panel_D_signed_radial_study_bars.tsv')
write_tsv(ring,'panel_D_observed_state_inner_ring.tsv')
write_tsv(data.frame(estimator=models,zero_radius=wheel_base[models],radial_units_per_RMSE=radial_scale), 'panel_D_radial_scale.tsv')
jsonlite::write_json(list(study_order='A order: observed study mean ascending, then study ID',
 clockwise=TRUE,angular_total=angular_total,angular_offset=angular_offset,positive='outward, baseline better',
 negative='inward, model better',all_three_tracks_same_linear_radial_scale=TRUE,radial_units_per_RMSE=radial_scale,
 no_clipping=TRUE,no_per_model_rescaling=TRUE,ring_states=states,
 negative_state='x < 0',weak_state='0 <= x < ln(1.10)',stronger_state='x >= ln(1.10)',state_cutoff=log1p(.10)),
 file.path(src,'display_definition.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)

theme_set(theme_classic(base_size=7,base_family=FONT)+theme(axis.line=element_line(linewidth=.25,colour=INK),
 axis.ticks=element_line(linewidth=.25,colour=INK),axis.text=element_text(size=6.2,colour=INK),axis.title=element_text(size=6.6),
 panel.grid=element_blank(),plot.margin=margin(2,3,2,2),legend.text=element_text(size=5.8),legend.title=element_blank()))
means$estimator<-factor(means$estimator,levels=models)
pA<-ggplot(means,aes(x=rank,group=estimator))+
 geom_ribbon(aes(ymin=pmin(0,observed_mean),ymax=pmax(0,observed_mean)),fill=OBS,alpha=.45,colour=NA)+
 geom_ribbon(aes(ymin=pmin(0,predicted_mean),ymax=pmax(0,predicted_mean),fill=estimator),alpha=.62,colour=NA)+
 geom_hline(yintercept=0,colour='#737E87',linewidth=.28)+
 facet_wrap(~estimator,nrow=1)+scale_fill_manual(values=MODEL_COLORS,guide='none')+
 scale_x_continuous(limits=c(1,169),breaks=c(1,50,100,169),expand=expansion(mult=0))+
 scale_y_continuous(limits=c(-.27,.78),breaks=c(-.2,0,.2,.4,.6),expand=expansion(mult=0))+
 labs(x='Study rank (observed mean ascending)',y='Study-mean logRR')+
 theme(strip.text=element_blank(),strip.background=element_blank(),panel.spacing=unit(4,'mm'))

pB<-ggplot(blodes,aes(x=axis,y=studies,stratum=state,alluvium=flow_id))+
 geom_alluvium(aes(fill=observed_state),width=.12,knot.pos=.45,alpha=.42,colour=NA)+
 geom_stratum(aes(fill=after_stat(stratum)),width=.13,colour='white',linewidth=.3)+
 stat_stratum(geom='text',aes(label=after_stat(ifelse(count>=8,as.character(count),''))),size=2.0,family=FONT,colour=INK)+
 facet_wrap(~estimator,nrow=1)+scale_fill_manual(values=STATE_COLORS,breaks=states,drop=FALSE,name=NULL)+
 scale_x_discrete(limits=c('Observed','Predicted'),expand=c(.12,.08))+
 scale_y_continuous(limits=c(0,169),expand=expansion(mult=0))+
 labs(x=NULL,y=NULL)+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.y=element_blank(),
 axis.text.x=element_text(size=6.2),strip.text=element_blank(),strip.background=element_blank(),panel.spacing=unit(4,'mm'),
 legend.position='bottom',legend.key.height=unit(2,'mm'),legend.key.width=unit(3,'mm'),legend.spacing.x=unit(1,'mm'),
 legend.margin=margin(0,0,0,0),legend.box.spacing=unit(1,'mm'))

intervals$opacity<-c(.15,.30,.55)[match(intervals$coverage,c(.95,.80,.50))]
crow_labels<-paste0(rep(unname(LABELS[models]),each=2),'\n',rep(c('vs train obs mean','vs train study mean'),3))
pC<-ggplot(intervals)+
 geom_rect(aes(xmin=lower,xmax=upper,ymin=row-.27,ymax=row+.27,fill=estimator,alpha=opacity),colour=NA)+
 geom_vline(xintercept=0,colour=INK,linewidth=.40)+
 geom_segment(data=cs,aes(x=point,xend=point,y=row-.33,yend=row+.33),colour=INK,linewidth=.55)+
 scale_fill_manual(values=MODEL_COLORS,guide='none')+scale_alpha_identity(guide='none')+
 scale_x_continuous(limits=c(-.002,.030),breaks=c(0,.01,.02,.03),expand=expansion(mult=0))+
 scale_y_continuous(limits=c(.5,6.5),breaks=6:1,labels=crow_labels,expand=expansion(mult=0))+
 labs(x='Δ study-equal RMSE (model − baseline)',y=NULL)+
 theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_text(size=6.4,lineheight=1.1),axis.title.x=element_text(size=6.4))

pD<-ggplot()+
 geom_rect(data=ring,aes(xmin=xmin,xmax=xmax,ymin=inner_radius,ymax=outer_radius,fill=state),colour=NA)+
 geom_path(data=zero_circles,aes(x=x,y=radius,group=estimator),colour=RULE,linewidth=.18)+
 geom_rect(data=d,aes(xmin=xmin,xmax=xmax,ymin=inner_radius,ymax=outer_radius,fill=estimator),colour=NA,alpha=.86)+
 geom_text(data=track_labels,aes(x=x,y=radius,label=label,colour=colour),size=2.05,family=FONT,fontface='bold')+
 geom_text(data=rank_labels,aes(x=x,y=radius,label=label),size=1.9,family=FONT,colour=MID)+
 annotate('text',x=0,y=0,label='169\nstudies',size=2.1,family=FONT,colour=INK,lineheight=1.05)+
 scale_fill_manual(values=c(MODEL_COLORS,STATE_COLORS),guide='none')+scale_colour_identity(guide='none')+
 scale_x_continuous(limits=c(0,angular_total),expand=expansion(mult=0))+
 scale_y_continuous(limits=c(0,radial_max),expand=expansion(mult=0))+
 coord_polar(theta='x',start=0,direction=1,clip='off')+
 theme_void(base_size=7,base_family=FONT)+theme(plot.margin=margin(0,0,0,0))

txt<-function(label,x,y,size=7,col=INK,face='plain',just='left')grid.text(label,x=unit(x,'npc'),y=unit(y,'npc'),just=just,gp=gpar(fontfamily=FONT,fontsize=size,col=col,fontface=face))
header<-function(tag,title){
 yy<-1-convertHeight(unit(2.4,'mm'),'npc',valueOnly=TRUE)
 xx<-.02+convertWidth(unit(4.4,'mm'),'npc',valueOnly=TRUE)
 txt(tag,.02,yy,9.5,face='bold');txt(title,xx,yy,8.5,face='bold')
}
embed<-function(p,x=.5,y=.45,w=.98,h=.76){pushViewport(viewport(x=x,y=y,width=w,height=h));grid.draw(ggplotGrob(p));popViewport()}
facet_headers<-function(y,centres=c(.205,.535,.865)){for(j in 1:3)txt(LABELS[models[j]],centres[j],y,7.1,col=MODEL_COLORS[j],face='bold',just='centre')}
draw_A<-function(){
 header('A','Ranked response envelope')
 grid.rect(x=unit(.78,'npc'),y=unit(.945,'npc'),width=unit(3,'mm'),height=unit(1.6,'mm'),gp=gpar(fill=OBS,col=NA,alpha=.45))
 txt('Observed',.795,.945,6.0,col=MID)
 facet_headers(.855);embed(pA,y=.415,h=.78)
}
draw_B<-function(){header('B','Response-state alluvial');facet_headers(.885,c(.17,.50,.83));embed(pB,y=.425,h=.82)}
draw_C<-function(){
 header('C','Nested bootstrap intervals')
 for(i in 1:3){
  xx<-c(.075,.28,.475)[i]
  grid.rect(x=unit(xx,'npc'),y=unit(.90,'npc'),width=unit(4,'mm'),height=unit(2.3,'mm'),gp=gpar(fill=INK,col=NA,alpha=c(.15,.30,.55)[i]))
  txt(c('95%','80%','50%')[i],xx+.030,.90,5.9,col=MID)
 }
 grid.lines(x=unit(c(.66,.66),'npc'),y=unit(c(.888,.912),'npc'),gp=gpar(col=INK,lwd=1.1))
 txt('Point estimate',.685,.90,5.9,col=MID)
 embed(pC,y=.425,h=.77)
 txt('Δ < 0: model better',.02,.018,5.8,col=MID)
}
draw_D<-function(){
 header('D','Study generalization wheel')
 embed(pD,y=.50,w=1.04,h=.80)
 txt('Inward: model better',.07,.115,5.9,col=MID)
 txt('Outward: baseline better',.07,.080,5.9,col=MID)
 txt('Inner ring: observed state',.07,.035,5.8,col=MID)
}
capture_panel<-function(fun,w,h){
 tmp<-tempfile(fileext='.pdf');on.exit(unlink(tmp))
 grid.grabExpr(fun(),width=w,height=h,wrap=TRUE,wrap.grobs=TRUE,
  device=function(width,height)grDevices::cairo_pdf(tmp,width=width,height=height,family=FONT))
}
width_mm<-183;height_mm<-218;hA<-60;hB<-60;hBottom<-98
gA<-capture_panel(draw_A,width_mm/25.4,hA/25.4);gB<-capture_panel(draw_B,width_mm/25.4,hB/25.4)
gC<-capture_panel(draw_C,106/25.4,hBottom/25.4);gD<-capture_panel(draw_D,77/25.4,hBottom/25.4)
wrap<-function(g)wrap_elements(full=g,clip=FALSE)
bottom<-(wrap(gC)|wrap(gD))+plot_layout(widths=c(106,77))
fig<-(wrap(gA)/wrap(gB)/bottom)+plot_layout(heights=c(hA,hB,hBottom))
fig<-fig & theme(plot.background=element_rect(fill='white',colour=NA),plot.margin=margin(0,0,0,0))
pdf_file<-file.path(out,'Figure3_v4.pdf');png_file<-file.path(out,'Figure3_v4.png')
grDevices::cairo_pdf(pdf_file,width=width_mm/25.4,height=height_mm/25.4,family=FONT);print(fig);invisible(dev.off())
ragg::agg_png(png_file,width=width_mm,height=height_mm,units='mm',res=600,background='white');print(fig);invisible(dev.off())
ragg::agg_png(file.path(src,'Figure3_v4_QA_preview.png'),width=width_mm,height=height_mm,units='mm',res=180,background='white');print(fig);invisible(dev.off())
for(it in list(list(name='A',plot=wrap(gA),height=hA),list(name='B',plot=wrap(gB),height=hB),list(name='C_D',plot=bottom,height=hBottom))){
 ragg::agg_png(file.path(src,paste0('Figure3_v4_QA_',it$name,'.png')),width=width_mm,height=it$height,units='mm',res=300,background='white')
 print(it$plot);invisible(dev.off())
}
after<-vapply(file.path(base,files),hash,character(1));names(after)<-names(files)
stopifnot(identical(before,after),identical(input_hash,setNames(vapply(file.path(oldsrc,inputs),hash,character(1)),names(inputs))))
con<-file(png_file,'rb');invisible(readBin(con,'raw',8));dims<-ppm<-NULL
repeat{
 len<-readBin(con,'integer',1,size=4,endian='big');if(!length(len))break
 typ<-rawToChar(readBin(con,'raw',4));raw<-readBin(con,'raw',len);readBin(con,'raw',4)
 if(typ=='IHDR')dims<-readBin(raw[1:8],'integer',2,size=4,endian='big')
 if(typ=='pHYs')ppm<-readBin(raw[1:8],'integer',2,size=4,endian='big')
 if(typ=='IEND')break
};close(con)
stopifnot(all(abs(dims-c(width_mm,height_mm)/25.4*600)<2),all(abs(ppm*.0254-600)<.02))
classification<-jsonlite::read_json(file.path(base,files[['classification']]),simplifyVector=TRUE)$A
supported<-sum(vapply(models,function(est)isTRUE(classification[[est]]$supported),logical(1)));stopifnot(supported==0)
qc<-list(task='A only',panels=4,studies=169,comparisons=1723,source_sha256=as.list(before),display_input_sha256=as.list(input_hash),source_unchanged=TRUE,
 new_model_fits=0,new_tests=0,new_bootstrap_draws=0,new_density_estimates=0,
 A='169 original study means in frozen observed-ascending order for each model; zero-based filled areas, not fitted lines',
 B='unchanged independent alluvial study counts, 169 per model',C='central 50/80/95% nested interval blocks of existing paired draws; frozen point and 95% endpoints',
 D='same study order as A clockwise; all 169 values per model; signed bars relative to own zero circle; common radial scale; no clipping',
 D_radial_units_per_RMSE=radial_scale,D_track_ranges_do_not_overlap=TRUE,
 existing_paired_draws_per_contrast=2000,C_frozen_95_CI_and_point=TRUE,visual_review='pending inspection',
 font=FONT,width_mm=width_mm,height_mm=height_mm,PNG_pixels=as.list(dims),PNG_dpi=as.list(ppm*.0254),
 output_sha256=list(pdf=hash(pdf_file),png=hash(png_file)),support_in_caption=paste0(supported,'/3 models supported'),
 libraries=as.list(setNames(vapply(required,function(p)as.character(packageVersion(p)),character(1)),required)))
jsonlite::write_json(qc,file.path(src,'Figure3_v4_QC.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
metric_text<-paste(vapply(c(models,'training_observation_mean','training_study_equal_mean','zero_response'),function(est){
 q<-metrics[metrics$estimator==est,];sprintf('%s: RMSE %.6f, R² %.6f',est,q$RMSE,q$R2)
},character(1)),collapse='; ')
caption<-c('Figure 3 | Response heterogeneity, compressed cross-study predictions and limited gains over mean baselines.',
 'All panels use only the frozen Phase1 Task A 169-fold LOSO results: 1,723 comparisons in 169 studies. RF is RandomForest; XGB is XGBoost. Each study contributes one arithmetic observed mean and one arithmetic OOF predicted mean per model. No model was retrained, no new random draw, statistical test or fitted trend was generated, and no scientific source file or frozen outcome was changed.',
  sprintf('A, Three ranked response envelopes use identical study order: observed study-mean logRR ascending, with study ID breaking ties. Grey is the observed mean, model color is the corresponding OOF predicted mean, and each area is filled to the same zero baseline on a common logRR scale. Filling joins neighboring ranks only and represents no fitted curve or uncertainty band. Predicted values are not sorted independently. The full study range and all extrema are retained. The observed response span is %.6f logRR; the Ridge, RF and XGB predicted spans are %.6f, %.6f and %.6f, respectively.',diff(range(means$observed_mean)),
    diff(range(means$predicted_mean[means$estimator=="Ridge"])),
    diff(range(means$predicted_mean[means$estimator=="RandomForest"])),
    diff(range(means$predicted_mean[means$estimator=="XGBoost"]))),
 'B, Three separate two-axis alluvials show Observed → Ridge, Observed → RF and Observed → XGB study-mean display states. States are negative (<0), weak positive (0≤x<ln(1.10)) and stronger positive (x≥ln(1.10)); exact zero belongs to the middle state. The fixed 10% boundary is descriptive, not a validated biological threshold or replacement label. Each branch sums to 169 studies. Flows are colored by observed state, destinations by predicted state. Counts below 8 remain in the plot but are unlabeled. These are not three stages of a single serial process.',
 'C, Nested central 50%, 80% and 95% intervals of six unchanged paired study-bootstrap ΔRMSE distributions (2,000 saved draws per contrast; original seed 20261006). ΔRMSE = model minus baseline RMSE under study-equal evaluation; negative values favor the model. Darker inner rectangles are the 50% intervals; intermediate and outer rectangles are 80% and 95%. The 50/80% endpoints are descriptive type-7 quantiles of the existing saved draws, not new resampling. The outer 95% endpoints and short vertical point-estimate markers are copied from the frozen Phase1 summary and verified against the same draw archive. The point estimate is the original full-OOF contrast, not a draw median. Intervals are conditional on fixed OOF predictions; the displayed 95% interval is not the multiplicity-adjusted support-gate boundary.',
 'D, A circular display of within-study RMSE(model) minus RMSE(training study-equal mean baseline) for all 169 studies per model, clockwise in exactly A’s order. The inner ring uses B’s observed response categories. Concentric zero circles identify Ridge (inner model track), RF (middle) and XGB (outer). Each colored radial bar extends inward for a negative difference (model better) or outward for a positive difference (baseline better); its radial length is proportional to absolute ΔRMSE. Every track shares one linear radial scale of 4 coordinate units per RMSE unit, with no model-specific normalization, clipping or study selection. The rank labels 1 and 169 delimit the plotting gap. Bar area and angular width are not error magnitude encodings; radial length alone encodes magnitude. The wheel is descriptive heterogeneity, with precise statistical judgment based on C and the frozen metrics.',
 paste0('Frozen full-OOF study-equal metrics: ',metric_text,'. Study-equal RMSE is sqrt(mean within-study MSE), not the mean of study RMSEs. Study-equal R² uses the weighted observed variance, not averaged within-study R².'),
 'The frozen overall gate supports 0/3 models. The original adjusted one-sided paired upper bound at 1−0.05/12, positivity of study-equal R² and other original baseline/MAE checks remain unchanged; no new gate is constructed from the display intervals. The claim concerns these frozen predictors, cohort, models and cross-study protocol, rather than biological efficacy or all possible prediction methods.',
 'A–D source data include every study ID, original mean, fixed state, saved draw ID, interval endpoint and signed radial bar coordinate. Original source hashes and export checks are supplied. All plotting and previews use R; the saved NPZ-to-TSV conversion was non-visual and generated no draws. No density, violin, ridgeline, scatter, heatmap, scorecard, SHAP or feature-importance panel is present. Figure 4 is not generated.')
writeLines(caption,file.path(src,'Figure3_v4_caption.txt'),useBytes=TRUE)
writeLines(c('# Figure 3 v4 source data',
 'shared_A_D_study_order.tsv fixes the identical order in A and D.',
 'panel_A_ranked_response_envelopes.tsv contains 169 observed/predicted study means per model.',
 'panel_B_study_state_assignments.tsv retains the individual study assignments; alluvial counts include zero-count cells.',
 'panel_C_existing_paired_bootstrap_draws.tsv retains 2,000 saved paired draws per contrast. panel_C_nested_interval_blocks.tsv documents 50/80/95% endpoints and their provenance.',
 'panel_D_signed_radial_study_bars.tsv retains all 507 signed differences and bar coordinates; the common scale is in panel_D_radial_scale.tsv.',
 'No raw scientific file was modified. Figure3_v4_caption.txt describes weighting, display states, intervals and radial encoding.',
 'Reproduce with figure3.R in the Figure3 directory, the original frozen Phase1 files and existing v3 source TSVs. All graphics use R.'),file.path(src,'README.md'),useBytes=TRUE)
writeLines(c('Official plotting references:',
 'https://ggplot2.tidyverse.org/reference/geom_ribbon.html',
 'https://ggplot2.tidyverse.org/reference/coord_polar.html',
 'https://corybrunson.github.io/ggalluvial/reference/geom_alluvium.html'),file.path(src,'plotting_documentation.txt'))
capture.output(sessionInfo(),file=file.path(src,'R_sessionInfo.txt'))
message('Figure 3 v4 rendered: ',paste(dims,collapse=' x '),' pixels, 600 dpi; no density, training or resampling.')
