#!/usr/bin/env Rscript
# Figure 4 final: visual revision of frozen Task B data only.
# No fitting, resampling, performance recalculation or statistical testing.
invisible(Sys.setlocale('LC_ALL','Chinese (Simplified)_China.utf8'))
sa<-grep('^--file=',commandArgs(FALSE),value=TRUE)
out<-dirname(normalizePath(sub('^--file=','',sa[1]),winslash='/'))
base<-normalizePath(file.path(out,'../..'),winslash='/')
root<-file.path(out,'figure4_plot_data');v3<-file.path(root,'v3');src<-file.path(root,'final')
dir.create(src,recursive=TRUE,showWarnings=FALSE)
required<-c('ggplot2','patchwork','ragg','digest','jsonlite')
stopifnot(all(vapply(required,requireNamespace,logical(1),quietly=TRUE)))
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(grid)})
read<-function(p)read.delim(p,check.names=FALSE,stringsAsFactors=FALSE)
rd<-function(f)read(file.path(root,f))
wt<-function(d,f)write.table(d,file.path(src,f),sep='\t',quote=TRUE,row.names=FALSE,na='NA',fileEncoding='UTF-8')
hash<-function(p)digest::digest(file=p,algo='sha256',serialize=FALSE)
oldqc<-jsonlite::fromJSON(file.path(root,'figure4_QC.json'))
files<-c(matrix='Schutz_2018_phase0/analysis_matrix.tsv',splits='Schutz_2018_phase0/split_manifest.tsv',
 crops='Schutz_2018_figures/Figure1/figure1_plot_data/crop_annotations_and_order.tsv',
 predictions='Schutz_2018_phase1/phase1_predictions.tsv',metrics='Schutz_2018_phase1/phase1_metrics.tsv',
 bootstrap='Schutz_2018_phase1/phase1_bootstrap.tsv',archive='Schutz_2018_phase1/run/bootstrap_draws.npz',
 protocol='Schutz_2018_phase1/evidence/protocol.json',
 Figure1_pdf='Schutz_2018_figures/Figure1/Figure1_v2.pdf',Figure1_png='Schutz_2018_figures/Figure1/Figure1_v2.png',
 Figure2_pdf='Schutz_2018_figures/Figure2/Figure2_final.pdf',Figure2_png='Schutz_2018_figures/Figure2/Figure2_final.png',
 Figure3_pdf='Schutz_2018_figures/Figure3/Figure3_final.pdf',Figure3_png='Schutz_2018_figures/Figure3/Figure3_final.png',
 Figure3_script='Schutz_2018_figures/Figure3/figure3.R',Figure3_state='Schutz_2018_figures/Figure3/FIGURE3_FINAL_STATE.md')
before<-setNames(vapply(file.path(base,files),hash,character(1)),names(files))
stopifnot(identical(unname(before),unname(unlist(oldqc$source_sha256[names(files)]))))
rootfiles<-list.files(root,pattern='\\.(tsv|json)$',full.names=TRUE)
roothash<-vapply(rootfiles,hash,character(1))
v3files<-list.files(v3,full.names=TRUE);v3hash<-vapply(v3files,hash,character(1))
oldfigs<-file.path(out,c('Figure4_v3.pdf','Figure4_v3.png','figure4_v3.R'));oldfigsha<-vapply(oldfigs,hash,character(1))
models<-c('Ridge','RandomForest','XGBoost')
baselines<-c('training_observation_mean','training_study_equal_mean')
LABELS<-c(Ridge='Ridge',RandomForest='RF',XGBoost='XGBoost')
MODEL_COLORS<-c(Ridge='#2C607C',RandomForest='#548578',XGBoost='#967B64')
levels_support<-c('High support','Low support','Singleton')
SUPPORT_COLORS<-setNames(c('#CAD5DD','#DDD5C6','#E2E5E6'),levels_support)
ROLE_COLORS<-c(Test='#59798F',Purge='#B4A18B')
FONT<-'Arial';INK<-'#242A30';MID<-'#66737D';RULE<-'#CBD2D7'
crops<-rd('crop_support_and_frozen_order.tsv');crops<-crops[order(crops$rank),]
anatomy<-rd('panel_B_whole_study_purge_anatomy.tsv')
ribbons<-rd('panel_B_proportional_ribbon_vertices.tsv')
contrasts<-rd('panel_C_frozen_cropwise_RMSE_contrasts.tsv')
hbands<-rd('panel_C_signed_folded_horizon_bands.tsv')
stopifnot(nrow(crops)==40,identical(crops$rank,1:40),sum(crops$comparisons)==1723,
 nrow(contrasts)==120)
# Scientific tables are byte-identical. Final D counts fixed crop-fold signs only.
copy_names<-c('crop_support_and_frozen_order.tsv','panel_B_whole_study_purge_anatomy.tsv',
 'panel_B_affected_row_membership.tsv','panel_B_proportional_ribbon_vertices.tsv',
 'panel_C_frozen_cropwise_RMSE_contrasts.tsv','model_palette.tsv','support_palette.tsv','support_summary.tsv')
invisible(file.copy(file.path(root,copy_names),file.path(src,copy_names),overwrite=TRUE))

# A/C: reuse the saved v3 coordinates byte-for-byte, with no repacking.
abcfiles<-c('panel_A_packed_circles.tsv','panel_A_circle_vertices.tsv',
 'crop_order_with_visual_gaps.tsv','panel_C_signed_folded_horizon_bands.tsv','panel_C_support_background_regions.tsv')
invisible(file.copy(file.path(v3,abcfiles),file.path(src,abcfiles),overwrite=TRUE))
circles<-read(file.path(v3,abcfiles[1]));vertices<-read(file.path(v3,abcfiles[2]))
crops<-read(file.path(v3,abcfiles[3]));hbands<-read(file.path(v3,abcfiles[4]));regions<-read(file.path(v3,abcfiles[5]))

# D: exact support-tier Marimekko areas count crop folds, not observations.
mosaic<-list();tiercounts<-list();cropstate<-contrasts
cropstate$response_state<-ifelse(cropstate$delta_RMSE<0,'Model better','Baseline better')
wt(cropstate,'panel_D_frozen_crop_fold_classification.tsv')
for(m in models){
 cursor<-0
 for(s in levels_support){
  q<-cropstate[cropstate$estimator==m&cropstate$support==s,]
  total<-nrow(q);better<-sum(q$delta_RMSE<0);other<-total-better
  width<-total/40;boundary<-better/total
  tiercounts[[length(tiercounts)+1]]<-data.frame(estimator=m,support=s,crop_folds=total,
   model_better=better,baseline_better=other,fraction_model_better=boundary)
  mosaic[[length(mosaic)+1]]<-data.frame(estimator=m,support=s,
   response_state=c('Model better','Baseline better'),n=c(better,other),tier_total=total,
   xmin=cursor,xmax=cursor+width,ymin=c(0,boundary),ymax=c(boundary,1),
   fill=c(MODEL_COLORS[m],'#E1E5E7'),text_colour=c('white',INK))
  cursor<-cursor+width
 }
 stopifnot(abs(cursor-1)<1e-12)
}
mosaic<-do.call(rbind,mosaic);tiercounts<-do.call(rbind,tiercounts)
mosaic$area<-(mosaic$xmax-mosaic$xmin)*(mosaic$ymax-mosaic$ymin)
mosaic$label<-paste0(mosaic$n,'/',mosaic$tier_total)
stopifnot(max(abs(mosaic$area-mosaic$n/40))<1e-12,all(mosaic$n>0))
wt(mosaic,'panel_D_support_tier_mosaic_rectangles.tsv');wt(tiercounts,'panel_D_support_tier_fold_counts.tsv')
theme_set(theme_classic(base_size=7,base_family=FONT)+theme(axis.line=element_line(linewidth=.25,colour=INK),
 axis.ticks=element_line(linewidth=.25,colour=INK),axis.text=element_text(size=5.6,colour=INK),axis.title=element_text(size=6.5),
 panel.grid=element_blank(),plot.margin=margin(2,3,2,2),legend.position='none'))
pA<-ggplot(vertices,aes(x,y,group=rank,fill=support))+geom_polygon(colour='white',linewidth=.30)+
 geom_text(data=circles,aes(x,y,label=label,size=fontsize),inherit.aes=FALSE,family=FONT,colour=INK,lineheight=1.03)+
 scale_fill_manual(values=SUPPORT_COLORS)+scale_size_identity()+coord_fixed()+
 theme_void(base_family=FONT)+theme(plot.margin=margin(0,0,0,0),legend.position='none')
pB<-ggplot(ribbons,aes(x,y,group=polygon,fill=role))+
 geom_vline(xintercept=0,colour='#DDE3E6',linewidth=.22)+geom_polygon(colour=NA)+
 scale_fill_manual(values=ROLE_COLORS)+
 scale_x_continuous(limits=c(0,max(anatomy$affected_n)*1.04),breaks=function(lim)pretty(lim,n=3),expand=expansion(mult=0))+
 scale_y_continuous(limits=c(.5,40.5),breaks=40:1,labels=crops$display_label,expand=expansion(mult=0))+
 labs(x='Affected comparison rows',y=NULL)+
 theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_text(size=5.6),axis.title.x=element_text(size=6.2))
hbands$estimator<-factor(hbands$estimator,levels=models)
pC<-ggplot()+geom_rect(data=transform(regions,fill=SUPPORT_COLORS[support]),
 aes(xmin=xmin,xmax=xmax,ymin=-.05,ymax=.05,fill=fill),alpha=.28,colour=NA)
for(k in 1:3)pC<-pC+geom_ribbon(data=hbands[hbands$band==k,],aes(x=plot_x,ymin=lower,ymax=upper,group=polygon,fill=fill),colour=NA)
pC<-pC+geom_hline(yintercept=0,colour='#58656E',linewidth=.32)+
 scale_fill_identity()+facet_grid(estimator~.,switch='y',labeller=labeller(estimator=LABELS))+
 scale_x_continuous(limits=c(.5,42.5),breaks=crops$plot_x,labels=crops$display_label,expand=expansion(mult=0))+
 scale_y_continuous(limits=c(-.05,.05),breaks=c(-.05,0,.05),labels=c('−.05','0','.05'),expand=expansion(mult=0))+
 labs(x=NULL,y='Folded ΔRMSE')+
 theme(axis.line.y=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=6),
 axis.text.y=element_text(size=5.3),axis.title.y=element_text(size=6.5),panel.spacing.y=unit(3,'mm'),
 strip.background=element_blank(),strip.text.y.left=element_text(angle=0,size=6.5,face='bold'),strip.placement='outside')
pD<-lapply(models,function(m){
 z<-mosaic[mosaic$estimator==m,]
 ggplot(z)+geom_rect(aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=fill),colour='white',linewidth=.45)+
  geom_text(aes(x=(xmin+xmax)/2,y=(ymin+ymax)/2,label=label,colour=text_colour),family=FONT,size=2.35)+
  scale_fill_identity()+scale_colour_identity()+
  scale_x_continuous(limits=c(0,1),expand=expansion(mult=0))+
  scale_y_continuous(limits=c(0,1),expand=expansion(mult=0))+
  theme_void(base_family=FONT)+theme(plot.margin=margin(0,0,0,0),legend.position='none')
})
txt<-function(label,x,y,size=7,col=INK,face='plain',just='left')grid.text(label,x=unit(x,'npc'),y=unit(y,'npc'),just=just,gp=gpar(fontfamily=FONT,fontsize=size,col=col,fontface=face))
header<-function(tag,title){yy<-1-convertHeight(unit(2.4,'mm'),'npc',valueOnly=TRUE);xx<-.02+convertWidth(unit(4.4,'mm'),'npc',valueOnly=TRUE);txt(tag,.02,yy,9.5,face='bold');txt(title,xx,yy,8.5,face='bold')}
embed<-function(p,x=.5,y=.45,w=.98,h=.76){pushViewport(viewport(x=x,y=y,width=w,height=h));grid.draw(ggplotGrob(p));popViewport()}
draw_A<-function(){
 header('A','Crop evidence field');txt('Circle area = comparisons; labels = number of studies.',.022,.928,5.8,col=MID);embed(pA,x=.50,y=.475,w=.96,h=.84)
 for(i in 1:3){xx<-c(.06,.40,.70)[i];grid.circle(x=unit(xx,'npc'),y=unit(.042,'npc'),r=unit(1.1,'mm'),gp=gpar(fill=SUPPORT_COLORS[i],col=NA));txt(levels_support[i],xx+.021,.042,5.8,col=MID)}
}
draw_B<-function(){
 header('B','Whole-study purge anatomy')
 for(i in 1:2){xx<-c(.05,.40)[i];grid.rect(x=unit(xx,'npc'),y=unit(.927,'npc'),width=unit(3,'mm'),height=unit(1.8,'mm'),gp=gpar(fill=ROLE_COLORS[i],col=NA));txt(c('Test crop','Other-crop purge')[i],xx+.027,.927,5.8,col=MID)}
 embed(pB,y=.442,h=.87)
}
draw_C<-function(){
 header('C','Crop-wise transfer horizon')
 # Compact magnitude legend alongside the heading, not a footer or dashboard.
 txt('|ΔRMSE| bands',.53,.970,5.6,col=MID)
 for(i in 1:3){xx<-c(.662,.772,.888)[i];grid.rect(x=unit(xx,'npc'),y=unit(.970,'npc'),width=unit(2.4,'mm'),height=unit(1.6,'mm'),gp=gpar(fill=c('#C8CACB','#888D90','#454E54')[i],col=NA));txt(c('0–.05','.05–.10','.10–.15')[i],xx+.012,.970,5.4,col=MID)}
 for(i in 1:3){xx<-.116+.866*(regions$centre[i]-.5)/42;txt(levels_support[i],xx,.895,6.4,col=MID,face='bold',just='centre')}
 embed(pC,y=.449,w=.985,h=.828)
 txt('Below 0: model better · Above 0: baseline better',.113,.018,6.0,col=MID)
}
draw_D<-function(){
 header('D','Support-tier transfer mosaic')
 txt('Model colour: Model better · Grey: Baseline better',.555,.960,5.8,col=MID)
 centres<-c(1/6,.5,5/6);w<-.302
 for(i in 1:3){
  x<-centres[i];txt(LABELS[models[i]],x,.847,7.2,col=MODEL_COLORS[i],face='bold',just='centre')
  embed(pD[[i]],x=x,y=.489,w=w,h=.553)
  t<-mosaic[mosaic$estimator==models[i]&mosaic$response_state=='Model better',]
  for(j in 1:3)txt(levels_support[j],x-w/2+w*(t$xmin[j]+t$xmax[j])/2,.160,5.8,col=MID,just='centre')
 }
 txt('Labels: crop folds / tier total',.022,.061,5.8,col=MID)
}
capture_panel<-function(fun,w,h){tmp<-tempfile(fileext='.pdf');on.exit(unlink(tmp));grid.grabExpr(fun(),width=w,height=h,wrap=TRUE,wrap.grobs=TRUE,device=function(width,height)grDevices::cairo_pdf(tmp,width=width,height=height,family=FONT))}
width_mm<-183;height_mm<-247;hTop<-107;hC<-80;hBottom<-60
gA<-capture_panel(draw_A,110/25.4,hTop/25.4);gB<-capture_panel(draw_B,73/25.4,hTop/25.4)
gC<-capture_panel(draw_C,width_mm/25.4,hC/25.4);gD<-capture_panel(draw_D,width_mm/25.4,hBottom/25.4)
wrap<-function(g)wrap_elements(full=g,clip=FALSE)
top<-(wrap(gA)|wrap(gB))+plot_layout(widths=c(110,73))
fig<-(top/wrap(gC)/wrap(gD))+plot_layout(heights=c(hTop,hC,hBottom))
fig<-fig & theme(plot.background=element_rect(fill='white',colour=NA),plot.margin=margin(0,0,0,0))
pdf_file<-file.path(out,'Figure4_final.pdf');png_file<-file.path(out,'Figure4_final.png')
grDevices::cairo_pdf(pdf_file,width=width_mm/25.4,height=height_mm/25.4,family=FONT);print(fig);invisible(dev.off())
ragg::agg_png(png_file,width=width_mm,height=height_mm,units='mm',res=600,background='white');print(fig);invisible(dev.off())
ragg::agg_png(file.path(src,'Figure4_final_QA_preview.png'),width=width_mm,height=height_mm,units='mm',res=180,background='white');print(fig);invisible(dev.off())
for(it in list(list(name='A_B',plot=top,height=hTop),list(name='C',plot=wrap(gC),height=hC),list(name='D',plot=wrap(gD),height=hBottom))){ragg::agg_png(file.path(src,paste0('Figure4_final_QA_',it$name,'.png')),width=width_mm,height=it$height,units='mm',res=300,background='white');print(it$plot);invisible(dev.off())}
after<-setNames(vapply(file.path(base,files),hash,character(1)),names(files))
stopifnot(identical(before,after),identical(roothash,vapply(rootfiles,hash,character(1))),identical(v3hash,vapply(v3files,hash,character(1))),identical(oldfigsha,vapply(oldfigs,hash,character(1))))
con<-file(png_file,'rb');invisible(readBin(con,'raw',8));dims<-ppm<-NULL
repeat{len<-readBin(con,'integer',1,size=4,endian='big');if(!length(len))break;typ<-rawToChar(readBin(con,'raw',4));raw<-readBin(con,'raw',len);readBin(con,'raw',4);if(typ=='IHDR')dims<-readBin(raw[1:8],'integer',2,size=4,endian='big');if(typ=='pHYs')ppm<-readBin(raw[1:8],'integer',2,size=4,endian='big');if(typ=='IEND')break};close(con)
stopifnot(all(abs(dims-c(width_mm,height_mm)/25.4*600)<2),all(abs(ppm*.0254-600)<.02))
caption<-c('Figure 4 | Cross-crop transfer under whole-study purge.',
 'All panels use frozen Phase1 Task B results (1,723 comparisons, 169 studies, 40 crops), with no retraining, new bootstrap draws or statistical tests. A/B/C scientific values, geometry and crop order are unchanged from v3.',
 'A, Crop evidence bubble field. Circle area = comparisons; labels = number of studies. Every crop is retained with circle area exactly proportional to its comparison count. Large circles additionally show crop names. The inherited support tiers are high support (≥3 studies and ≥20 comparisons), low support (other nonsingleton crops), and singleton (one study), containing 15/11/14 crops and 1,492/164/67 comparisons. Small crops retain their correct area; complete names and counts are in source data. Positions are visual layout only. Studies can support multiple crops, so study counts across crops are not additive.',
 'B, Forty independent tapered ribbons in the frozen support order. Total horizontal width encodes affected rows (test + purge); segments encode held-out crop test rows and other-crop rows from the same studies purged from training. Width, not polygon area, encodes count. The faint line gives a common origin. Train/test and train/purge study overlaps are zero. Purged rows/studies may recur across crop folds.',
 'C, Crop-wise study-equal ΔRMSE = RMSE(model) − RMSE(training study-equal mean baseline), retaining all 40 crops in the unchanged support order and original tie ordering. Rank is not chronological or selected by performance. Gaps separate the support tiers. The |ΔRMSE| bands are three common 0.05-wide magnitude bands (0–0.05, 0.05–0.10, 0.10–0.15), overlaid progressively darker, with negative contributions below zero (model better) and positive above (baseline better). The ordinate is folded within-band height, not untransformed ΔRMSE. Full signed contrasts are retained in source data. Connections between adjacent ranks form the envelope only; there is no fitted trend or crop-level significance test.',
 'D, Support-tier transfer Marimekko plots for Ridge, RF and XGBoost, derived by counting the same frozen crop-level contrasts used in C. Within each model, horizontal widths are proportional to tier crop counts (15/40, 11/40, 14/40). Each tier is split vertically into model better (ΔRMSE<0; model color, below) and baseline better (ΔRMSE≥0; pale grey, above). Ties are included in the latter category. Rectangle area encodes the number of crop folds divided by 40. Labels are n/tier total. Model-better counts for high/low/singleton tiers are Ridge 4/15, 6/11, 6/14; RF 6/15, 6/11, 6/14; XGBoost 3/15, 7/11, 7/14. Counts describe local advantages by evidence tier and are not probabilities, significance tests or independent experiment counts. Crop folds can share studies, and a majority of favorable folds would not by itself establish reliable generalization. No bootstrap panel is included in Figure 4.',
 'All graphics, exports and previews use R. Scientific input files, v3 source files and Figures 1–3 remain unchanged. No Figure 5 work is performed.')
writeLines(caption,file.path(src,'Figure4_final_caption.txt'),useBytes=TRUE)
writeLines(c('# Figure 4 final source data','',
 'A/B/C: copied scientific/geometry TSVs are identical to v3; no repacking or new horizon values.',
 'D: panel_D_frozen_crop_fold_classification.tsv retains every frozen model/crop contrast and adds its sign category.',
 'panel_D_support_tier_fold_counts.tsv has all nine model/tier counts. panel_D_support_tier_mosaic_rectangles.tsv has all 18 exact count-proportional rectangles.',
 'Area = crop-fold count / 40; width = tier crop count / 40; height = category count / tier count. No statistical tests or model fitting.',
 'Run ../../figure4.R from this source directory for R-only reproduction. Caption, session info, QC and previews included.'),file.path(src,'README.md'),useBytes=TRUE)
capture.output(sessionInfo(),file=file.path(src,'R_sessionInfo.txt'))
qc<-list(status='RENDERED_PENDING_VISUAL_QC',source_sha256=as.list(before),root_plot_sources_unchanged=TRUE,
 v3_files_unchanged=TRUE,v3_artifacts_sha256=as.list(oldfigsha),
 panels=4,crops=40,comparisons=1723,studies=169,ABC_scientific_geometry_unchanged=TRUE,
 mosaic_models=3,tier_crop_counts=c(15,11,14),mosaic_rectangles=18,area_equals_crop_folds_divided_by_40=TRUE,
 model_better_definition='delta_RMSE < 0',baseline_better_definition='delta_RMSE >= 0',
 models_fitted=0,new_bootstrap_draws=0,new_statistical_tests=0,bootstrap_panel_removed=TRUE,
 font=FONT,width_mm=width_mm,height_mm=height_mm,PNG_pixels=as.list(dims),PNG_dpi=as.list(ppm*.0254),
 output_sha256=list(pdf=hash(pdf_file),png=hash(png_file)),script_sha256=hash(file.path(out,'figure4.R')))
jsonlite::write_json(qc,file.path(src,'figure4_final_QC.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
message('Figure 4 final rendered: ',paste(dims,collapse=' x '),' pixels; 600 dpi; frozen results unchanged.')
