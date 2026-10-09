#!/usr/bin/env Rscript
# Figure 4 v3: visual revision of frozen Task B data only.
# No fitting, resampling, performance recalculation or statistical testing.
invisible(Sys.setlocale('LC_ALL','Chinese (Simplified)_China.utf8'))
sa<-grep('^--file=',commandArgs(FALSE),value=TRUE)
out<-dirname(normalizePath(sub('^--file=','',sa[1]),winslash='/'))
base<-normalizePath(file.path(out,'../..'),winslash='/')
root<-file.path(out,'figure4_plot_data');src<-file.path(root,'v3')
dir.create(src,recursive=TRUE,showWarnings=FALSE)
required<-c('ggplot2','patchwork','ragg','digest','jsonlite','ggraph','tidygraph')
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
draws<-rd('panel_D_existing_paired_bootstrap_draws.tsv')
ds<-rd('panel_D_bootstrap_fraction_and_frozen_point.tsv')
stopifnot(nrow(crops)==40,identical(crops$rank,1:40),sum(crops$comparisons)==1723,
 nrow(contrasts)==120,nrow(draws)==12000,nrow(ds)==6)
# Copies of scientific tables are byte-identical; v3 adds only drawing geometry.
copy_names<-c('crop_support_and_frozen_order.tsv','panel_B_whole_study_purge_anatomy.tsv',
 'panel_B_affected_row_membership.tsv','panel_B_proportional_ribbon_vertices.tsv',
 'panel_C_frozen_cropwise_RMSE_contrasts.tsv','panel_D_existing_paired_bootstrap_draws.tsv',
 'panel_D_bootstrap_fraction_and_frozen_point.tsv','model_palette.tsv','support_palette.tsv','support_summary.tsv')
invisible(file.copy(file.path(root,copy_names),file.path(src,copy_names),overwrite=TRUE))

# A: the graph is used only for deterministic circle packing, not evidence clustering.
# The visual-layout seed has no connection to the frozen bootstrap seed.
nodes<-data.frame(name=c('layout_root',crops$crop),weight=c(0,crops$comparisons),rank=0:40)
graph<-tidygraph::tbl_graph(nodes=nodes,edges=data.frame(from=rep(1,40),to=2:41),directed=TRUE)
set.seed(4003)
layout<-as.data.frame(ggraph::create_layout(graph,layout='circlepack',weight=weight,sort.by=rank))
circles<-layout[layout$leaf,c('name','x','y','r')]
circles<-circles[match(crops$crop,circles$name),];rownames(circles)<-NULL
circles<-cbind(crops,circles[,c('x','y','r')])
# ggraph stores packed radii at single precision; restore exact analytic area.
# A microscopic separation of centres avoids roundoff-only contact overlaps.
circles$r<-sqrt(circles$comparisons/pi);circles$x<-circles$x*1.000001;circles$y<-circles$y*1.000001
circles$area<-pi*circles$r^2
stopifnot(max(abs(circles$area/circles$comparisons-1))<1e-10)
pair_gap<-c()
for(i in 1:39)for(j in (i+1):40)pair_gap<-c(pair_gap,sqrt((circles$x[i]-circles$x[j])^2+(circles$y[i]-circles$y[j])^2)-circles$r[i]-circles$r[j])
stopifnot(min(pair_gap)>-1e-8)
circles$label<-ifelse(circles$comparisons>=20,paste0(gsub(' ','\n',circles$display_label),'\n',circles$study_support,' studies'),'')
circles$fontsize<-ifelse(circles$comparisons>=150,2.70,ifelse(circles$comparisons>=60,2.35,2.05))
vertices<-do.call(rbind,lapply(1:40,function(i){a<-seq(0,2*pi,length.out=257);data.frame(rank=circles$rank[i],support=circles$support[i],x=circles$x[i]+circles$r[i]*cos(a),y=circles$y[i]+circles$r[i]*sin(a))}))
wt(circles,'panel_A_packed_circles.tsv');wt(vertices,'panel_A_circle_vertices.tsv')

# C: purely geometric gaps between frozen support groups; no crop/metric selection.
gap<-1
crops$plot_x<-crops$rank+gap*(match(crops$support,levels_support)-1)
hbands$plot_x<-crops$plot_x[match(hbands$crop,crops$crop)]
hbands$polygon<-paste(hbands$estimator,hbands$band,hbands$sign,hbands$support,sep=':')
regions<-do.call(rbind,lapply(levels_support,function(s){x<-crops$plot_x[crops$support==s];data.frame(support=s,xmin=min(x)-.5,xmax=max(x)+.5,centre=mean(range(x)))}))
wt(crops,'crop_order_with_visual_gaps.tsv');wt(hbands,'panel_C_signed_folded_horizon_bands.tsv');wt(regions,'panel_C_support_background_regions.tsv')

# D: sort all existing draws; never generate draws or refit a model.
sorted<-geometry<-crossings<-vector('list',6)
scale_delta<-.58/max(abs(draws$delta_RMSE)) # same vertical scale for all six rows
halfwidth<-.038
for(i in 1:6){
 q<-draws[draws$estimator==ds$estimator[i]&draws$baseline==ds$baseline[i],]
 q<-q[order(q$delta_RMSE,q$draw_id),];q$sorted_rank<-1:2000;q$rank_fraction<-(q$sorted_rank-.5)/2000
 q$contrast_id<-i;q$row<-7-i;q$y<-q$row+scale_delta*q$delta_RMSE
 q$sign<-ifelse(q$delta_RMSE<0,'Model favoured','Baseline favoured')
 stopifnot(identical(sort(q$draw_id),1:2000),sum(q$delta_RMSE<0)==ds$negative_draws[i],
  sum(q$delta_RMSE<0)/2000==ds$fraction_negative[i])
 k<-sum(q$delta_RMSE<0);stopifnot(k>0,k<2000)
 xzero<-q$rank_fraction[k]+(0-q$delta_RMSE[k])/(q$delta_RMSE[k+1]-q$delta_RMSE[k])*(q$rank_fraction[k+1]-q$rank_fraction[k])
 anchor<-q[k,];anchor$rank_fraction<-xzero;anchor$delta_RMSE<-0;anchor$y<-q$row[1];anchor$draw_id<-NA_integer_;anchor$sorted_rank<-NA_integer_
 gn<-rbind(q[1:k,],anchor);gp<-rbind(anchor,q[(k+1):2000,]);gn$sign<-'Model favoured';gp$sign<-'Baseline favoured'
 gg<-rbind(gn,gp);gg$lower<-gg$y-halfwidth;gg$upper<-gg$y+halfwidth
 gg$fill<-ifelse(gg$sign=='Model favoured',MODEL_COLORS[ds$estimator[i]],'#DDE2E5')
 gg$polygon<-paste(i,gg$sign)
 sorted[[i]]<-q;geometry[[i]]<-gg
 crossings[[i]]<-data.frame(contrast_id=i,estimator=ds$estimator[i],baseline=ds$baseline[i],row=7-i,x_zero=xzero,
  last_negative_draw=q$draw_id[k],first_nonnegative_draw=q$draw_id[k+1],negative_draws=k,
  fraction_negative=ds$fraction_negative[i],frozen_point_delta=ds$point[i])
}
sorted<-do.call(rbind,sorted);geometry<-do.call(rbind,geometry);crossings<-do.call(rbind,crossings)
sorted_export<-sorted
for(nm in names(sorted_export))if(is.double(sorted_export[[nm]]))sorted_export[[nm]]<-sprintf('%.17g',sorted_export[[nm]])
wt(sorted_export,'panel_D_sorted_bootstrap_draws.tsv');wt(geometry,'panel_D_dominance_ribbon_geometry.tsv');wt(crossings,'panel_D_zero_crossings.tsv')
wt(data.frame(vertical_scale=scale_delta,halfwidth=halfwidth,draws_per_row=2000),'panel_D_visual_scale.tsv')

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
pD<-ggplot()+geom_segment(data=crossings,aes(x=0,xend=1,y=row,yend=row),colour='#D7DDE1',linewidth=.18)+
 geom_ribbon(data=geometry,aes(x=rank_fraction,ymin=lower,ymax=upper,fill=fill,group=polygon),colour=NA)+
 geom_segment(data=crossings,aes(x=x_zero,xend=x_zero,y=row-.105,yend=row+.105),colour=INK,linewidth=.28)+
 scale_fill_identity()+scale_x_continuous(limits=c(0,1),expand=expansion(mult=0))+
 scale_y_continuous(limits=c(.78,6.68),expand=expansion(mult=0))+
 theme_void(base_family=FONT)+theme(plot.margin=margin(0,0,0,0))

txt<-function(label,x,y,size=7,col=INK,face='plain',just='left')grid.text(label,x=unit(x,'npc'),y=unit(y,'npc'),just=just,gp=gpar(fontfamily=FONT,fontsize=size,col=col,fontface=face))
header<-function(tag,title){yy<-1-convertHeight(unit(2.4,'mm'),'npc',valueOnly=TRUE);xx<-.02+convertWidth(unit(4.4,'mm'),'npc',valueOnly=TRUE);txt(tag,.02,yy,9.5,face='bold');txt(title,xx,yy,8.5,face='bold')}
embed<-function(p,x=.5,y=.45,w=.98,h=.76){pushViewport(viewport(x=x,y=y,width=w,height=h));grid.draw(ggplotGrob(p));popViewport()}
draw_A<-function(){
 header('A','Crop evidence field');embed(pA,x=.50,y=.475,w=.96,h=.84)
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
 txt('Bands',.59,.970,5.6,col=MID)
 for(i in 1:3){xx<-c(.662,.772,.888)[i];grid.rect(x=unit(xx,'npc'),y=unit(.970,'npc'),width=unit(2.4,'mm'),height=unit(1.6,'mm'),gp=gpar(fill=c('#C8CACB','#888D90','#454E54')[i],col=NA));txt(c('0–.05','.05–.10','.10–.15')[i],xx+.012,.970,5.4,col=MID)}
 for(i in 1:3){xx<-.116+.866*(regions$centre[i]-.5)/42;txt(levels_support[i],xx,.895,6.4,col=MID,face='bold',just='centre')}
 embed(pC,y=.449,w=.985,h=.828)
 txt('Below 0: model better · Above 0: baseline better',.113,.018,6.0,col=MID)
}
draw_D<-function(){
 header('D','Bootstrap dominance ribbons')
 txt('Fraction of bootstrap draws\nfavouring model',.745,.862,6.1,col=MID)
 xD<-.506;yD<-.423;wD<-.465;hD<-.687
 embed(pD,x=xD,y=yD,w=wD,h=hD)
 for(i in 1:6){
  yy<-yD-hD/2+hD*((7-i)-.78)/(6.68-.78)
  txt(LABELS[ds$estimator[i]],.022,yy+.005,6.8,col=MODEL_COLORS[ds$estimator[i]],face='bold')
  txt(ifelse(ds$baseline[i]==baselines[1],'vs train observation mean','vs train study-equal mean'),.022,yy-.033,5.8,col=MID)
  txt(sprintf('%.2f%%',100*ds$fraction_negative[i]),.760,yy+.005,7.0,col=MODEL_COLORS[ds$estimator[i]],face='bold')
  txt(sprintf('Point ΔRMSE %+0.4f',ds$point[i]),.760,yy-.032,5.7,col=MID)
 }
 txt('Sorted paired bootstrap ΔRMSE',.273,.026,5.8,col=MID)
 # Zero crossing is a short line at the actual sign boundary, never an inflated segment.
 grid.lines(x=unit(c(.624,.624),'npc'),y=unit(c(.019,.044),'npc'),gp=gpar(col=INK,lwd=.7));txt('Zero crossing',.636,.032,5.5,col=MID)
}
capture_panel<-function(fun,w,h){tmp<-tempfile(fileext='.pdf');on.exit(unlink(tmp));grid.grabExpr(fun(),width=w,height=h,wrap=TRUE,wrap.grobs=TRUE,device=function(width,height)grDevices::cairo_pdf(tmp,width=width,height=height,family=FONT))}
width_mm<-183;height_mm<-247;hTop<-107;hC<-80;hBottom<-60
gA<-capture_panel(draw_A,110/25.4,hTop/25.4);gB<-capture_panel(draw_B,73/25.4,hTop/25.4)
gC<-capture_panel(draw_C,width_mm/25.4,hC/25.4);gD<-capture_panel(draw_D,width_mm/25.4,hBottom/25.4)
wrap<-function(g)wrap_elements(full=g,clip=FALSE)
top<-(wrap(gA)|wrap(gB))+plot_layout(widths=c(110,73))
fig<-(top/wrap(gC)/wrap(gD))+plot_layout(heights=c(hTop,hC,hBottom))
fig<-fig & theme(plot.background=element_rect(fill='white',colour=NA),plot.margin=margin(0,0,0,0))
pdf_file<-file.path(out,'Figure4_v3.pdf');png_file<-file.path(out,'Figure4_v3.png')
grDevices::cairo_pdf(pdf_file,width=width_mm/25.4,height=height_mm/25.4,family=FONT);print(fig);invisible(dev.off())
ragg::agg_png(png_file,width=width_mm,height=height_mm,units='mm',res=600,background='white');print(fig);invisible(dev.off())
ragg::agg_png(file.path(src,'Figure4_v3_QA_preview.png'),width=width_mm,height=height_mm,units='mm',res=180,background='white');print(fig);invisible(dev.off())
for(it in list(list(name='A_B',plot=top,height=hTop),list(name='C',plot=wrap(gC),height=hC),list(name='D',plot=wrap(gD),height=hBottom))){ragg::agg_png(file.path(src,paste0('Figure4_v3_QA_',it$name,'.png')),width=width_mm,height=it$height,units='mm',res=300,background='white');print(it$plot);invisible(dev.off())}
after<-setNames(vapply(file.path(base,files),hash,character(1)),names(files))
stopifnot(identical(before,after),identical(roothash,vapply(rootfiles,hash,character(1))))
con<-file(png_file,'rb');invisible(readBin(con,'raw',8));dims<-ppm<-NULL
repeat{len<-readBin(con,'integer',1,size=4,endian='big');if(!length(len))break;typ<-rawToChar(readBin(con,'raw',4));raw<-readBin(con,'raw',len);readBin(con,'raw',4);if(typ=='IHDR')dims<-readBin(raw[1:8],'integer',2,size=4,endian='big');if(typ=='pHYs')ppm<-readBin(raw[1:8],'integer',2,size=4,endian='big');if(typ=='IEND')break};close(con)
stopifnot(all(abs(dims-c(width_mm,height_mm)/25.4*600)<2),all(abs(ppm*.0254-600)<.02))
caption<-c('Figure 4 | Cross-crop transfer under whole-study purge.',
 'All panels use frozen Phase1 Task B results (1,723 comparisons, 169 studies, 40 crops), with no retraining, new bootstrap draws or statistical tests. Crop support/order, model colors, individual crop RMSE contrasts and global bootstrap estimates are unchanged.',
 'A, Packed-circle evidence field. Every crop is retained and analytic circle area is exactly proportional to comparison count. Color identifies the inherited support tiers: high support (≥3 studies and ≥20 comparisons), low support (all other nonsingleton crops), and singleton (one study); the tiers contain 15/11/14 crops and 1,492/164/67 comparisons. Large circles show crop name and study support. Small crops remain at their correct area without forced labels; all 40 crop names and counts are in the source data. Packing takes the same fixed crop-rank metadata as B/C; spatial positions are visual layout only and do not encode quantitative axes, taxonomy or experiment independence. Study counts across crops are not additive.',
 'B, Forty independent tapered ribbons, in the same frozen crop-support order as C and the A layout input. Total horizontal width is the affected row count (test + purge), divided into held-out crop test rows and other-crop rows from those same studies that must be purged. Width, not polygon area, represents count; the faint vertical line provides a common origin. Train/test and train/purge study overlap are zero. Counts across folds can repeat purged rows/studies. Full membership IDs are provided.',
 'C, Frozen crop-wise study-equal ΔRMSE = RMSE(model) − RMSE(training study-equal mean baseline), with all 40 crops ordered by unchanged study support and original tie ordering. Rank is not chronological or selected by performance. Small gaps separate the three inherited support tiers. Three common 0.05-wide folded magnitude bands (0–0.05, 0.05–0.10, 0.10–0.15) overlay progressively darker contributions; negative is below zero (model better), positive above zero (baseline better). The ordinate is folded within-band height, not untransformed ΔRMSE. Full signed values remain in source data. Neighboring crop ranks are connected only for the envelope; no fitted trend or crop-level significance test is implied.',
 'D, Continuous sorted paired-bootstrap ribbons for three models against two training-mean baselines. Each ribbon uses all 2,000 saved global Task B study-equal ΔRMSE draws, ordered by value with original draw IDs retained. Its vertical displacement from a faint zero reference encodes the sorted signed values at a common scale; constant narrow thickness is a display device. Model color denotes ΔRMSE<0 and pale grey ΔRMSE≥0. A short line marks the linearly interpolated zero crossing between adjacent opposite-sign draws. The exact fraction of draws favouring the model (including 0.10% and 0.05% for Ridge) and the original full-OOF point ΔRMSE are printed at right; small negative portions are not enlarged. Frequencies are descriptive fractions of saved conditional paired study-bootstrap draws, not p-values, posterior probabilities or significance probabilities. Bootstrap seed 20261006, fixed predictions; no models were refitted. These draws describe the global study-equal crop-out OOF metric, not an average of 40 crop RMSEs or separate crop-specific bootstraps.',
 'All plots, exports and previews use R; frozen scientific files and Figures 1–3 are unchanged. Source data include exact circle areas, complete crop identifiers, original split membership, unchanged contrasts, visual gap coordinates, sorted bootstrap draw IDs and zero-crossing geometry.')
writeLines(caption,file.path(src,'Figure4_v3_caption.txt'),useBytes=TRUE)
writeLines(c('# Figure 4 v3 source data','',
 'Frozen scientific source tables were copied byte-for-byte from ../. No outcome, prediction, metric or draw was changed.',
 'A: panel_A_packed_circles.tsv (40 crops; analytic area = comparisons), circle vertices (257 per crop). No minimum radius. Positions from ggraph circle packing, visual seed 4003.',
 'B: original purge counts, ribbons and full split membership. Original crop rank is unchanged.',
 'C: original 120 contrasts; horizon bands preserve all original lower/upper/fill values. Only plot_x gaps and polygon grouping by support are added.',
 'D: original draws TSV (byte-identical); sorted_draws includes all 12,000 original paired IDs. Ribbon_geometry adds interpolated zero anchors (not new bootstrap draws). Zero crossings and exact fractions are stored separately.',
 'Reproduce with ../../figure4.R. R performs packing and plotting only. The root source folder preserves the earlier draft for provenance.',
 'Caption, visual contract, QC, session info and R-generated previews accompany the source files.'),file.path(src,'README.md'),useBytes=TRUE)
writeLines(c('https://ggraph.data-imaginist.com/reference/layout_tbl_graph_circlepack.html','https://ggplot2.tidyverse.org/reference/geom_polygon.html','https://ggplot2.tidyverse.org/reference/geom_ribbon.html'),file.path(src,'plotting_documentation.txt'))
capture.output(sessionInfo(),file=file.path(src,'R_sessionInfo.txt'))
qc<-list(status='RENDERED_PENDING_VISUAL_QC',source_sha256=as.list(before),root_plot_sources_unchanged=TRUE,
 panels=4,crops=40,comparisons=1723,studies=169,exact_circle_area=TRUE,minimum_circle_gap=min(pair_gap),
 crop_order_unchanged=TRUE,horizon_contrasts_unchanged=TRUE,horizon_group_gaps=gap,
 paired_draws_per_contrast=2000,paired_draw_ids_retained=TRUE,frequencies_exact=TRUE,
 models_fitted=0,new_bootstrap_draws=0,new_statistical_tests=0,
 font=FONT,width_mm=width_mm,height_mm=height_mm,PNG_pixels=as.list(dims),PNG_dpi=as.list(ppm*.0254),
 output_sha256=list(pdf=hash(pdf_file),png=hash(png_file)),script_sha256=hash(file.path(out,'figure4.R')))
jsonlite::write_json(qc,file.path(src,'figure4_v3_QC.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
message('Figure 4 v3 rendered: ',paste(dims,collapse=' x '),' pixels; 600 dpi; frozen results unchanged.')
