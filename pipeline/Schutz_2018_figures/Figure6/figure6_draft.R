#!/usr/bin/env Rscript
# Plot frozen results only. No model training or inferential analysis.
invisible(Sys.setlocale('LC_ALL','Chinese (Simplified)_China.utf8'))
sa<-grep('^--file=',commandArgs(FALSE),value=TRUE)
out<-dirname(normalizePath(sub('^--file=','',sa[1]),winslash='/'))
base<-normalizePath(file.path(out,'../..'),winslash='/')
src<-file.path(out,'figure6_plot_data');dir.create(src,showWarnings=FALSE)
.libPaths(c(file.path(base,'Schutz_2018_figures/R_lib'),.libPaths()))
pk<-c('grid','patchwork','ggplot2','ragg','digest','jsonlite')
stopifnot(all(vapply(pk,requireNamespace,logical(1),quietly=TRUE)))
suppressPackageStartupMessages({library(grid);library(patchwork);library(ggplot2)})
read<-function(p)read.delim(p,check.names=FALSE,stringsAsFactors=FALSE,na.strings='NA')
wt<-function(d,f)write.table(d,file.path(src,f),sep='\t',row.names=FALSE,quote=TRUE,na='NA',fileEncoding='UTF-8')
hash<-function(p)digest::digest(file=p,algo='sha256',serialize=FALSE)
files<-c(phase1_metrics='Schutz_2018_phase1/phase1_metrics.tsv',
 phase2_metrics='Schutz_2018_phase2/phase2_sensitivity_metrics.tsv',
 external_predictions='Schutz_2018_external_Vermeersch_validation/external_predictions.tsv',
 external_metrics='Schutz_2018_external_Vermeersch_validation/external_metrics.tsv',
 external_QC='Schutz_2018_external_Vermeersch_validation/evidence/final_qc.json',
 external_summary='Schutz_2018_external_Vermeersch_validation/evidence/evaluation_summary.json',
 dictionary='Schutz_2018_phase0/data_dictionary.tsv')
protected<-c(file.path(base,files),file.path(base,c('Schutz_2018_phase0/analysis_matrix.tsv',
 'Schutz_2018_phase0/split_manifest.tsv','Schutz_2018_phase1/phase1_predictions.tsv',
 'Schutz_2018_phase1/phase1_bootstrap.tsv',
 'Schutz_2018_figures/Figure1/Figure1_v2.pdf','Schutz_2018_figures/Figure1/Figure1_v2.png',
 'Schutz_2018_figures/Figure2/Figure2_final.pdf','Schutz_2018_figures/Figure2/Figure2_final.png',
 'Schutz_2018_figures/Figure3/Figure3_final.pdf','Schutz_2018_figures/Figure3/Figure3_final.png',
 'Schutz_2018_figures/Figure4/Figure4_final.pdf','Schutz_2018_figures/Figure4/Figure4_final.png',
 'Schutz_2018_figures/Figure5/Figure5_v2.pdf','Schutz_2018_figures/Figure5/Figure5_v2.png')))
before<-vapply(protected,hash,character(1))
rd<-function(k)read(file.path(base,files[[k]]))
m1<-rd('phase1_metrics');m2<-rd('phase2_metrics');ext<-rd('external_predictions');em<-rd('external_metrics')
eqc<-jsonlite::read_json(file.path(base,files[['external_QC']]),simplifyVector=TRUE)
stopifnot(hash(file.path(base,files[['external_predictions']]))==eqc$output_sha256[['external_predictions.tsv']],
 hash(file.path(base,files[['external_metrics']]))==eqc$output_sha256[['external_metrics.tsv']])
models<-c('Ridge','RandomForest','XGBoost')
baselines<-c('training_observation_mean','training_study_equal_mean')
methods<-c(models,baselines,'zero_response')
LABELS<-c(Ridge='Ridge',RandomForest='RF',XGBoost='XGBoost',
 training_observation_mean='Train observation mean',training_study_equal_mean='Train study-equal mean',zero_response='Zero')
MODEL_COLORS<-c(Ridge='#2C607C',RandomForest='#548578',XGBoost='#967B64')
FONT<-'Arial';INK<-'#242A30';MID<-'#66737D';RULE<-'#CBD2D7'
wt(data.frame(method=methods,label=LABELS[methods],colour=c(MODEL_COLORS,rep('#879198',3))),'model_palette.tsv')

# A: select study-equal primary RMSE; contrast already evaluated frozen metrics.
columns<-data.frame(column=1:6,scenario=rep(c('primary','expanded_1758_12','main_1723_crop_blind_11'),each=2),
 task=rep(c('A','B'),3),group=rep(c('Primary','Expanded cohort','Crop-blind'),each=2),
 label=rep(c('LOSO','Crop-out'),3))
a<-do.call(rbind,lapply(1:6,function(i)do.call(rbind,lapply(models,function(m)do.call(rbind,lapply(baselines,function(b){
 if(i<=2){
  z<-m1[m1$scope=='overall_oof'&m1$task==columns$task[i]&m1$weighting=='study_equal',]
  rm<-z$RMSE[z$estimator==m];rb<-z$RMSE[z$estimator==b]
  delta<-rm-rb;source<-'Phase1 metrics';scenario<-'phase1_main_1723_12_reference'
 }else{
  z<-m2[m2$scenario==columns$scenario[i]&m2$task==columns$task[i]&m2$weighting=='study_equal'&m2$metric=='RMSE',]
  rm<-z$point[z$kind=='absolute'&z$estimator==m];rb<-z$point[z$kind=='absolute'&z$estimator==b]
  delta<-z$point[z$kind=='paired_delta_model_minus_baseline'&z$estimator==m&z$baseline==b&!is.na(z$baseline)]
  source<-'Phase2 metrics';scenario<-columns$scenario[i]
 }
 stopifnot(length(delta)==1,length(rm)==1,length(rb)==1,abs(delta-(rm-rb))<1e-12)
 data.frame(column=i,scenario=scenario,task=columns$task[i],estimator=m,baseline=b,weighting='study_equal',
 model_RMSE=rm,baseline_RMSE=rb,delta_RMSE=delta,half=ifelse(b==baselines[1],'upper','lower'),source=source)
}))))))
stopifnot(nrow(a)==36,all(is.finite(a$delta_RMSE)))
wt(a,'panel_A_frozen_RMSE_contrasts.tsv');wt(columns,'panel_A_column_order.tsv')
delta_limit<-ceiling(max(abs(a$delta_RMSE))/.005)*.005

# B: immutable prediction rows with purely graphical ordering.
dictionary<-rd('dictionary');predictors<-dictionary$column[dictionary$role=='main_predictor']
stopifnot(nrow(ext)==57,!anyDuplicated(ext$comparison_id),length(unique(ext$control_cluster_id))==13,
 length(unique(ext$field_year_id))==5,length(predictors)==12,nrow(unique(ext[predictors]))==1)
known<-colSums(!is.na(ext[predictors]));stopifnot(sum(known==57)==7,sum(known==0)==5)
stopifnot(all(vapply(paste0('prediction_',methods),function(n)length(unique(ext[[n]]))==1,logical(1))))
fields<-c('Vermeersch_2022_Bottelare','Vermeersch_2022_Merelbeke','Vermeersch_2023_Melle',
 'Vermeersch_2023_Merelbeke','Vermeersch_2023_Poperinge')
field_labels<-sub('Vermeersch_','',fields);field_labels<-gsub('_',' ',field_labels)
ext<-ext[order(match(ext$field_year_id,fields),ext$comparison_id),]
ext$display_rank<-seq_len(nrow(ext));ext$field_rank<-match(ext$field_year_id,fields)
ext$barcode_x<-ext$display_rank+(ext$field_rank-1)*.8
bc<-do.call(rbind,lapply(c('Observed',methods),function(m)data.frame(comparison_id=ext$comparison_id,
 field_year_id=ext$field_year_id,control_cluster_id=ext$control_cluster_id,display_rank=ext$display_rank,
 barcode_x=ext$barcode_x,method=m,logRR=if(m=='Observed')ext$external_logRR else ext[[paste0('prediction_',m)]])))
wt(bc,'panel_B_external_barcode.tsv')
wt(ext[,c('comparison_id','field_year_id','control_cluster_id','display_rank','field_rank','barcode_x')],
 'panel_B_external_dependency_and_order.tsv')
wt(data.frame(predictor=predictors,confirmed=known==57,missing=known==0),'external_predictor_profile.tsv')
response_limit<-ceiling(max(bc$logRR)/.05)*.05
stopifnot(min(bc$logRR)>=0)
response_rgb<-grDevices::colorRamp(c('#F2F1ED','#B8A78F','#876A50'),space='Lab')
response_colour<-function(v){z<-response_rgb(v/response_limit);grDevices::rgb(z[,1],z[,2],z[,3],maxColorValue=255)}
wt(data.frame(logRR=seq(0,response_limit,length.out=101),colour=response_colour(seq(0,response_limit,length.out=101))),
 'shared_response_colour_scale.tsv')

# C: all displayed summaries copied directly from frozen external metrics.
fy<-em[em$scope=='field_year'&em$method=='Ridge',]
fy<-fy[match(fields,fy$field_year_id),]
overall<-em[em$scope=='overall'&em$weighting=='field_year_equal',];overall<-overall[match(methods,overall$method),]
comparison<-em[em$scope=='overall'&em$weighting=='comparison'&em$method=='Ridge',]
obs<-data.frame(id=c(fields,'observed_comparison_mean','observed_field_year_equal_mean'),
 label=c(field_labels,'Observed · comparison mean','Observed · field-year equal mean'),
 logRR=c(fy$observed_mean,comparison$observed_mean,overall$observed_mean[1]),kind=c(rep('field-year',5),'comparison mean','field-year equal mean'))
pred<-data.frame(method=methods,label=LABELS[methods],logRR=overall$prediction_mean,field_year_equal_RMSE=overall$RMSE)
ranking<-pred[order(pred$field_year_equal_RMSE,pred$method),];ranking$display_rank<-1:nrow(ranking)
stopifnot(nrow(obs)==7,nrow(pred)==6)
wt(obs,'panel_C_frozen_observed_means.tsv');wt(pred,'panel_C_frozen_prediction_means.tsv');wt(ranking,'panel_C_frozen_RMSE_order.tsv')
wt(fy[,c('field_year_id','n_comparisons','n_control_clusters','observed_mean')],'external_field_year_structure.tsv')

# Shared typographic and graphical helpers; all coordinates are display-only.
txt<-function(label,x,y,size=7,col=INK,face='plain',just='left')grid.text(label,unit(x,'npc'),unit(y,'npc'),just=just,
 gp=gpar(fontfamily=FONT,fontsize=size,col=col,fontface=face))
line<-function(x0,y0,x1,y1,col=RULE,lwd=.6)grid.lines(x=unit(c(x0,x1),'npc'),y=unit(c(y0,y1),'npc'),gp=gpar(col=col,lwd=lwd))
header<-function(tag,title){txt(tag,.02,.964,9.5,face='bold');txt(title,.047,.964,8.5,face='bold')}
capsule_half<-function(x0,x1,y,upper,colour){
 # Half of a pill shape with exact data extent; radius limited for short glyphs.
 xx<-sort(c(x0,x1));left<-convertX(unit(xx[1],'npc'),'mm',valueOnly=TRUE);right<-convertX(unit(xx[2],'npc'),'mm',valueOnly=TRUE)
 yy<-convertY(unit(y,'npc'),'mm',valueOnly=TRUE);r<-min(1.25,(right-left)/2)
 if(r<=0)return(invisible(NULL))
 angle<-seq(pi,0,length.out=40);sign<-if(upper)1 else -1
 lx<-left+r+r*cos(angle[angle>=pi/2]);ly<-yy+sign*r*sin(angle[angle>=pi/2])
 rx<-right-r+r*cos(angle[angle<pi/2]);ry<-yy+sign*r*sin(angle[angle<pi/2])
 grid.polygon(unit(c(left,lx,rx,right),'mm'),unit(c(yy,ly,ry,yy),'mm'),gp=gpar(fill=colour,col=NA))
}
centres<-seq(.19,.91,length.out=6);halfwidth<-.061
draw_A<-function(){
 header('A','Robustness across prespecified sensitivities');txt('Study-equal ΔRMSE',.98,.964,6,col=MID,just='right')
 txt('Upper: vs train observation mean · Lower: vs train study-equal mean',.02,.90,6.2,col=MID)
 for(g in 1:3){ii<-(2*g-1):(2*g);txt(columns$group[ii[1]],mean(centres[ii]),.818,7,face='bold',just='centre');line(centres[ii[1]]-halfwidth,.79,centres[ii[2]]+halfwidth,.79)}
 for(i in 1:6){
  x<-centres[i];txt(columns$label[i],x,.744,6.5,just='centre')
  for(j in 1:3){y<-c(.61,.46,.31)[j];line(x,y-.04,x,y+.04,col='#ADB8BE',lwd=.65)
   for(b in 1:2){d<-a$delta_RMSE[a$column==i&a$estimator==models[j]&a$baseline==baselines[b]]
    capsule_half(x,x+d/delta_limit*halfwidth,y+ifelse(b==1,.002,-.002),b==1,MODEL_COLORS[j])}
  }
  line(x-halfwidth,.18,x+halfwidth,.18,col=RULE,lwd=.55)
  for(v in c(-delta_limit,0,delta_limit)){xx<-x+v/delta_limit*halfwidth;line(xx,.18,xx,.17,col=RULE,lwd=.55);txt(format(v,nsmall=2),xx,.127,5.3,col=MID,just=if(v<0)'left'else if(v>0)'right'else'centre')}
 }
 for(j in 1:3)txt(LABELS[models[j]],.02,c(.61,.46,.31)[j],7.1,col=MODEL_COLORS[j],face='bold')
 txt('Left of 0: model better · Right of 0: baseline better',.02,.041,6.1,col=MID)
}
draw_B<-function(){
 header('B','Independent external response barcode')
 txt('57 comparisons · 13 control clusters · 5 field-years',.98,.914,6.2,col=MID,just='right')
 txt('7 confirmed / 5 missing predictors · 1 unique profile; predictions constant within each model',.02,.866,6.1,col=MID)
 xmap<-function(v).205+(v-.5)/(max(ext$barcode_x)+.5)*.773
 for(i in 1:5){xx<-ext$barcode_x[ext$field_rank==i];lb<-sub(' ','\n',field_labels[i]);txt(lb,xmap(mean(range(xx))),.772,5.8,col=MID,just='centre')}
 display_methods<-c('Observed',methods);ys<-seq(.657,.22,length.out=7)
 for(i in 1:7){m<-display_methods[i];z<-bc[bc$method==m,];cc<-response_colour(z$logRR)
  for(k in 1:nrow(z)){lo<-xmap(z$barcode_x[k]-.5);hi<-xmap(z$barcode_x[k]+.5)
   grid.rect(unit((lo+hi)/2,'npc'),unit(ys[i],'npc'),width=unit(hi-lo+.00005,'npc'),height=unit(2.8,'mm'),gp=gpar(fill=cc[k],col=NA))}
  col<-if(m%in%models)MODEL_COLORS[m]else if(m=='Observed')INK else MID
  txt(if(m=='Observed')'Observed'else LABELS[m],.02,ys[i],6.2,col=col,face=if(m%in%models||m=='Observed')'bold'else'plain')
 }
 # Continuous vector legend; no embedded raster gradient.
 steps<-seq(0,response_limit,length.out=201);cc<-response_colour(head(steps,-1))
 for(i in 1:200)grid.rect(unit(.48+(i-.5)/200*.28,'npc'),unit(.097,'npc'),width=unit(.28/200+.00005,'npc'),height=unit(1.7,'mm'),gp=gpar(fill=cc[i],col=NA))
 txt('logRR',.436,.097,6.2,col=MID,just='right')
 for(v in seq(0,response_limit,length.out=4))txt(sprintf('%.2f',v),.48+v/response_limit*.28,.052,5.7,col=MID,just='centre')
}
ruler_limit<-ceiling(max(obs$logRR)/.05)*.05
ymap<-function(v).125+v/ruler_limit*.70
obs$true_y<-ymap(obs$logRR);obs$label_y<-obs$true_y
# Deterministic, measurement-independent minimum typographic spacing.
ix<-order(obs$true_y);gap<-.054
for(k in 2:length(ix))obs$label_y[ix[k]]<-max(obs$true_y[ix[k]],obs$label_y[ix[k-1]]+gap)
pred$true_y<-ymap(pred$logRR)
pred_order<-order(pred$logRR,pred$method);pred$label_y<-NA_real_;pred$label_y[pred_order]<-seq(.15,.565,length.out=6)
pred$tab_x<-c(.496,.47,.496,.47,.522,.47)
wt(obs,'panel_C_observed_ruler_coordinates.tsv');wt(pred,'panel_C_prediction_ruler_coordinates.tsv')
draw_C<-function(){
 header('C','External calibration and field-year-equal error')
 txt('Observed means',.02,.883,6.5,face='bold');txt('Predicted / training means',.52,.883,6.5,face='bold')
 x<-.445;line(x,ymap(0),x,ymap(ruler_limit),col=INK,lwd=.8)
 for(v in seq(0,ruler_limit,by=.05)){yy<-ymap(v);big<-abs(v/.1-round(v/.1))<1e-9
  line(x-.007,yy,x+.007,yy,col=INK,lwd=if(big).7 else .45)
  if(big)txt(sprintf('%.1f',v),x-.04,yy,5.7,col=MID,just='right')
 }
 txt('logRR',x,.074,6.2,col=MID,just='centre')
 for(i in 1:nrow(obs)){isoverall<-obs$kind[i]!='field-year';col<-if(isoverall)INK else '#7D929E'
  grid.rect(unit(x-.017,'npc'),unit(obs$true_y[i],'npc'),width=unit(.029,'npc'),height=unit(if(isoverall).75 else .55,'mm'),gp=gpar(fill=col,col=NA))
  line(.356,obs$label_y[i],x-.033,obs$true_y[i],col=RULE,lwd=.55)
  txt(obs$label[i],.02,obs$label_y[i],6.3,col=col,face=if(isoverall)'bold'else'plain')
 }
 for(i in 1:nrow(pred)){col<-if(pred$method[i]%in%models)MODEL_COLORS[pred$method[i]]else '#879198'
  line(x+.007,pred$true_y[i],pred$tab_x[i]-.009,pred$true_y[i],col=col,lwd=.35)
  grid.rect(unit(pred$tab_x[i],'npc'),unit(pred$true_y[i],'npc'),width=unit(.018,'npc'),height=unit(.55,'mm'),gp=gpar(fill=col,col=NA))
  line(pred$tab_x[i]+.01,pred$true_y[i],.566,pred$label_y[i],col=col,lwd=.55)
  txt(pred$label[i],.574,pred$label_y[i],6.2,col=col,face=if(pred$method[i]%in%models)'bold'else'plain')
 }
 # Keep ruler numerals readable where a thin observed-label leader passes.
 for(v in seq(0,ruler_limit,by=.1)){yy<-ymap(v)
  grid.rect(unit(x-.051,'npc'),unit(yy,'npc'),width=unit(3.6,'mm'),height=unit(2.0,'mm'),gp=gpar(fill='white',col=NA))
  txt(sprintf('%.1f',v),x-.04,yy,5.7,col=MID,just='right')}
 # Small text list, sorted existing RMSE values; all methods retained.
 txt('Field-year-equal RMSE',.80,.78,6.1,col=MID)
 short<-c(Ridge='Ridge',RandomForest='RF',XGBoost='XGBoost',training_observation_mean='Train obs. mean',training_study_equal_mean='Train study mean',zero_response='Zero')
 for(i in 1:6){m<-ranking$method[i];col<-if(m%in%models)MODEL_COLORS[m]else MID
  yy<-.717-(i-1)*.061;txt(short[m],.80,yy,5.8,col=col);txt(sprintf('%.4f',ranking$field_year_equal_RMSE[i]),.978,yy,5.8,col=col,just='right')}
 txt('Means copied from frozen external metrics; no recalibration or significance test',.02,.034,5.9,col=MID)
}
capture_panel<-function(fun,w,h){tmp<-tempfile(fileext='.pdf');on.exit(unlink(tmp));grid.grabExpr(fun(),width=w,height=h,wrap=TRUE,wrap.grobs=TRUE,
 device=function(width,height)grDevices::cairo_pdf(tmp,width=width,height=height,family=FONT))}
width_mm<-183;height_mm<-225;heights<-c(68,70,87)
grobs<-lapply(1:3,function(i)capture_panel(list(draw_A,draw_B,draw_C)[[i]],width_mm/25.4,heights[i]/25.4))
wrap<-function(g)wrap_elements(full=g,clip=TRUE)
fig<-(wrap(grobs[[1]])/wrap(grobs[[2]])/wrap(grobs[[3]]))+plot_layout(heights=heights)
fig<-fig & theme(plot.background=element_rect(fill='white',colour=NA),plot.margin=margin(0,0,0,0))
pdf_file<-file.path(out,'Figure6_draft.pdf');png_file<-file.path(out,'Figure6_draft.png')
grDevices::cairo_pdf(pdf_file,width=width_mm/25.4,height=height_mm/25.4,family=FONT);print(fig);invisible(dev.off())
ragg::agg_png(png_file,width=width_mm,height=height_mm,units='mm',res=600,background='white');print(fig);invisible(dev.off())
for(i in 0:3){ragg::agg_png(file.path(src,if(i==0)'Figure6_QA_preview.png'else paste0('Figure6_QA_',LETTERS[i],'.png')),
 width=width_mm,height=if(i==0)height_mm else heights[i],units='mm',res=if(i==0)180 else 300,background='white')
 print(if(i==0)fig else wrap(grobs[[i]]));invisible(dev.off())}
after<-vapply(protected,hash,character(1));stopifnot(identical(before,after))
wt(data.frame(file=sub(paste0('^',base,'/'),'',protected),SHA256=before),'protected_source_sha256.tsv')
con<-file(png_file,'rb');invisible(readBin(con,'raw',8));dims<-ppm<-NULL
repeat{len<-readBin(con,'integer',1,size=4,endian='big');if(!length(len))break;typ<-rawToChar(readBin(con,'raw',4));raw<-readBin(con,'raw',len);readBin(con,'raw',4)
 if(typ=='IHDR')dims<-readBin(raw[1:8],'integer',2,size=4,endian='big');if(typ=='pHYs')ppm<-readBin(raw[1:8],'integer',2,size=4,endian='big');if(typ=='IEND')break};close(con)
stopifnot(all(abs(dims-c(width_mm,height_mm)/25.4*600)<2),all(abs(ppm*.0254-600)<.02))
caption<-c('Figure 6 | Robustness and independent external validation.',
 'All plotted results are frozen Phase1, Phase2 and Vermeersch external-validation outputs. No model was fitted, tuned or recalibrated for this figure. No new bootstrap, statistical test, feature analysis or model selection was performed. The figure summarizes prespecified robustness settings and a limited external evidence domain; it does not establish a causal explanation of failure or universal non-transferability.',
 'A, Split capsule field of study-equal RMSE differences, model minus each training-mean baseline. The upper half compares the training observation mean; the lower half compares the training study-equal mean. Both baselines were defined within the corresponding training folds. All six columns and all models use the same signed linear scale, with an individual zero reference. Glyph horizontal extent alone encodes absolute difference; left denotes model improvement and right denotes baseline improvement. No uncertainty interval is encoded. Primary is the frozen 1,723-comparison/169-study, 12-predictor cohort; expanded cohort is the prespecified 1,758-comparison/171-study, 12-predictor sensitivity; crop-blind retains the main cohort with the frozen 11 predictors. LOSO holds out whole studies; crop-out holds out the crop and purges all other rows of its studies. Primary differences are obtained by subtraction of frozen Phase1 metric values; sensitivity differences are copied from frozen Phase2 paired point estimates and verified against their absolute metrics. The displayed metric weighting is study-equal, not observation-equal. All 36 point differences are positive; neither sensitivity restores improvement over the two mean baselines.',
 'B, One-dimensional response barcode. Each narrow rectangle is one of the 57 frozen external comparisons; field-years are ordered 2022 Bottelare, 2022 Merelbeke, 2023 Melle, 2023 Merelbeke and 2023 Poperinge. Comparisons within each field-year follow their frozen IDs, not outcome order. Small blank gaps delimit field-years. All seven strips use exactly the same zero-anchored continuous logRR color scale without clipping; observed external_logRR and all six predictions are copied from the frozen prediction file. Seven of 12 external predictors are confirmed and five remain missing, with exactly one unique predictor profile; each deterministic model consequently has one constant prediction across all comparisons. Cultivar, strain, site and year remain provenance and are not newly introduced model features.',
 'C, Vertical calibration ruler. Left notches show the five within-field-year observed means, the comparison-weighted observed mean and the field-year-equal observed mean. Right tabs show each model prediction mean, the two frozen Schütz training means and zero. Measurement positions use one common linear logRR scale; thin leaders separate crowded labels and do not change values. All means are copied directly from frozen external metrics. The adjacent text list gives field-year-equal RMSE in ascending order for all six methods; ranking is display order only, not a new model-selection step. Field-year-equal RMSE is the square root of the field-year-equal average squared comparison errors, not the mean of the five RMSEs. Ridge slightly improves on the training study-equal mean but does not improve on the training observation mean; RF and XGBoost improve on neither training mean. No significance or equivalence claim is made.',
 'External dependence: 57 comparisons belong to 13 shared-control clusters nested in five field-years from one publication and four locations; Merelbeke repeats across years. They are not 57 independent experiments or five independent publications. Field-year equal weighting limits dominance by larger trials but does not remove dependence. Public reported treatment/control harvest-yield means use common 15% moisture standardization; the dry-matter factor cancels in their frozen log ratio. Predictor completeness and single-profile scope constrain interpretation. Full comparison/cluster memberships, displayed metrics, placement coordinates and color scales are available in source data.')
writeLines(caption,file.path(out,'Figure6_draft_caption.txt'),useBytes=TRUE)
invisible(file.copy(file.path(out,'Figure6_draft_caption.txt'),file.path(src,'Figure6_draft_caption.txt'),overwrite=TRUE))
writeLines(c('# Figure 6 plot data','',
 'Reproduce exports with ../figure6.R using the frozen project files. Only R drawing/export code is executed.',
 'A: 36 frozen study-equal RMSE contrasts, 18 cells, two baseline halves. B: 399 rectangles, 57 comparisons × 7 strips.',
 'C: seven observed means and six model/baseline means copied from existing metrics; exact and label-only ruler positions saved separately.',
 'Caption explains dependence, field-year weighting and the single predictor profile. Source hashes, graphical scales and session info accompany tables.'),file.path(src,'README.md'),useBytes=TRUE)
capture.output(sessionInfo(),file=file.path(src,'R_sessionInfo.txt'))
qc<-list(status='RENDERED_PENDING_VISUAL_QC',sources_unchanged=TRUE,source_sha256=as.list(setNames(before,protected)),
 A_contrasts=nrow(a),A_all_positive=all(a$delta_RMSE>0),A_delta_scale=c(-delta_limit,delta_limit),
 B_comparisons=nrow(ext),B_rectangles=nrow(bc),field_years=5,control_clusters=13,
 confirmed_predictors=7,missing_predictors=5,unique_predictor_profiles=nrow(unique(ext[predictors])),
 C_observed_means=nrow(obs),C_prediction_means=nrow(pred),new_models_fitted=0,new_inferential_analyses=0,
 font=FONT,width_mm=width_mm,height_mm=height_mm,PNG_pixels=as.list(dims),PNG_dpi=as.list(ppm*.0254),
 response_scale=c(0,response_limit),ruler_scale=c(0,ruler_limit),
 output_sha256=list(pdf=hash(pdf_file),png=hash(png_file)),script_sha256=hash(file.path(out,'figure6.R')))
jsonlite::write_json(qc,file.path(src,'figure6_QC.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
message('Figure 6 draft rendered from frozen results: ',paste(dims,collapse=' x '),' px, 600 dpi.')
