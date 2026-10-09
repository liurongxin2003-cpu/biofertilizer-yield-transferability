#!/usr/bin/env Rscript
# Appearance-only revision; immutable quantitative tables copied from verified draft.
invisible(Sys.setlocale('LC_ALL','Chinese (Simplified)_China.utf8'))
sa<-grep('^--file=',commandArgs(FALSE),value=TRUE)
out<-dirname(normalizePath(sub('^--file=','',sa[1]),winslash='/'))
base<-normalizePath(file.path(out,'../..'),winslash='/');old<-file.path(out,'figure6_plot_data');src<-file.path(old,'final')
.libPaths(c(file.path(base,'Schutz_2018_figures/R_lib'),.libPaths()))
pk<-c('grid','patchwork','ggplot2','ragg','digest','jsonlite','png')
stopifnot(all(vapply(pk,requireNamespace,logical(1),quietly=TRUE)))
suppressPackageStartupMessages({library(grid);library(patchwork);library(ggplot2)})
dir.create(src,recursive=TRUE,showWarnings=FALSE)
hash<-function(p)digest::digest(file=p,algo='sha256',serialize=FALSE)
read<-function(p)read.delim(p,check.names=FALSE,stringsAsFactors=FALSE,na.strings='NA')
rd<-function(f)read(file.path(old,f))
wt<-function(d,f)write.table(d,file.path(src,f),sep='\t',row.names=FALSE,quote=TRUE,na='NA',fileEncoding='UTF-8')
qc0<-jsonlite::fromJSON(file.path(old,'figure6_QC.json'))
stopifnot(qc0$status=='DRAFT_COMPLETE_QC_PASS',
 hash(file.path(out,'Figure6_draft.pdf'))==qc0$output_sha256$pdf,hash(file.path(out,'Figure6_draft.png'))==qc0$output_sha256$png,
 hash(file.path(out,'figure6_draft.R'))==qc0$script_sha256)
protected<-unique(c(names(qc0$source_sha256),list.files(old,full.names=TRUE,recursive=FALSE)[!dir.exists(list.files(old,full.names=TRUE,recursive=FALSE))],
 file.path(out,'Figure6_v2.pdf'),file.path(out,'Figure6_v2.png'),file.path(out,'figure6_v2.R'),
 file.path(out,c('Figure6_draft.pdf','Figure6_draft.png','figure6_draft.R','Figure6_draft_caption.txt'))))
protected<-unique(c(protected,list.files(file.path(old,'v3'),recursive=TRUE,full.names=TRUE),
 file.path(out,c('Figure6_v3.pdf','Figure6_v3.png','figure6_v3.R'))))
before<-vapply(protected,hash,character(1))
stopifnot(all(vapply(names(qc0$source_sha256),hash,character(1))==unlist(qc0$source_sha256)))
a<-rd('panel_A_frozen_RMSE_contrasts.tsv');columns<-rd('panel_A_column_order.tsv')
bc<-rd('panel_B_external_barcode.tsv');order<-rd('panel_B_external_dependency_and_order.tsv')
obs<-rd('panel_C_frozen_observed_means.tsv');pred<-rd('panel_C_frozen_prediction_means.tsv');ranking<-rd('panel_C_frozen_RMSE_order.tsv')
palette<-rd('model_palette.tsv')
for(f in c('panel_A_frozen_RMSE_contrasts.tsv','panel_A_column_order.tsv','panel_B_external_barcode.tsv',
 'panel_B_external_dependency_and_order.tsv','panel_C_frozen_observed_means.tsv','panel_C_frozen_prediction_means.tsv',
 'panel_C_frozen_RMSE_order.tsv','external_predictor_profile.tsv','external_field_year_structure.tsv','model_palette.tsv')){
 invisible(file.copy(file.path(old,f),file.path(src,f),overwrite=TRUE))
 stopifnot(hash(file.path(old,f))==hash(file.path(src,f)))
}
models<-c('Ridge','RandomForest','XGBoost');baselines<-c('training_observation_mean','training_study_equal_mean')
methods<-c(models,baselines,'zero_response');display_methods<-c('Observed',methods)
LABELS<-setNames(palette$label,palette$method);MODEL_COLORS<-setNames(palette$colour[match(models,palette$method)],models)
FONT<-'Arial';INK<-'#242A30';MID<-'#66737D';RULE<-'#CBD2D7'
LTYPE<-c(Ridge='solid',RandomForest='22',XGBoost='44',training_observation_mean='solid',training_study_equal_mean='22',zero_response='solid')
delta_limit<-qc0$A_delta_scale[2];response_limit<-qc0$response_scale[2];ruler_limit<-qc0$ruler_scale[2]
stopifnot(nrow(a)==36,nrow(bc)==399,nrow(obs)==7,nrow(pred)==6)
response_rgb<-grDevices::colorRamp(c('#F4F2EE','#C9B58D','#9E794E','#543B2C'),space='Lab')
response_colour<-function(v){z<-response_rgb(v/response_limit);grDevices::rgb(z[,1],z[,2],z[,3],maxColorValue=255)}
wt(data.frame(logRR=seq(0,response_limit,length.out=101),colour=response_colour(seq(0,response_limit,length.out=101))),
 'shared_response_colour_scale.tsv')
txt<-function(label,x,y,size=7,col=INK,face='plain',just='left')grid.text(label,unit(x,'npc'),unit(y,'npc'),just=just,
 gp=gpar(fontfamily=FONT,fontsize=size,col=col,fontface=face))
line<-function(x0,y0,x1,y1,col=RULE,lwd=.6,lty='solid')grid.lines(x=unit(c(x0,x1),'npc'),y=unit(c(y0,y1),'npc'),gp=gpar(col=col,lwd=lwd,lty=lty,lineend='butt'))
header<-function(tag,title){txt(tag,.02,.964,9.5,face='bold');txt(title,.047,.964,8.5,face='bold')}
centres<-seq(.19,.91,length.out=6);halfwidth<-.061
ap<-a;ap$x_zero<-centres[ap$column];ap$x_marker<-ap$x_zero+ap$delta_RMSE/delta_limit*halfwidth
ap$y<-c(.61,.46,.31)[match(ap$estimator,models)]+ifelse(ap$half=='upper',.016,-.016)
ap$shape<-ifelse(ap$half=='upper','solid rectangular bar','lighter rectangular bar');wt(ap,'panel_A_marker_coordinates.tsv')
draw_A<-function(){
 header('A','Robustness across prespecified sensitivities');txt('Study-equal ΔRMSE',.98,.964,6,col=MID,just='right')
 txt('Upper: vs train observation mean · Lower: vs train study-equal mean',.02,.90,6.2,col=MID)
 for(g in 1:3){ii<-(2*g-1):(2*g);txt(columns$group[ii[1]],mean(centres[ii]),.818,7,face='bold',just='centre');line(centres[ii[1]]-halfwidth,.79,centres[ii[2]]+halfwidth,.79)}
 for(i in 1:6){x<-centres[i];txt(columns$label[i],x,.744,6.5,just='centre')
  for(j in 1:3){y<-c(.61,.46,.31)[j]
   line(x-halfwidth,y,x+halfwidth,y,col='#D8DEE1',lwd=.45);line(x,y-.042,x,y+.042,col='#85949D',lwd=.65)
   for(b in 1:2){z<-ap[ap$column==i&ap$estimator==models[j]&ap$baseline==baselines[b],]
    col<-MODEL_COLORS[j]
    fill<-if(b==1)col else grDevices::adjustcolor(col,alpha.f=.48)
    grid.rect(unit((x+z$x_marker)/2,'npc'),unit(z$y,'npc'),width=unit(abs(z$x_marker-x),'npc'),height=unit(1.1,'mm'),gp=gpar(fill=fill,col=col,lwd=.35))
   }
  }
  line(x-halfwidth,.18,x+halfwidth,.18,col=RULE,lwd=.55)
  for(v in c(-delta_limit,0,delta_limit)){xx<-x+v/delta_limit*halfwidth;line(xx,.18,xx,.17,col=RULE,lwd=.55)
   txt(format(v,nsmall=2),xx,.127,5.3,col=MID,just=if(v<0)'left'else if(v>0)'right'else'centre')}
 }
 for(j in 1:3)txt(LABELS[models[j]],.02,c(.61,.46,.31)[j],7.1,col=MODEL_COLORS[j],face='bold')
 txt('Left of 0: model better · Right of 0: baseline better',.02,.041,6.1,col=MID)
}
fields<-unique(order$field_year_id);field_labels<-gsub('_',' ',sub('Vermeersch_','',fields))
bp<-bc;bp$field_rank<-match(bp$field_year_id,fields);bp$display_x<-bp$display_rank+(bp$field_rank-1)*1.8
xmap<-function(v).025+(v-.5)/(max(bp$display_x)+.5)*.657
bp$x_left<-xmap(bp$display_x-.5);bp$x_right<-xmap(bp$display_x+.5)
bp$height_mm<-ifelse(bp$method=='Observed',3.8,2.2);bp$fill<-response_colour(bp$logRR)
wt(bp,'panel_B_barcode_display_coordinates.tsv')
draw_B<-function(){
 header('B','Independent external response barcode')
 txt('57 comparisons · 13 control clusters · 5 field-years',.98,.914,6.2,col=MID,just='right')
 txt('7 confirmed / 5 missing predictors · 1 unique profile; predictions constant within each model',.02,.866,6.1,col=MID)
 for(i in 1:5){xx<-bp$display_x[bp$field_rank==i&bp$method=='Observed'];txt(sub(' ','\n',field_labels[i]),xmap(mean(range(xx))),.772,5.3,col=MID,just='centre')}
 ys<-seq(.655,.22,length.out=7)
 for(i in 1:4){z<-bp[bp$field_rank==i&bp$method=='Observed',];z2<-bp[bp$field_rank==i+1&bp$method=='Observed',]
  xx<-(max(z$x_right)+min(z2$x_left))/2;line(xx,.177,xx,.707,col='#E3E6E8',lwd=.45)}
 txt('logRR',.978,.762,5.8,col=MID,just='right')
 for(i in 1:7){m<-display_methods[i];z<-bp[bp$method==m,]
  for(k in 1:nrow(z))grid.rect(unit((z$x_left[k]+z$x_right[k])/2,'npc'),unit(ys[i],'npc'),width=unit(z$x_right[k]-z$x_left[k]+.00004,'npc'),height=unit(z$height_mm[k],'mm'),gp=gpar(fill=z$fill[k],col=NA))
  if(m!='Observed')for(f in 1:5){zz<-z[z$field_rank==f,];lo<-min(zz$x_left);hi<-max(zz$x_right)
   grid.rect(unit((lo+hi)/2,'npc'),unit(ys[i],'npc'),width=unit(hi-lo,'npc'),height=unit(2.2,'mm'),gp=gpar(fill=NA,col='#B9C0C4',lwd=.35))}
  col<-if(m%in%models)MODEL_COLORS[m]else if(m=='Observed')INK else MID
  txt(if(m=='Observed')'Observed'else LABELS[m],.710,ys[i],6,col=col,face=if(m=='Observed'||m%in%models)'bold'else'plain')
  if(m!='Observed')txt(sprintf('%.4f',pred$logRR[match(m,pred$method)]),.978,ys[i],5.5,col=col,just='right')
 }
 steps<-seq(0,response_limit,length.out=201);cc<-response_colour(head(steps,-1))
 for(i in 1:200)grid.rect(unit(.25+(i-.5)/200*.28,'npc'),unit(.095,'npc'),width=unit(.28/200+.00004,'npc'),height=unit(1.7,'mm'),gp=gpar(fill=cc[i],col=NA))
 txt('logRR',.214,.095,6.2,col=MID,just='right')
 for(v in seq(0,response_limit,length.out=4))txt(sprintf('%.2f',v),.25+v/response_limit*.28,.051,5.7,col=MID,just='centre')
}
# Orthogonal label ladders leave the exact continuous ruler positions unchanged.
ymap<-function(v).125+v/ruler_limit*.70
obs$true_y<-ymap(obs$logRR);obs$label_y<-obs$true_y;ix<-order(obs$true_y)
for(k in 2:length(ix))obs$label_y[ix[k]]<-max(obs$true_y[ix[k]],obs$label_y[ix[k-1]]+.056)
obs$elbow_x<-NA_real_;obs$elbow_x[ix]<-seq(.342,.389,length.out=nrow(obs))
pred$true_y<-ymap(pred$logRR);ix<-order(pred$logRR,pred$method)
pred$label_y<-pmax(pred$true_y,.15)
# The same orthogonal ladder with shorter bridges and adequate label spacing.
for(k in 2:length(ix))pred$label_y[ix[k]]<-max(pred$label_y[ix[k]],pred$label_y[ix[k-1]]+.055)
pred$elbow_x<-NA_real_;pred$elbow_x[ix]<-seq(.507,.45,length.out=nrow(pred))
# Semantic C colors: location hue and year shade, fixed models, grey baselines.
C_LOCATION_HUES<-c(Bottelare='#5F3178',Merelbeke='#A23F65',Melle='#426D98',Poperinge='#9C782A')
tint<-function(col,f=.28){v<-grDevices::col2rgb(col);v<-round(v*(1-f)+255*f);grDevices::rgb(v[1,],v[2,],v[3,],maxColorValue=255)}
field<-obs$kind=='field-year'
location<-ifelse(field,sub('Vermeersch_[0-9]{4}_','',obs$id),NA_character_)
year<-rep(NA_integer_,nrow(obs));year[field]<-as.integer(substr(obs$label[field],1,4))
C_OBS_COLORS<-setNames(rep(NA_character_,nrow(obs)),obs$id)
C_OBS_COLORS[field]<-C_LOCATION_HUES[location[field]]
C_OBS_COLORS[field&year==2022]<-tint(C_OBS_COLORS[field&year==2022])
C_OBS_COLORS[obs$kind=='comparison mean']<-'#1A1A1A'
C_OBS_COLORS[obs$kind=='field-year equal mean']<-'#555555'
C_PRED_COLORS<-c(MODEL_COLORS,training_observation_mean='#777777',training_study_equal_mean='#999999',zero_response='#ADADAD')
stopifnot(!anyNA(C_OBS_COLORS),identical(unname(C_PRED_COLORS[models]),unname(MODEL_COLORS)))
wt(data.frame(id=c(obs$id,pred$method),domain=c(rep('Observed',nrow(obs)),rep('Prediction / baseline',nrow(pred))),
 colour=c(unname(C_OBS_COLORS[obs$id]),unname(C_PRED_COLORS[pred$method])),
 location=c(location,rep(NA_character_,nrow(pred))),year=c(year,rep(NA_integer_,nrow(pred))),
 semantic_role=c(ifelse(field,'Location hue / year shade',ifelse(obs$kind=='comparison mean','Observed aggregate: near-black','Observed aggregate: deep-grey')),
 ifelse(pred$method%in%models,'Fixed manuscript model color','Neutral training baseline'))),'panel_C_semantic_colour_palette.tsv')
wt(obs,'panel_C_observed_ladder_coordinates.tsv');wt(pred,'panel_C_prediction_ladder_coordinates.tsv')
draw_C<-function(){
 header('C','External calibration ladder and field-year-equal error')
 txt('Observed means',.02,.89,6.5,face='bold');txt('Predicted / training means',.54,.89,6.5,face='bold')
 x<-.43;line(x,ymap(0),x,ymap(ruler_limit),col=INK,lwd=.8)
 for(v in seq(0,ruler_limit,by=.05)){yy<-ymap(v);big<-abs(v/.1-round(v/.1))<1e-9
  line(x-.005,yy,x+.005,yy,col=INK,lwd=if(big).7 else .45)}
 txt('logRR',x,.074,6.2,col=MID,just='centre')
 obs_labels<-obs$label;obs_labels[obs$kind=='comparison mean']<-'Obs. comparison mean';obs_labels[obs$kind=='field-year equal mean']<-'Obs. field-year-equal mean'
 for(i in 1:nrow(obs)){overall<-obs$kind[i]!='field-year';col<-C_OBS_COLORS[obs$id[i]];lw<-if(overall)1.6 else .8
  # Left measurement ticks, followed by orthogonal bridges to label rungs.
  line(x-.031,obs$true_y[i],x-.006,obs$true_y[i],col=col,lwd=lw)
  line(obs$elbow_x[i],obs$true_y[i],x-.033,obs$true_y[i],col=col,lwd=.6)
  line(obs$elbow_x[i],obs$label_y[i],obs$elbow_x[i],obs$true_y[i],col=col,lwd=.6)
  line(.322,obs$label_y[i],obs$elbow_x[i],obs$label_y[i],col=col,lwd=.6)
  line(.306,obs$label_y[i],.321,obs$label_y[i],col=col,lwd=lw)
  txt(obs_labels[i],.02,obs$label_y[i],6.0,col=col,face=if(overall)'bold'else'plain')
  txt(sprintf('%.4f',obs$logRR[i]),.289,obs$label_y[i],5.7,col=col,just='right')
 }
 for(i in 1:nrow(pred)){m<-pred$method[i];col<-C_PRED_COLORS[m];lt<-LTYPE[m]
  line(x+.006,pred$true_y[i],x+.018,pred$true_y[i],col=col,lwd=1.0,lty=lt)
  line(x+.018,pred$true_y[i],pred$elbow_x[i],pred$true_y[i],col=col,lwd=.65,lty=lt)
  line(pred$elbow_x[i],pred$true_y[i],pred$elbow_x[i],pred$label_y[i],col=col,lwd=.65,lty=lt)
  line(pred$elbow_x[i],pred$label_y[i],.524,pred$label_y[i],col=col,lwd=.65,lty=lt)
  line(.524,pred$label_y[i],.543,pred$label_y[i],col=col,lwd=if(m=='zero_response').8 else 1.4,lty=lt)
  txt(pred$label[i],.55,pred$label_y[i],5.9,col=col,face=if(m%in%models)'bold'else'plain')
  txt(sprintf('%.4f',pred$logRR[i]),.768,pred$label_y[i],5.6,col=col,just='right')
 }
 # Axis labels on a small white backing remain separate from orthogonal rungs.
 for(v in seq(0,ruler_limit,by=.1)){yy<-ymap(v)
  grid.rect(unit(.397,'npc'),unit(yy,'npc'),width=unit(3.4,'mm'),height=unit(2.0,'mm'),gp=gpar(fill='white',col=NA))
  txt(sprintf('%.1f',v),.405,yy,5.7,col=MID,just='right')}
 txt('Field-year-equal',.804,.83,6.1,col=MID);txt('RMSE',.804,.792,6.1,col=MID)
 short<-c(Ridge='Ridge',RandomForest='RF',XGBoost='XGBoost',training_observation_mean='Train obs. mean',training_study_equal_mean='Train study mean',zero_response='Zero')
 for(i in 1:6){m<-ranking$method[i];col<-C_PRED_COLORS[m];yy<-.72-(i-1)*.060
  txt(short[m],.804,yy,5.6,col=col);txt(sprintf('%.4f',ranking$field_year_equal_RMSE[i]),.980,yy,5.6,col=col,just='right')}
 txt('Observed colours denote field location; shade denotes year.',.02,.034,5.9,col=MID)
}
capture_panel<-function(fun,w,h){tmp<-tempfile(fileext='.pdf');on.exit(unlink(tmp));grid.grabExpr(fun(),width=w,height=h,wrap=TRUE,wrap.grobs=TRUE,
 device=function(width,height)grDevices::cairo_pdf(tmp,width=width,height=height,family=FONT))}
width_mm<-183;height_mm<-225;heights<-c(63,70,92)
grobs<-lapply(1:3,function(i)capture_panel(list(draw_A,draw_B,draw_C)[[i]],width_mm/25.4,heights[i]/25.4))
wrap<-function(g)wrap_elements(full=g,clip=TRUE)
fig<-(wrap(grobs[[1]])/wrap(grobs[[2]])/wrap(grobs[[3]]))+plot_layout(heights=heights)
fig<-fig & theme(plot.background=element_rect(fill='white',colour=NA),plot.margin=margin(0,0,0,0))
pdf_file<-file.path(out,'Figure6_final.pdf');png_file<-file.path(out,'Figure6_final.png')
grDevices::cairo_pdf(pdf_file,width=width_mm/25.4,height=height_mm/25.4,family=FONT);print(fig);invisible(dev.off())
ragg::agg_png(png_file,width=width_mm,height=height_mm,units='mm',res=600,background='white');print(fig);invisible(dev.off())
for(i in 0:3){ragg::agg_png(file.path(src,if(i==0)'Figure6_QA_preview.png'else paste0('Figure6_QA_',LETTERS[i],'.png')),
 width=width_mm,height=if(i==0)height_mm else heights[i],units='mm',res=if(i==0)180 else 300,background='white')
 print(if(i==0)fig else wrap(grobs[[i]]));invisible(dev.off())}
# R-only grayscale conversion of the R-rendered QA image, never used as data.
rgb<-png::readPNG(file.path(src,'Figure6_QA_preview.png'))
grey<-.2126*rgb[,,1]+.7152*rgb[,,2]+.0722*rgb[,,3]
png::writePNG(grey,file.path(src,'Figure6_QA_grayscale.png'))
after<-vapply(protected,hash,character(1));stopifnot(identical(before,after))
wt(data.frame(file=protected,SHA256=before),'protected_source_sha256.tsv')
cap<-readLines(file.path(out,'Figure6_draft_caption.txt'),encoding='UTF-8')
cap[3]<-sub('A, Split capsule field','A, Paired horizontal bars',cap[3],fixed=TRUE)
cap[3]<-sub('The upper half compares','The upper solid bar compares',cap[3],fixed=TRUE)
cap[3]<-sub('the lower half compares','the lower lighter bar compares',cap[3],fixed=TRUE)
cap[3]<-sub('Glyph horizontal extent alone encodes absolute difference','Bar signed endpoint and horizontal extent from zero encode the difference',cap[3],fixed=TRUE)
cap[4]<-paste(cap[4],'The observed strip is 3.8 mm tall (+35.7% relative to the draft); other strips are 2.2 mm. An enhanced, identical continuous color mapping is used for every row, with no row-specific normalization. Thin outlines delimit prediction groups and faint vertical lines mark larger field-year gaps. Direct right-hand labels give the frozen constant predictions.')
cap[5]<-sub('C, Vertical calibration ruler.','C, External calibration ladder.',cap[5],fixed=TRUE)
cap[5]<-sub('Right tabs show','Right label rungs show',cap[5],fixed=TRUE)
cap[5]<-sub('Measurement positions use one common linear logRR scale; thin leaders separate crowded labels and do not change values.',
 'Measurement ticks on the central ruler retain one common linear logRR scale. Direct name/value label rungs are spaced for legibility and joined to the exact measurement ticks by orthogonal horizontal/vertical steps only; label-rung height is typographic placement, not an additional measurement position. No diagonal connector is used. The two overall observed means use thicker reference ticks. Observed colours denote field location; shade denotes year. Bottelare, Merelbeke, Melle and Poperinge have distinct location hues; 2022 observations use a lighter tint, and Merelbeke 2023 uses the deeper tone of the same hue as Merelbeke 2022. Year changes shade only, without new symbols or line patterns. Observed comparison and field-year-equal aggregate means are near-black and deep-grey. Ridge, RF and XGBoost use the fixed manuscript blue, green and brown, respectively; training baselines use grey shades. RMSE entries follow the same model/baseline color semantics. Prediction labels are placed closer to the exact ruler values with fixed minimum spacing to shorten orthogonal bridges; ruler measurements, observed-label positions and all numeric values remain unchanged. A and B are identical to v3.',cap[5],fixed=TRUE)
writeLines(cap,file.path(out,'Figure6_final_caption.txt'),useBytes=TRUE);writeLines(cap,file.path(src,'Figure6_final_caption.txt'),useBytes=TRUE)
capture.output(sessionInfo(),file=file.path(src,'R_sessionInfo.txt'))
con<-file(png_file,'rb');invisible(readBin(con,'raw',8));dims<-ppm<-NULL
repeat{len<-readBin(con,'integer',1,size=4,endian='big');if(!length(len))break;typ<-rawToChar(readBin(con,'raw',4));raw<-readBin(con,'raw',len);readBin(con,'raw',4)
 if(typ=='IHDR')dims<-readBin(raw[1:8],'integer',2,size=4,endian='big');if(typ=='pHYs')ppm<-readBin(raw[1:8],'integer',2,size=4,endian='big');if(typ=='IEND')break};close(con)
stopifnot(identical(dims,c(4322L,5314L)),all(abs(ppm*.0254-600)<.02))
qc<-list(status='FINAL_RENDERED_PENDING_QC',scientific_tables_byte_identical=TRUE,protected_files_unchanged=TRUE,
 A_contrasts=nrow(a),B_values=nrow(bc),C_means=nrow(obs)+nrow(pred),C_RMSE=nrow(ranking),
 no_new_fitting_or_analysis=TRUE,A_B_model_palette_unchanged=TRUE,C_semantic_colour_mapping=TRUE,A_B_unchanged=TRUE,A_delta_scale=qc0$A_delta_scale,
 response_scale=qc0$response_scale,ruler_scale=qc0$ruler_scale,observed_strip_height_mm=3.8,prediction_strip_height_mm=2.2,
 width_mm=width_mm,height_mm=height_mm,PNG_pixels=as.list(dims),PNG_dpi=as.list(ppm*.0254),
 output_sha256=list(pdf=hash(pdf_file),png=hash(png_file)),script_sha256=hash(file.path(out,'figure6_final.R')))
jsonlite::write_json(qc,file.path(src,'figure6_QC.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
message('Figure 6 final rendered. Scientific values, original source data, draft and Figures1–5 unchanged.')

