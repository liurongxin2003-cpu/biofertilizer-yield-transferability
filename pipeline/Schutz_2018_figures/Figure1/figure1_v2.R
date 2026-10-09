#!/usr/bin/env Rscript
# Visual refinement only: all panel source tables and original draft stay read-only.
invisible(Sys.setlocale('LC_ALL','Chinese (Simplified)_China.utf8'))
req<-c('ggplot2','patchwork','ragg','digest','jsonlite')
missing<-req[!vapply(req,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))stop('Missing R packages: ',paste(missing,collapse=', '))
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(grid)})
args<-commandArgs(trailingOnly=FALSE);sa<-grep('^--file=',args,value=TRUE)
out<-if(length(sa))dirname(normalizePath(sub('^--file=','',sa[1]),winslash='/')) else normalizePath('Schutz_2018_figures/Figure1',winslash='/')
base<-normalizePath(file.path(out,'../..'),winslash='/');src<-file.path(out,'figure1_plot_data')
qa<-file.path(src,'v2');dir.create(qa,recursive=TRUE,showWarnings=FALSE)
hash<-function(p)digest::digest(file=p,algo='sha256',serialize=FALSE)
oldqc<-jsonlite::read_json(file.path(src,'figure1_QC.json'),simplifyVector=TRUE)
scientific_paths<-file.path(base,unlist(oldqc$source_paths))
science_before<-vapply(scientific_paths,hash,character(1))
stopifnot(identical(unname(science_before),unname(unlist(oldqc$source_sha256))))
source_files<-list.files(src,pattern='\\.tsv$',full.names=TRUE)
protected<-c(source_files,file.path(out,c('figure1.R','Figure1_draft.pdf','Figure1_draft.png')))
protected_before<-vapply(protected,hash,character(1))
rd<-function(f)read.delim(file.path(src,f),stringsAsFactors=FALSE,check.names=FALSE)
atlas<-rd('panel_A_study_crop_counts.tsv');crops<-rd('crop_annotations_and_order.tsv')
studies<-rd('study_annotations_and_order.tsv');cf<-rd('panel_B_crop_function_counts.tsv')
ddata<-rd('panel_D_LOSO_study_membership.tsv');edata<-rd('panel_E_crop_out_study_bands.tsv')
detail_data<-rd('panel_E_mixed_crop_study_detail.tsv');fexample<-rd('panel_F_representative_candidate_cluster.tsv')
clusters<-rd('panel_F_all_candidate_control_membership.tsv')
stopifnot(nrow(atlas)==169*40,sum(atlas$comparisons)==1723,nrow(cf)==40*6,sum(cf$comparisons)==1723,
          nrow(studies)==169,length(unique(studies$study_id))==169,nrow(crops)==40,
          nrow(clusters)==1117,length(unique(clusters$candidate_cluster))==315,
          length(unique(clusters$study_id))==90)
studies<-studies[order(studies$rank),];crops<-crops[order(crops$rank),]
study_order<-studies$study_id;crop_order<-crops$crop
d_folds<-unique(ddata$fold_id);d_studies<-unique(ddata$held_out_study);d_ranks<-match(d_studies,study_order)
e_folds<-unique(edata$fold_id);e_crops<-unique(edata$held_out_crop)
stopifnot(identical(as.integer(d_studies),as.integer(unlist(oldqc$representative_LOSO_held_studies))),
          identical(e_crops,as.character(unlist(oldqc$representative_crop_holdouts))),
          unique(detail_data$study_id)==100,sum(detail_data$membership=='Test')==6,
          sum(detail_data$membership=='Purge')==64,
          identical(as.character(unique(fexample$candidate_cluster)),as.character(oldqc$shared_control_example)))
fexample<-fexample[order(fexample$pair_id),]
country_levels<-c('India','Iran','Argentina','Pakistan','Egypt','Other','Multiple')
country_col<-setNames(c('#536F80','#8E939A','#93AFAA','#B59B7A','#91879B','#DFE3E5','#555A60'),country_levels)
stopifnot(all(studies$country_display %in% country_levels))
FONT<-'Arial';INK<-'#242A30';MID<-'#66737D';RULE<-'#CDD4D9';AXIS_PT<-.7
ROLECOL<-c(Train='#D8E2E8',Test='#2C607C',Purge='#BB8751',Absent='#FFFFFF')
count_max<-max(atlas$comparisons,cf$comparisons)
count_palette<-grDevices::colorRampPalette(c('#F1F3F4','#DCE8ED','#94B2C2','#527B94','#183F58'))(256)
count_color<-function(n)count_palette[pmin(256L,pmax(1L,1L+as.integer(round(log1p(n)/log1p(count_max)*255))))]
txt<-function(label,x,y,size=7,col=INK,face='plain',just='left',rot=0,lineheight=1.0){
 grid.text(label,x=unit(x,'npc'),y=unit(y,'npc'),just=just,rot=rot,
           gp=gpar(fontfamily=FONT,fontsize=size,col=col,fontface=face,lineheight=lineheight))
}
rect<-function(x,y,w,h,fill,border=NA,lwd=.2){grid.rect(x=unit(x,'npc'),y=unit(y,'npc'),width=unit(w,'npc'),height=unit(h,'npc'),gp=gpar(fill=fill,col=border,lwd=lwd))}
line<-function(x0,y0,x1,y1,col=RULE,lwd=AXIS_PT,lty=1){grid.segments(x0=unit(x0,'npc'),y0=unit(y0,'npc'),x1=unit(x1,'npc'),y1=unit(y1,'npc'),gp=gpar(col=col,lwd=lwd,lty=lty))}
point<-function(x,y,size=1,col=INK,fill=col){grid.points(x=unit(x,'npc'),y=unit(y,'npc'),pch=21,size=unit(size,'mm'),gp=gpar(col=col,fill=fill,lwd=.4))}
header<-function(tag,title,subtitle){
  yt<-1-convertHeight(unit(2.4,'mm'),'npc',valueOnly=TRUE)
  ys<-yt-convertHeight(unit(4.2,'mm'),'npc',valueOnly=TRUE)
  tx<-.02+convertWidth(unit(4.4,'mm'),'npc',valueOnly=TRUE)
  txt(tag,.02,yt,9.5,face='bold');txt(title,tx,yt,8.5,face='bold');txt(subtitle,.02,ys,6.4,col=MID)
}
capture_panel<-function(fun,w,h){
 tmp<-tempfile(fileext='.pdf');on.exit(unlink(tmp))
 grid.grabExpr(fun(),width=w,height=h,wrap=TRUE,
   device=function(width,height)grDevices::cairo_pdf(tmp,width=width,height=height,family=FONT))
}
draw_country_key<-function(){
 txt('Country',.025,.080,6.2,col=MID)
 for(i in seq_along(country_levels)){
   j<-if(i<=4)i else i-4;x<-.035+(j-1)*.24;y<-if(i<=4).052 else .025
   rect(x,y,.011,.012,country_col[[country_levels[i]]]);txt(country_levels[i],x+.013,y,5.9)
 }
}
draw_count_key<-function(){
 txt('Comparisons (A, B)',.36,.080,6.2)
 lx<-.36;lw<-.61;nc<-150L;vals<-expm1(seq(0,log1p(count_max),length.out=nc))
 rect(lx+(seq_len(nc)-.5)*lw/nc,rep(.052,nc),lw/nc,.012,count_color(vals))
 for(v in c(0,1,5,20,count_max)){
   x<-lx+lw*log1p(v)/log1p(count_max);line(x,.044,x,.038,INK,.45);txt(v,x,.022,5.8,just='center')
 }
}
draw_roles<-function(x,y,include_purge=FALSE){
 labs<-if(include_purge)c('Train','Test','Purge','Absent') else c('Train','Test')
 spacing<-if(include_purge).235 else .21
 for(i in seq_along(labs)){
   xx<-x+(i-1)*spacing
   rect(xx,y,.017,convertHeight(unit(1.5,'mm'),'npc',valueOnly=TRUE),ROLECOL[[labs[i]]],if(labs[i]=='Absent')RULE else NA)
   txt(labs[i],xx+.017,y,6.3)
 }
}

draw_A<-function(){
 header('A','Study × crop evidence atlas','1,723 comparisons · 169 studies · 40 crop groups')
 x0<-.115;x1<-.86;y0<-.25;y1<-.79;dx<-(x1-x0)/40;dy<-(y1-y0)/169
 xx<-x0+(seq_len(40)-.5)*dx;yy<-y1-(seq_len(169)-.5)*dy
 rect((x0+x1)/2,(y0+y1)/2,x1-x0,y1-y0,count_color(0))
 rect(xx[atlas$crop_rank],yy[atlas$study_rank],dx*.97,dy*.94,count_color(atlas$comparisons))
 rect(rep(.098,169),yy,.010,dy,country_col[studies$country_display])
 line(.89,y0,.89,y1);point(.89+.095*studies$comparisons/80,yy,.55,MID)
 txt('Country',.098,.804,5.9,just='center');txt('Study size',.939,.804,6.2,just='center')
 txt('0',.89,y0-.014,5.8,just='center');txt('80',.985,y0-.014,5.8,just='center')
 txt('n',.94,y0-.033,6.1,just='center')
 for(r in c(1,25,50,100,150,169)){txt(r,.078,yy[r],5.8,just='right');line(.082,yy[r],.088,yy[r],INK,.45)}
 txt('Study rank (size ↓)',.029,(y0+y1)/2,6.4,just='center',rot=90)
 # Enlarged track with real support-count axis, not a rank scale.
 sy0<-.827;sy1<-.894;maxsup<-max(crops$study_support)
 txt('Crop support (unique studies)',x0,.916,6.7)
 txt(paste0('max ',maxsup),x1,.916,6.0,just='right',col=MID)
 for(v in c(0,15,30)){
   y<-sy0+(sy1-sy0)*v/maxsup
   line(x0,y,x1,y,if(v==0)RULE else '#E7ECEF',if(v==0).65 else .35)
   txt(v,x0-.018,y,5.9,just='right')
 }
 line(xx,rep(sy0,40),xx,sy0+(sy1-sy0)*crops$study_support/maxsup,MID,.65)
 point(xx,sy0+(sy1-sy0)*crops$study_support/maxsup,.72,INK)
 # Display-only wrapping. Source crop categories and order are unchanged.
 wrapnames<-c('Phaseolus beans'='Phaseolus\nbeans','Finger millet'='Finger\nmillet','Pearl millet'='Pearl\nmillet')
 for(j in seq_len(40)){
   name<-crops$display_label[j];label<-if(name %in% names(wrapnames))wrapnames[[name]] else name
   y<-y0-.010
   if(name=='Phaseolus beans')y<-y0-.045
   if(name=='Finger millet')y<-y0-.060
   txt(label,xx[j],y,5.6,just='right',rot=67,lineheight=.98)
 }
 txt('Crops ordered by study support',.49,.112,6.3,just='center')
 draw_country_key()
}

draw_B<-function(){
 header('B','Crop × function evidence','Comparisons; crop order matches A')
 x0<-.36;x1<-.98;y0<-.25;y1<-.895;dx<-(x1-x0)/6;dy<-(y1-y0)/40
 xx<-x0+(seq_len(6)-.5)*dx;yy<-y1-(seq_len(40)-.5)*dy
 rect(xx[as.integer(cf$function_code)+1],yy[cf$crop_rank],dx*.93,dy*.9,count_color(cf$comparisons))
 contrast<-ifelse(log1p(cf$comparisons)/log1p(count_max)>.65,'#FFFFFF',INK)
 for(i in which(cf$comparisons>0))txt(cf$comparisons[i],xx[as.integer(cf$function_code[i])+1],yy[cf$crop_rank[i]],5.5,col=contrast[i],just='center')
 for(j in seq_len(40))txt(crops$display_label[j],x0-.018,yy[j],5.6,just='right')
 flabels<-c('N-fixing','P-solubilizing','N + P','Unclear','K-solubilizing','AMF')
 for(j in seq_len(6))txt(flabels[j],xx[j],y0-.011,5.9,just='right',rot=60)
 # Preserve the frozen biological scope of source code 4 while following the requested label.
 txt('Original code 4 retains S/K solubilization.',.025,.116,5.6,col=MID)
 draw_count_key()
}

theme_set(theme_classic(base_size=7,base_family=FONT)+theme(
 axis.line=element_line(linewidth=.25,colour=INK),axis.ticks=element_line(linewidth=.25,colour=INK),
 axis.text=element_text(size=6.5,colour=INK),axis.title=element_text(size=6.9),
 plot.margin=margin(3,4,3,3),panel.grid.minor=element_blank()))
profile<-studies;profile$fill_display<-ifelse(profile$rank<=5,'top five','other')
pC<-ggplot(profile)+
 geom_rect(aes(xmin=rank-.5,xmax=rank+.5,ymin=1,ymax=comparisons,fill=fill_display),colour=NA)+
 scale_fill_manual(values=c('top five'='#BCCDD6','other'='#EDF1F3'),guide='none')+
 geom_step(aes(rank,comparisons),direction='hv',colour='#586B78',linewidth=.4)+
 geom_rug(aes(y=comparisons),sides='r',alpha=.18,colour=MID,length=unit(.021,'npc'),linewidth=.3)+
 scale_x_continuous(breaks=c(1,50,100,169),limits=c(.5,174),expand=expansion(mult=c(.015,.01)))+
 scale_y_log10(breaks=c(1,2,5,10,20,40,80),labels=c(1,2,5,10,20,40,80),limits=c(.9,95))+
 labs(x='Study rank by comparison count',y='Comparisons / study (log scale)')+
 theme(panel.grid.major.y=element_line(colour='#EEF0F2',linewidth=.2))
draw_C<-function(){
 header('C','Study evidence concentration',sprintf('Largest study: %.1f%% · top five: %.1f%% of comparisons',100*max(studies$comparison_share),100*sum(studies$comparison_share[1:5])))
 pushViewport(viewport(x=.50,y=.435,width=.99,height=.78))
 grid.draw(ggplotGrob(pC));popViewport()
}

draw_D<-function(){
 header('D','Whole-study LOSO membership','6 representative folds of 169 · columns = studies')
 x0<-.31;x1<-.985;y0<-.25;y1<-.832;dx<-(x1-x0)/169;dy<-(y1-y0)/6
 xx<-x0+(seq_len(169)-.5)*dx;yy<-y1-(seq_len(6)-.5)*dy
 for(i in seq_along(d_folds)){
   dd<-ddata[ddata$fold_id==d_folds[i],];rect(xx[dd$study_rank],rep(yy[i],169),dx*.95,dy*.68,ROLECOL[dd$membership])
   txt(sprintf('S%03d · n=%d',d_studies[i],studies$comparisons[d_ranks[i]]),x0-.022,yy[i],6.8,just='right')
 }
 for(r in c(1,50,100,169)){txt(r,xx[r],y0-.038,6.3,just='center');line(xx[r],y0-.012,xx[r],y0-.024,INK,.45)}
 txt('Study rank (same order as A)',.65,.130,6.8,just='center')
 draw_roles(.045,.064,FALSE)
 txt('Study overlap = 0 (all 169 folds)',.49,.064,6.1)
}

draw_E<-function(){
 header('E','Crop-out with whole-study purge','5 representative folds of 40 · same 169 study columns')
 x0<-.225;x1<-.985;y0<-.445;y1<-.867;dx<-(x1-x0)/169;dy<-(y1-y0)/length(e_crops)
 xx<-x0+(seq_len(169)-.5)*dx
 for(i in seq_along(e_crops)){
   yc<-y1-(i-.5)*dy
   for(b in c('Held-out crop','Other crops')){
     dd<-edata[edata$fold_id==e_folds[i] & edata$band==b,];yb<-yc+if(b=='Held-out crop')dy*.18 else -dy*.18
     rect(xx[dd$study_rank],rep(yb,169),dx*.95,dy*.32,ROLECOL[dd$membership],border=ifelse(dd$membership=='Absent','#EDF0F2',NA),lwd=.12)
     txt(if(b=='Held-out crop')'H' else 'O',x0-.012,yb,5.8,just='right',col=MID)
   }
   txt(crops$display_label[match(e_crops[i],crop_order)],x0-.045,yc,6.8,just='right')
 }
 for(r in c(1,50,100,169))txt(r,xx[r],y0-.018,6.2,just='center')
 txt('H: held-out crop; O: other crops within each study',x0,.389,6.1)
 line(.02,.352,.985,.352,RULE,.5)
 txt('Mixed-crop study S100 (Ortas 2012)',.025,.318,7.3,face='bold')
 txt('Holding out maize excludes all 70 rows from training',.025,.277,6.4)
 bx0<-.025;bw<-.945;by<-.177;bh<-.068;frac<-6/70
 txt('Test: maize (n=6)',bx0,.231,6.3);txt('Purge: other crops (n=64)',.47,.231,6.3)
 rect(bx0+bw*frac/2,by,bw*frac,bh,ROLECOL[['Test']]);rect(bx0+bw*frac+bw*(1-frac)/2,by,bw*(1-frac),bh,ROLECOL[['Purge']])
 txt('6',bx0+bw*frac/2,by,6.3,col='white',just='center');txt('64',bx0+bw*frac+bw*(1-frac)/2,by,6.3,col='white',just='center')
 txt('Other crops from test studies never enter training.',.025,.110,6.3)
 draw_roles(.043,.048,TRUE)
}

draw_F<-function(){
 header('F','Potential shared-control dependence','Repeated control signatures; physical sharing unverified')
 # Plain typographic summary, not cards or an infographic score display.
 xs<-c(.065,.405,.785)
 for(i in 1:3)txt(c('315','1,117','90')[i],xs[i],.822,10,face='bold')
 for(i in 1:3)txt(c('Suspected clusters','Comparisons','Studies')[i],xs[i],.775,6.3)
 txt('Candidate signatures',xs[1],.739,5.7,col=MID)
 txt('64.8% of 1,723',xs[2],.739,5.9,col=MID)
 txt('Of 169',xs[3],.739,5.9,col=MID)
 line(.02,.702,.985,.702,RULE,.5)
 txt('Shared-control example: Mäder 2011',.025,.658,6.9,face='bold')
 txt('Control signature',.06,.606,6.1,col=MID)
 txt('Coded C row',.55,.606,6.1,col=MID,just='center')
 txt('Comparison',.825,.606,6.1,col=MID,just='center')
 ys<-seq(.55,.382,length.out=nrow(fexample));cy<-mean(ys)
 txt('Same\nsignature',.17,cy,6.6,just='center',lineheight=1.15)
 # A bracket groups matching records; dashed ties encode suspected dependence only.
 line(.345,min(ys),.345,max(ys),MID,.65);line(.33,min(ys),.345,min(ys),MID,.65);line(.33,max(ys),.345,max(ys),MID,.65)
 for(i in seq_along(ys)){
   line(.345,ys[i],.465,ys[i],MID,.5,lty=2)
   txt(fexample$control_excel_row[i],.55,ys[i],6.8,just='center')
   txt(paste0('Pair ',fexample$pair_id[i]),.825,ys[i],6.8,just='center')
 }
 txt('Signature match ≠ confirmed physical shared control.',.025,.306,6.1,col=MID)
 line(.02,.252,.985,.252,RULE,.5)
 txt('Minimum safe isolation unit: whole study',.025,.196,7.1,face='bold')
}

# Physical panel heights grow 17–18% for D/E; the top row stays dominant.
gA<-capture_panel(draw_A,122/25.4,125/25.4);gB<-capture_panel(draw_B,61/25.4,125/25.4)
gC<-capture_panel(draw_C,91.5/25.4,54/25.4);gD<-capture_panel(draw_D,91.5/25.4,54/25.4)
gE<-capture_panel(draw_E,91.5/25.4,67/25.4);gF<-capture_panel(draw_F,91.5/25.4,67/25.4)
wrap<-function(g)wrap_elements(full=g,clip=FALSE)
top<-(wrap(gA)|wrap(gB))+plot_layout(widths=c(2,1))
middle<-(wrap(gC)|wrap(gD));bottom<-(wrap(gE)|wrap(gF))
fig<-(top/middle/bottom)+plot_layout(heights=c(125,54,67))
fig<-fig & theme(plot.background=element_rect(fill='white',colour=NA),plot.margin=margin(0,0,0,0))
width_mm<-183;height_mm<-246
pdf_file<-file.path(out,'Figure1_v2.pdf');png_file<-file.path(out,'Figure1_v2.png')
grDevices::cairo_pdf(pdf_file,width=width_mm/25.4,height=height_mm/25.4,family=FONT,onefile=TRUE)
print(fig);invisible(dev.off())
ragg::agg_png(png_file,width=width_mm,height=height_mm,units='mm',res=600,background='white')
print(fig);invisible(dev.off())
ragg::agg_png(file.path(qa,'Figure1_v2_QA_preview.png'),width=width_mm,height=height_mm,units='mm',res=180,background='white')
print(fig);invisible(dev.off())
ragg::agg_png(file.path(qa,'Figure1_v2_QA_A_B.png'),width=width_mm,height=125,units='mm',res=300,background='white')
print(top);invisible(dev.off())
ragg::agg_png(file.path(qa,'Figure1_v2_QA_C_D_E_F.png'),width=width_mm,height=121,units='mm',res=300,background='white')
print((middle/bottom)+plot_layout(heights=c(54,67)));invisible(dev.off())
stopifnot(identical(science_before,vapply(scientific_paths,hash,character(1))),
          identical(protected_before,vapply(protected,hash,character(1))))
con<-file(png_file,'rb');sig<-readBin(con,'raw',8);stopifnot(all(as.integer(sig)==c(137,80,78,71,13,10,26,10)))
dims<-NULL;ppm<-NULL
repeat{
 len<-readBin(con,'integer',1,size=4,endian='big');if(!length(len))break
 type<-rawToChar(readBin(con,'raw',4));raw<-readBin(con,'raw',len);readBin(con,'raw',4)
 if(type=='IHDR')dims<-readBin(raw[1:8],'integer',2,size=4,endian='big')
 if(type=='pHYs')ppm<-readBin(raw[1:8],'integer',2,size=4,endian='big')
 if(type=='IEND')break
};close(con)
stopifnot(all(abs(dims-c(width_mm,height_mm)/25.4*600)<2),all(abs(ppm*.0254-600)<.02))
qc<-list(scientific_content_unchanged=TRUE,six_panels_preserved=TRUE,original_draft_and_plot_tables_unchanged=TRUE,
         source_sha256=as.list(setNames(protected_before,basename(protected))),no_model_training=TRUE,no_performance_read=TRUE,
         panel_heights_mm=list(top=125,D=54,E=67),D_height_increase_percent=(54/46-1)*100,E_height_increase_percent=(67/57-1)*100,
         backend='R + ggplot2/grid + patchwork',font=FONT,width_mm=183,height_mm=246,
         comparison_legend='B; shared with A',country_legend='A',membership_legend='D/E',
         function_display_labels=c('N-fixing','P-solubilizing','N + P','Unclear','K-solubilizing','AMF'),
         function_code_4='S_or_K_solubilizer source scope retained; no reclassification',
         PNG_pixels=as.list(dims),PNG_dpi=as.list(ppm*.0254),output_sha256=list(pdf=hash(pdf_file),png=hash(png_file)),
         visual_review='pending inspection')
jsonlite::write_json(qc,file.path(qa,'Figure1_v2_QC.json'),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(qa,'R_sessionInfo.txt'))
writeLines(c('Figure 1 v2: appearance-only refinement; original scientific content and source data unchanged.',
 'A: study–crop atlas with larger support-count track, axis ticks, display-only wrapping and adjacent country key.',
 'B: exact requested function labels, frozen code/counts retained, original S/K code-4 scope disclosed; shared comparison key adjacent.',
 'C: ranked empirical step profile with light area starting at comparison count=1 on a log scale and rug; no smoothing or fitted trend.',
 'D/E: original representative memberships and S100 maize example; physical row heights increased 17–18%; adjacent Train/Test/Purge keys.',
 'F: same four-pair candidate-control example, bracket/table depiction and unboxed summary; sharing remains suspected.',
 'All other scientific legend statements in the original caption remain valid. No Figure 2 work.'),file.path(qa,'Figure1_v2_refinement_notes.txt'))
message('Figure1_v2 exported: ',paste(dims,collapse=' × '),' pixels, 600 dpi; scientific source tables unchanged.')
