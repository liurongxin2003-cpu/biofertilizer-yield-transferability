#!/usr/bin/env Rscript
# Descriptive outcome landscape only; no model fits, pooled effects or hypothesis tests.
invisible(Sys.setlocale('LC_ALL','Chinese (Simplified)_China.utf8'))
required<-c('ggplot2','patchwork','ragg','digest','jsonlite')
missing<-required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))stop('Missing R packages: ',paste(missing,collapse=', '))
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(grid)})
args<-commandArgs(trailingOnly=FALSE);sa<-grep('^--file=',args,value=TRUE)
out<-if(length(sa))dirname(normalizePath(sub('^--file=','',sa[1]),winslash='/')) else normalizePath('Schutz_2018_figures/Figure2',winslash='/')
base<-normalizePath(file.path(out,'../..'),winslash='/');src<-file.path(out,'figure2_plot_data','final')
dir.create(src,recursive=TRUE,showWarnings=FALSE)
files<-c(matrix='Schutz_2018_phase0/analysis_matrix.tsv',dictionary='Schutz_2018_phase0/data_dictionary.tsv',
         crop_order='Schutz_2018_figures/Figure1/figure1_plot_data/crop_annotations_and_order.tsv',
         Figure1_pdf='Schutz_2018_figures/Figure1/Figure1_v2.pdf',Figure1_png='Schutz_2018_figures/Figure1/Figure1_v2.png',
         draft_pdf='Schutz_2018_figures/Figure2/Figure2_draft.pdf',draft_png='Schutz_2018_figures/Figure2/Figure2_draft.png',
         v2_pdf='Schutz_2018_figures/Figure2/Figure2_v2.pdf',v2_png='Schutz_2018_figures/Figure2/Figure2_v2.png')
hash<-function(p)digest::digest(file=p,algo='sha256',serialize=FALSE)
before<-vapply(file.path(base,files),hash,character(1));names(before)<-names(files)
stopifnot(before[['matrix']]=='66fd6e7ddd65426225167bfa3016f518cad26ff836e0748decdd3afa39209094',
          before[['dictionary']]=='f2e49b119be38ac84214f33cb6bc7abac0808402fb7fa9aa815f9b420fb72cfa')
original<-read.delim(file.path(base,files[['matrix']]),check.names=FALSE,stringsAsFactors=FALSE)
# Use the frozen outcome directly; do not calculate a new effect size from source yields.
dat<-original[,c('observation_id','study_id','pair_id','logRR','crop','inoculant_function')];rm(original)
stopifnot(nrow(dat)==1723,length(unique(dat$study_id))==169,length(unique(dat$crop))==40,
          all(is.finite(dat$logRR)),!anyDuplicated(dat$observation_id),setequal(unique(dat$inoculant_function),0:5))
frozen_values<-dat$logRR
write_tsv<-function(d,name)write.table(d,file.path(src,name),sep='\t',row.names=FALSE,quote=FALSE,na='NA',fileEncoding='UTF-8')
summarize_values<-function(x){
 q<-quantile(x,c(.25,.5,.75),type=7,names=FALSE)
 data.frame(n=length(x),mean=mean(x),q1=q[1],median=q[2],q3=q[3],minimum=min(x),maximum=max(x),
            negative=sum(x<0),zero=sum(x==0),positive=sum(x>0))
}
overall<-summarize_values(dat$logRR)
stopifnot(overall$negative==220,overall$zero==5,overall$positive==1498)
study_ids<-sort(unique(dat$study_id))
study<-do.call(rbind,lapply(study_ids,function(id){
 dd<-dat[dat$study_id==id,];ss<-summarize_values(dd$logRR)
 data.frame(study_id=id,comparisons=nrow(dd),study_mean_logRR=ss$mean,study_median_logRR=ss$median,
            q1=ss$q1,q3=ss$q3,minimum=ss$minimum,maximum=ss$maximum,crop_groups=length(unique(dd$crop)))
}))
stopifnot(nrow(study)==169,sum(study$comparisons)==1723)
study$display_size_mm<-.80+.25*sqrt((study$comparisons-1)/(max(study$comparisons)-1))
study_spread<-summarize_values(study$study_mean_logRR)
crops<-do.call(rbind,lapply(sort(unique(dat$crop)),function(crop){
 dd<-dat[dat$crop==crop,];data.frame(crop=crop,comparisons=nrow(dd),studies=length(unique(dd$study_id)),summarize_values(dd$logRR)[,-1])
}))
crops<-crops[order(-crops$studies,-crops$comparisons,crops$crop),];crops$support_rank<-seq_len(nrow(crops))
crops$display_label<-gsub('_',' ',crops$crop,fixed=TRUE)
crops$display_label<-paste0(toupper(substr(crops$display_label,1,1)),substring(crops$display_label,2))
top12<-crops[1:12,];crop_order<-top12$crop
prior_order<-read.delim(file.path(base,files[['crop_order']]),stringsAsFactors=FALSE)
stopifnot(identical(crop_order,prior_order$crop[1:12]),sum(top12$comparisons)==1301L)
dictionary<-read.delim(file.path(base,files[['dictionary']]),check.names=FALSE,stringsAsFactors=FALSE)
function_def<-dictionary$definition[dictionary$column=='inoculant_function']
code_names<-strsplit(function_def,';',fixed=TRUE)[[1]]
stopifnot(identical(as.integer(sub('=.*','',code_names)),0:5))
function_names<-sub('^[0-5]=','',code_names)
function_labels<-c('N-fixing','P-solubilizing','N + P traits','Unclear','S/K-solubilizing','AMF')
# Display-only names; exact frozen function_names and 0–5 membership are unchanged.
functions<-do.call(rbind,lapply(0:5,function(code){
 dd<-dat[dat$inoculant_function==code,]
 data.frame(function_code=code,function_group=function_names[code+1],display_label=function_labels[code+1],
            comparisons=nrow(dd),studies=length(unique(dd$study_id)),summarize_values(dd$logRR)[,-1])
}))
dat$function_group<-function_names[dat$inoculant_function+1]
write_tsv(dat,'panel_A_all_comparison_logRR.tsv');write_tsv(overall,'panel_A_descriptive_summary.tsv')
write_tsv(study,'panel_B_E_study_response_summary.tsv');write_tsv(study_spread,'panel_B_study_mean_spread.tsv')
write_tsv(crops,'all_crop_support_and_response_summary.tsv');write_tsv(top12,'panel_C_F_top12_crop_summary.tsv')
write_tsv(functions,'panel_D_function_summary.tsv')

# Descriptive Gaussian KDE only; each curve is normalized to its own maximum for display.
# Restrict the visible curve to observed endpoints; raw points retain all extremes.
density_curve<-function(x,weights=NULL){
 den<-density(x,bw=bw.nrd0(x),kernel='gaussian',n=512,from=min(x),to=max(x),weights=weights)
 data.frame(logRR=den$x,density=den$y,relative_density=den$y/max(den$y),bandwidth=den$bw)
}
da<-density_curve(dat$logRR)
study$density_weight<-rep(1/nrow(study),nrow(study))
db<-density_curve(study$study_mean_logRR,weights=study$density_weight)
write_tsv(da,'panel_A_density.tsv');write_tsv(db,'panel_B_study_mean_density.tsv')
half_polygon<-function(curve,base,height,id){
 data.frame(logRR=c(curve$logRR,rev(curve$logRR)),y=c(base+height*curve$relative_density,rep(base,nrow(curve))),group=id)
}
polyA<-half_polygon(da,.235,.54,'comparison');polyB<-half_polygon(db,.36,.45,'study')
set.seed(20261007)
study$display_y<-.20+runif(nrow(study),-.055,.055)
craw<-dat[dat$crop %in% crop_order,];craw$group_row<-13-match(craw$crop,crop_order)
craw$display_y<-craw$group_row-.16+runif(nrow(craw),-.08,.08)
drawdat<-dat;drawdat$group_row<-6-drawdat$inoculant_function
drawdat$display_y<-drawdat$group_row-.16+runif(nrow(drawdat),-.08,.08)
c_density<-list();d_density<-list();c_poly<-list();d_poly<-list()
for(i in seq_along(crop_order)){
 k<-crop_order[i];curve<-density_curve(dat$logRR[dat$crop==k]);curve$crop<-k;curve$group_row<-13-i
 c_density[[i]]<-curve;c_poly[[i]]<-half_polygon(curve,13-i+.035,.38,k)
}
for(code in 0:5){
 curve<-density_curve(dat$logRR[dat$inoculant_function==code]);curve$function_code<-code;curve$group_row<-6-code
 d_density[[code+1]]<-curve;d_poly[[code+1]]<-half_polygon(curve,6-code+.035,.38,as.character(code))
}
cpoly<-do.call(rbind,c_poly);dpoly<-do.call(rbind,d_poly)
write_tsv(study,'panel_B_study_points_display.tsv');write_tsv(craw,'panel_C_raw_points_display.tsv')
write_tsv(drawdat,'panel_D_raw_points_display.tsv');write_tsv(do.call(rbind,c_density),'panel_C_density.tsv')
write_tsv(do.call(rbind,d_density),'panel_D_density.tsv')

# Display bins: not outcome recoding, model labels or validated biological categories.
positive_yield_boundary<-0.10
cutoff<-log1p(positive_yield_boundary)
dat$response_bin<-ifelse(dat$logRR<0,'negative',ifelse(dat$logRR<cutoff,'weak_positive_including_zero','stronger_positive'))
bin_order<-c('negative','weak_positive_including_zero','stronger_positive')
heat<-do.call(rbind,lapply(crop_order,function(k){
 dd<-dat[dat$crop==k,];do.call(rbind,lapply(seq_along(bin_order),function(j){
   count<-sum(dd$response_bin==bin_order[j])
   data.frame(crop=k,support_rank=match(k,crop_order),response_bin=bin_order[j],bin_column=j,
              comparisons=count,crop_comparisons=nrow(dd),proportion=count/nrow(dd),
              exact_zero_comparisons=sum(dd$logRR==0 & dd$response_bin==bin_order[j]))
 }))
}))
stopifnot(nrow(heat)==36,sum(heat$comparisons)==1301,
          all(abs(tapply(heat$proportion,heat$crop,sum)-1)<1e-12),
          all(tapply(heat$comparisons,heat$crop,sum)[crop_order]==top12$comparisons),
          identical(dat$logRR,frozen_values))
write_tsv(heat,'panel_F_crop_response_bins.tsv')
write_tsv(dat[,c('observation_id','study_id','crop','logRR','response_bin')],'panel_F_display_bin_assignments.tsv')
jsonlite::write_json(list(negative='logRR < 0',middle='0 <= logRR < ln(1.10); includes exact zeros',
   stronger='logRR >= ln(1.10)',cutoff_logRR=cutoff,cutoff_yield_increase=positive_yield_boundary,
   role='descriptive display only; unchanged frozen outcome',selection='top12 by unique-study support; ties comparison count then crop label',
   choice='10% descriptive display threshold retained in the user-authorized v2 refinement'),
   file.path(src,'response_bin_definition.json'),pretty=TRUE,auto_unbox=TRUE)

FONT<-'Arial';INK<-'#242A30';MID<-'#66737D';RULE<-'#CDD4D9';BLUE<-'#2C607C';CLOUD<-'#DCE8ED';LIGHT<-'#EEF2F4'
LW_AXIS<-.25;LW_REF<-.30;LW_DENSITY<-.40
XLIM<-c(-1.15,2.60);XBREAKS<-c(-1,0,1,2,2.5)
BLIM<-c(-.30,.80);BBREAKS<-c(-.2,0,.2,.4,.6,.8)
theme_set(theme_classic(base_size=7,base_family=FONT)+theme(
 axis.line=element_line(linewidth=LW_AXIS,colour=INK),axis.ticks=element_line(linewidth=LW_AXIS,colour=INK),
 axis.text=element_text(size=6.3,colour=INK),axis.title=element_text(size=6.8),
 plot.margin=margin(3,4,3,3),panel.grid=element_blank(),legend.title=element_text(size=6),legend.text=element_text(size=5.9)))
zero_ref<-geom_vline(xintercept=0,linetype='dashed',colour='#909BA4',linewidth=LW_REF)
strip_theme<-theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank())
pA<-ggplot()+
 geom_polygon(data=polyA,aes(logRR,y),fill=CLOUD,colour=NA)+
 geom_line(data=data.frame(logRR=da$logRR,y=.235+.54*da$relative_density),aes(logRR,y),colour='#527B94',linewidth=LW_DENSITY)+
 zero_ref+
 geom_segment(data=dat,aes(x=logRR,xend=logRR),y=.12,yend=.17,linewidth=.2,alpha=.13,colour=MID)+
 geom_segment(aes(x=overall$q1,xend=overall$q3),y=.045,yend=.045,colour=BLUE,linewidth=.9)+
 geom_point(aes(x=overall$median,y=.045),size=1.6,colour=BLUE)+
 geom_segment(aes(x=overall$median,xend=overall$median),y=.19,yend=.80,colour=BLUE,linewidth=.3,linetype='dotted')+
 scale_x_continuous(limits=XLIM,breaks=XBREAKS,expand=expansion(mult=0))+
 scale_y_continuous(limits=c(-.01,.84),expand=expansion(mult=0))+
 labs(x='Comparison logRR',y=NULL)+strip_theme
pB<-ggplot()+
 geom_polygon(data=polyB,aes(logRR,y),fill=CLOUD,colour=NA)+
 geom_line(data=data.frame(logRR=db$logRR,y=.36+.45*db$relative_density),aes(logRR,y),colour='#527B94',linewidth=LW_DENSITY)+
 zero_ref+
 geom_point(data=study,aes(study_mean_logRR,display_y,size=display_size_mm),shape=21,fill=BLUE,colour=BLUE,stroke=.12,alpha=.50)+
 scale_size_identity()+
 geom_segment(aes(x=study_spread$q1,xend=study_spread$q3),y=.065,yend=.065,colour=BLUE,linewidth=.9)+
 geom_point(aes(x=study_spread$median,y=.065),size=1.6,colour=BLUE)+
 scale_x_continuous(limits=BLIM,breaks=BBREAKS,expand=expansion(mult=0))+
 scale_y_continuous(limits=c(.01,.85),expand=expansion(mult=0))+
 labs(x='Study mean logRR',y=NULL)+strip_theme
csummary<-top12;csummary$group_row<-13-seq_len(12)
dsummary<-functions;dsummary$group_row<-6-dsummary$function_code
clabel<-sprintf('%s\nn=%d; S=%d',top12$display_label,top12$comparisons,top12$studies)
dlabel<-sprintf('%s\nn=%d; S=%d',functions$display_label,functions$comparisons,functions$studies)
raincloud<-function(poly,raw,summary,labels,groups){
 ggplot()+
  geom_polygon(data=poly,aes(logRR,y,group=group),fill=CLOUD,colour='#91AFC0',linewidth=LW_REF)+
  zero_ref+
  geom_point(data=raw,aes(logRR,display_y),size=.40,alpha=.14,colour='#536F80')+
  geom_segment(data=summary,aes(x=q1,xend=q3,y=group_row-.36,yend=group_row-.36),colour=BLUE,linewidth=.6)+
  geom_point(data=summary,aes(median,group_row-.36),size=.9,colour=BLUE)+
  scale_x_continuous(limits=XLIM,breaks=XBREAKS,expand=expansion(mult=c(.015,.025)))+
  scale_y_continuous(limits=c(.48,groups+.55),breaks=groups:1,labels=labels,expand=expansion(mult=0))+
  labs(x='Comparison logRR',y=NULL)+
  theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_text(size=6.3,lineheight=1.08))
}
pC<-raincloud(cpoly,craw,csummary,clabel,12)
pD_dictionary<-raincloud(dpoly,drawdat,dsummary,
    sprintf('%s\nn=%d; S=%d',functions$function_group,functions$comparisons,functions$studies),6)
pD<-raincloud(dpoly,drawdat,dsummary,dlabel,6)
pE<-ggplot(study,aes(comparisons,study_mean_logRR))+
 geom_hline(yintercept=0,linetype='dashed',colour='#909BA4',linewidth=LW_REF)+
 geom_point(shape=21,fill='#536F80',colour='#536F80',stroke=.12,size=1.3,alpha=.6)+
 scale_x_log10(limits=c(.90,90),breaks=c(1,2,5,10,20,40,80),labels=c('1','2','5','10','20','40','80'),expand=expansion(mult=0))+
 scale_y_continuous(limits=BLIM,breaks=BBREAKS,expand=expansion(mult=0))+
 labs(x='Comparisons per study (log scale)',y='Study mean logRR')
heat$display_row<-13-heat$support_rank;heat$cell_label<-sprintf('%.1f%%',100*heat$proportion)
heat$text_colour<-ifelse(heat$proportion>.65,'white',INK)
pF<-ggplot(heat,aes(bin_column,display_row))+
 geom_tile(aes(fill=proportion),colour='white',linewidth=.35,width=.96,height=.93)+
 geom_text(aes(label=cell_label,colour=text_colour),size=2.25,family=FONT)+scale_colour_identity()+
 scale_fill_gradientn(colours=c('#F1F3F4','#DCE8ED','#94B2C2','#527B94','#183F58'),limits=c(0,1),
    breaks=c(0,.5,1),labels=c('0%','50%','100%'),name='Within-crop comparisons',
    guide=guide_colourbar(title.position='top',barwidth=unit(30,'mm'),barheight=unit(1.5,'mm'),ticks=TRUE))+
 scale_x_continuous(limits=c(.5,3.5),breaks=1:3,
    labels=c('Negative\n<0','Weak\npositive†\n0–ln(1.10)','Stronger\npositive\n≥ln(1.10)'),position='top',expand=expansion(mult=0))+
 scale_y_continuous(limits=c(.5,12.5),breaks=12:1,labels=top12$display_label,expand=expansion(mult=0))+
 labs(x=NULL,y=NULL)+theme(axis.line=element_blank(),axis.ticks=element_blank(),
    axis.text.x=element_text(size=6,lineheight=1.05),axis.text.y=element_text(size=6.3),
    legend.position='bottom',legend.margin=margin(2,0,0,0),legend.box.spacing=unit(1,'mm'))

# Uniform Figure 1 header grammar, with all plots captured/exported by R.
txt<-function(label,x,y,size=7,col=INK,face='plain',just='left')grid.text(label,x=unit(x,'npc'),y=unit(y,'npc'),just=just,gp=gpar(fontfamily=FONT,fontsize=size,col=col,fontface=face))
header<-function(tag,title,subtitle){
 yt<-1-convertHeight(unit(2.4,'mm'),'npc',valueOnly=TRUE);ys<-yt-convertHeight(unit(4.2,'mm'),'npc',valueOnly=TRUE)
 tx<-.02+convertWidth(unit(4.4,'mm'),'npc',valueOnly=TRUE)
 txt(tag,.02,yt,9.5,face='bold');txt(title,tx,yt,8.5,face='bold');txt(subtitle,.02,ys,6.4,col=MID)
}
capture_panel<-function(fun,w,h){
 tmp<-tempfile(fileext='.pdf');on.exit(unlink(tmp))
 grid.grabExpr(fun(),width=w,height=h,wrap=TRUE,
    device=function(width,height)grDevices::cairo_pdf(tmp,width=width,height=height,family=FONT))
}
embed<-function(p,x=.50,y=.45,w=.99,h=.76,freeze_widths_from=NULL){
 pushViewport(viewport(x=x,y=y,width=w,height=h))
 g<-ggplotGrob(p)
 # Build both layouts within the active Arial-capable Cairo device.
 if(!is.null(freeze_widths_from))g$widths<-ggplotGrob(freeze_widths_from)$widths
 grid.draw(g);popViewport()
}
draw_A<-function(){
 header('A','Overall response distribution','All 1,723 comparisons · 169 studies')
 txt(sprintf('Median %.3f · IQR %.3f–%.3f',overall$median,overall$q1,overall$q3),.02,.825,6.3,col=BLUE)
 embed(pA,y=.43,h=.75)
}
draw_B<-function(){
 header('B','Study-level response distribution','169 within-study means · equal-study density')
 txt(sprintf('Median %.3f · IQR %.3f–%.3f',study_spread$median,study_spread$q1,study_spread$q3),.02,.825,6.3,col=BLUE)
 embed(pB,y=.45,h=.72)
 txt('Comparisons/study',.02,.052,6.1,col=MID)
 for(i in seq_along(c(1,20,80))){
  n<-c(1,20,80)[i];x<-.59+(i-1)*.14;size<-.80+.25*sqrt((n-1)/79)
  grid.points(x=unit(x,'npc'),y=unit(.052,'npc'),pch=21,size=unit(size,'mm'),gp=gpar(fill=BLUE,col=BLUE,lwd=.3))
  txt(n,x+.02,.052,5.9)
 }
}
draw_C<-function(){
 header('C','Crop-specific response landscape','Top 12 crops by unique-study support')
 embed(pC,y=.435,h=.79)
 txt('n: comparisons; S: studies · line: IQR; dot: median',.02,.027,5.9,col=MID)
}
draw_D<-function(){
 header('D','Inoculant-function landscape','All six frozen function groups')
 embed(pD,y=.435,h=.79,freeze_widths_from=pD_dictionary)
 txt('Raw points + median/IQR; no inferential intervals',.02,.027,5.8,col=MID)
}
draw_E<-function(){
 header('E','Support versus response','One point per study · no fitted trend')
 embed(pE,y=.42,h=.78)
}
draw_F<-function(){
 header('F','Crop × response summary','Comparison proportions; same crops/order as C')
 embed(pF,y=.425,h=.77)
 txt('† Middle bin includes exact zero; cutoff = 10% yield increase.',.02,.025,5.5,col=MID)
}

height_top<-60;height_middle<-105;height_bottom<-77;width_mm<-183;height_mm<-242
gA<-capture_panel(draw_A,91.5/25.4,height_top/25.4);gB<-capture_panel(draw_B,91.5/25.4,height_top/25.4)
gC<-capture_panel(draw_C,103/25.4,height_middle/25.4);gD<-capture_panel(draw_D,80/25.4,height_middle/25.4)
gE<-capture_panel(draw_E,103/25.4,height_bottom/25.4);gF<-capture_panel(draw_F,80/25.4,height_bottom/25.4)
wrap<-function(g)wrap_elements(full=g,clip=FALSE)
top<-(wrap(gA)|wrap(gB));middle<-(wrap(gC)|wrap(gD))+plot_layout(widths=c(103,80))
bottom<-(wrap(gE)|wrap(gF))+plot_layout(widths=c(103,80))
fig<-(top/middle/bottom)+plot_layout(heights=c(height_top,height_middle,height_bottom))
fig<-fig & theme(plot.background=element_rect(fill='white',colour=NA),plot.margin=margin(0,0,0,0))
pdf_file<-file.path(out,'Figure2_final.pdf');png_file<-file.path(out,'Figure2_final.png')
grDevices::cairo_pdf(pdf_file,width=width_mm/25.4,height=height_mm/25.4,family=FONT,onefile=TRUE)
print(fig);invisible(dev.off())
ragg::agg_png(png_file,width=width_mm,height=height_mm,units='mm',res=600,background='white')
print(fig);invisible(dev.off())
ragg::agg_png(file.path(src,'Figure2_final_QA_preview.png'),width=width_mm,height=height_mm,units='mm',res=180,background='white')
print(fig);invisible(dev.off())
for(item in list(list(name='A_B',plot=top,height=height_top),list(name='C_D',plot=middle,height=height_middle),list(name='E_F',plot=bottom,height=height_bottom))){
 ragg::agg_png(file.path(src,paste0('Figure2_final_QA_',item$name,'.png')),width=width_mm,height=item$height,units='mm',res=300,background='white')
 print(item$plot);invisible(dev.off())
}
after<-vapply(file.path(base,files),hash,character(1));names(after)<-names(files)
stopifnot(identical(before,after),identical(dat$logRR,frozen_values))
con<-file(png_file,'rb');sig<-readBin(con,'raw',8);stopifnot(all(as.integer(sig)==c(137,80,78,71,13,10,26,10)))
dims<-NULL;ppm<-NULL
repeat{
 len<-readBin(con,'integer',1,size=4,endian='big');if(!length(len))break
 typ<-rawToChar(readBin(con,'raw',4));raw<-readBin(con,'raw',len);readBin(con,'raw',4)
 if(typ=='IHDR')dims<-readBin(raw[1:8],'integer',2,size=4,endian='big')
 if(typ=='pHYs')ppm<-readBin(raw[1:8],'integer',2,size=4,endian='big')
 if(typ=='IEND')break
};close(con)
stopifnot(all(abs(dims-c(width_mm,height_mm)/25.4*600)<2),all(abs(ppm*.0254-600)<.02))
qc<-list(cohort=list(comparisons=1723,studies=169),source_sha256=as.list(before),input_and_Figure1_unchanged=TRUE,
   original_logRR_unchanged=TRUE,missing_outcome=0,nonfinite_outcome=0,study_summary='unweighted arithmetic mean per study',
   crop_selection='unique-study support descending; comparison count and crop name tie breaks',top12_comparisons=1301,
   all_1723_comparisons_in_A_D=TRUE,all_169_studies_in_B_E=TRUE,all_1301_top12_comparisons_in_C_F=TRUE,
   no_outlier_removal=TRUE,raw_point_outcome_jitter=FALSE,display_y_jitter_seed=20261007,
   IQR='descriptive sample quartiles type7, not confidence interval',KDE='Gaussian bw=nrd0; shape normalized per group; trim display to observed endpoints',
   top12_crop_groups=as.list(crop_order),function_codes_unchanged=TRUE,function_names_exact_dictionary=as.list(function_names),
   function_display_labels=as.list(function_labels),D_v2_plotting_area_widths_preserved=TRUE,
   B_density_weight_per_study=1/169,B_display_diameter_mm=c(.80,1.05),C_D_raw_point_alpha=.14,
   E_x_axis='log10 transformed count; no regression',F_weighting='comparison proportions; not study-equal',
   no_models_or_performance_files=TRUE,no_meta_analysis_fit=TRUE,no_regression_or_tests=TRUE,
   positive_display_boundary_yield_change=.10,positive_display_boundary_logRR=cutoff,
   middle_bin_includes_zeros=TRUE,top12_exact_zero_comparisons=sum(heat$exact_zero_comparisons),
   backend='R: ggplot2 + grid + patchwork',font=FONT,width_mm=width_mm,height_mm=height_mm,
   PNG_pixels=as.list(dims),PNG_dpi=as.list(ppm*.0254),output_sha256=list(pdf=hash(pdf_file),png=hash(png_file)),visual_review='pending inspection')
jsonlite::write_json(qc,file.path(src,'Figure2_final_QC.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
capture.output(sessionInfo(),file=file.path(src,'R_sessionInfo.txt'))
caption<-c('Figure 2 | Biofertilizer yield-response landscape.',
 'All panels use the frozen primary cohort of 1,723 comparisons from 169 studies. The existing outcome is the natural logarithm of treatment dry yield divided by the original coded control dry yield; coded comparators can contain active inoculants and are not uniformly untreated controls. No yields/effect sizes were recalculated, no observations were removed and no winsorization was applied.',
 'A, All comparison logRR values as a descriptive Gaussian kernel density with raw-value rug. The point and horizontal segment show the comparison median and IQR; the dotted line marks the median and the dashed line zero. Densities are trimmed to observed endpoints and normalized to their maximum for display; the curve is not a fitted treatment-effect distribution.',
 'B, One unweighted arithmetic mean of the recorded comparison logRR values per study. The density gives each of the 169 study means exactly equal weight (1/169), regardless of study comparison count. All 169 points are shown; the reduced display-only point-diameter range of 0.80–1.05 mm indicates comparison counts and does not affect the density, median or IQR. The point/segment below the point cloud are the equally weighted median/IQR across these study means. This is not a pooled/inverse-variance/random-effects estimator. The study-mean axis range reflects its observed span and differs from the individual-comparison range.',
 'C, Raw values, per-crop normalized half-densities and descriptive median/IQR for the 12 crops with the highest number of unique supporting studies. Ties use comparison count and crop label, matching Figure 1. All 1,301 comparisons in these crops are displayed. Labels give comparison counts n and unique-study counts S; the same study can contribute to multiple crops.',
 'D, The same distribution grammar for all six frozen nominal inoculant-function codes. Human-readable display labels N-fixing, P-solubilizing, N + P traits, Unclear, S/K-solubilizing and AMF correspond respectively to frozen N_fixing, P_solubilizing, N_and_P_traits, unclear, S_or_K_solubilizer and AMF (codes 0–5). Only displayed text changes: function codes, membership and scientific scope remain frozen. All 1,723 comparisons and all extremes are retained; C/D point opacity is unchanged from v2.',
 'E, Exact comparison count per study on a log10-transformed count axis versus its unweighted mean logRR. Tick labels show original counts. All 169 studies contribute one point, including the smallest and largest studies; neither counts nor response values are jittered. No regression, smoother, correlation test, uncertainty band or formal support-versus-response inference is fitted.',
 'F, Within-crop percentages of comparison rows, calculated as the number of comparisons in each bin divided by the total comparisons for that crop; these are comparison-based proportions, not study-equal proportions. The unchanged bins are negative (logRR<0), weak-positive-including-zero (0≤logRR<ln(1.10)) and stronger-positive (logRR≥ln(1.10)). The 10% yield-increase cutoff is a descriptive plotting choice, not a validated biological threshold, new outcome or model label. Exact zeros are retained in the middle bin. Percentages are rounded to one decimal and may not sum to exactly 100% visually; exact counts/proportions are supplied as source data.',
 'All distributions and percentages are descriptive summaries of potentially dependent comparison records, not 1,723 independent experiments. IQR denotes sample spread rather than confidence/credible intervals; no p-values, significance stars, pooled diamonds, random-effects meta-analysis or new predictive analysis is performed. Points are jittered only perpendicular to the response axis, with the original logRR unchanged. No Figure 3 is produced.')
writeLines(caption,file.path(src,'Figure2_final_caption.txt'),useBytes=TRUE)
message('Figure 2 final exported: ',paste(dims,collapse=' × '),' pixels, 600 dpi; sources, draft, v2 and logRR unchanged.')
