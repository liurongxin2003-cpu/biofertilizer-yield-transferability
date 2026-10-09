#!/usr/bin/env Rscript
# Appearance-only final refinement of Figure3 v4; existing quantitative inputs are read-only.
invisible(Sys.setlocale('LC_ALL','Chinese (Simplified)_China.utf8'))
sa<-grep('^--file=',commandArgs(trailingOnly=FALSE),value=TRUE)
out<-if(length(sa))dirname(normalizePath(sub('^--file=','',sa[1]),winslash='/'))else normalizePath('Schutz_2018_figures/Figure3',winslash='/')
base<-normalizePath(file.path(out,'../..'),winslash='/');src<-file.path(out,'figure3_plot_data','v4');qa<-file.path(out,'figure3_plot_data','final')
.libPaths(c(file.path(base,'Schutz_2018_figures','R_lib'),.libPaths()))
required<-c('ggplot2','ggalluvial','patchwork','ragg','digest','jsonlite','png')
missing<-required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))stop('Missing R packages: ',paste(missing,collapse=', '))
suppressPackageStartupMessages({library(ggplot2);library(ggalluvial);library(patchwork);library(grid)})
dir.create(qa,recursive=TRUE,showWarnings=FALSE)
hash<-function(p)digest::digest(file=p,algo='sha256',serialize=FALSE)
read<-function(f)read.delim(file.path(src,f),check.names=FALSE,stringsAsFactors=FALSE)
v4qc<-jsonlite::read_json(file.path(src,'Figure3_v4_QC.json'),simplifyVector=TRUE)
template<-file.path(out,'figure3_v4.R')
# CODE-ONLY DERIVATION: validate the edited template against its public identity
# and the recorded original identity, without changing or spoofing the frozen QC.
template_actual_sha<-hash(template)
if(template_actual_sha!=v4qc$reproducible_script_sha256){
  identity<-read.delim(file.path(dirname(base),'provenance','source_code_map.tsv'),check.names=FALSE,stringsAsFactors=FALSE)
  identity<-identity[identity$file=='pipeline/Schutz_2018_figures/Figure3/figure3_v4.R',]
  stopifnot(nrow(identity)==1,identity$original_sha256==v4qc$reproducible_script_sha256,
    identity$final_sha256==template_actual_sha)
}
protected<-c(list.files(src,full.names=TRUE,recursive=TRUE),template,file.path(out,c('Figure3_v4.pdf','Figure3_v4.png')),
 file.path(base,c('Schutz_2018_phase1/phase1_predictions.tsv','Schutz_2018_phase1/phase1_metrics.tsv',
 'Schutz_2018_phase1/phase1_bootstrap.tsv','Schutz_2018_phase1/run/bootstrap_draws.npz','Schutz_2018_phase1/evidence/protocol.json',
 'Schutz_2018_figures/Figure1/Figure1_v2.pdf','Schutz_2018_figures/Figure1/Figure1_v2.png',
 'Schutz_2018_figures/Figure2/Figure2_final.pdf','Schutz_2018_figures/Figure2/Figure2_final.png')))
before<-vapply(protected,hash,character(1))
stopifnot(hash(file.path(out,'Figure3_v4.pdf'))==v4qc$output_sha256$pdf,hash(file.path(out,'Figure3_v4.png'))==v4qc$output_sha256$png)
models<-c('Ridge','RandomForest','XGBoost');LABELS<-c(Ridge='Ridge',RandomForest='RF',XGBoost='XGB')
MODEL_COLORS<-c(Ridge='#2C607C',RandomForest='#548578',XGBoost='#967B64')
states<-c('Negative','Weak positive','Stronger positive');STATE_COLORS<-setNames(c('#748E9D','#BEC4C8','#B39C7D'),states)
FONT<-'Arial';INK<-'#242A30';MID<-'#66737D';RULE<-'#CBD2D7';OBS<-'#AAB3BA'
means<-read('panel_A_ranked_response_envelopes.tsv')
intervals<-read('panel_C_nested_interval_blocks.tsv');cs<-read('panel_C_frozen_point_and_CI.tsv')
blodes<-read('panel_B_alluvial_lodes.tsv');blodes$axis<-factor(blodes$axis,levels=c('Observed','Predicted'))
blodes$state<-factor(blodes$state,levels=states);blodes$estimator<-factor(blodes$estimator,levels=models)
d<-read('panel_D_signed_radial_study_bars.tsv');ring<-read('panel_D_observed_state_inner_ring.tsv')
scales<-read('panel_D_radial_scale.tsv');wheel_base<-setNames(scales$zero_radius,scales$estimator)
display<-jsonlite::read_json(file.path(src,'display_definition.json'),simplifyVector=TRUE)
angular_total<-display$angular_total;radial_max<-6.3
zero_circles<-do.call(rbind,lapply(models,function(est)data.frame(estimator=est,x=seq(10.5,179.5,length.out=600),radius=unname(wheel_base[est]))))
track_labels<-data.frame(estimator=models,x=0,radius=wheel_base[models],label=LABELS[models],colour=MODEL_COLORS[models])
rank_labels<-data.frame(x=c(11,179),radius=6.05,label=c('1','169'))
stopifnot(nrow(means)==507,nrow(d)==507,nrow(ring)==169,nrow(intervals)==18,nrow(cs)==6)

# Reuse only the plotting/composition definitions, excluding all v4 data creation and export code.
l<-readLines(template,encoding='UTF-8')
first<-which(grepl('^theme_set\\(',l));last<-which(grepl('^pdf_file<-',l))-1L
stopifnot(length(first)==1,length(last)==1,first<last)
plot_code<-paste(l[first:last],collapse='\n')
replace_exact<-function(old,new){
 stopifnot(length(gregexpr(old,plot_code,fixed=TRUE)[[1]])==1,grepl(old,plot_code,fixed=TRUE))
 plot_code<<-sub(old,new,plot_code,fixed=TRUE)
}
replace_exact("x='Study rank (observed mean ascending)'","x='Studies ranked by observed mean logRR'")
replace_exact("label='169\\nstudies',size=2.1","label='169\\nstudies',size=1.85")
replace_exact('xx<-c(.075,.28,.475)[i]','xx<-c(.075,.225,.375)[i]')
replace_exact("width=unit(4,'mm'),height=unit(2.3,'mm')","width=unit(3.4,'mm'),height=unit(2.0,'mm')")
replace_exact("x=unit(c(.66,.66),'npc')","x=unit(c(.545,.545),'npc')")
replace_exact("txt('Point estimate',.685,.90","txt('Point estimate',.570,.90")
replace_exact('embed(pD,y=.50,w=1.04,h=.80)','embed(pD,y=.50,w=1.1648,h=.896)')
replace_exact("txt('Inward: model better',.07,.115","txt('Inward: model better',.07,.170")
replace_exact("txt('Outward: baseline better',.07,.080","txt('Outward: baseline better',.07,.135")
eval(parse(text=plot_code))
pdf_file<-file.path(out,'Figure3_final.pdf');png_file<-file.path(out,'Figure3_final.png')
grDevices::cairo_pdf(pdf_file,width=width_mm/25.4,height=height_mm/25.4,family=FONT);print(fig);invisible(dev.off())
ragg::agg_png(png_file,width=width_mm,height=height_mm,units='mm',res=600,background='white');print(fig);invisible(dev.off())
ragg::agg_png(file.path(qa,'Figure3_final_QA_preview.png'),width=width_mm,height=height_mm,units='mm',res=180,background='white');print(fig);invisible(dev.off())
for(it in list(list(name='A',plot=wrap(gA),height=hA),list(name='B',plot=wrap(gB),height=hB),list(name='C_D',plot=bottom,height=hBottom))){
 ragg::agg_png(file.path(qa,paste0('Figure3_final_QA_',it$name,'.png')),width=width_mm,height=it$height,units='mm',res=300,background='white')
 print(it$plot);invisible(dev.off())
}
stopifnot(identical(png::readPNG(file.path(src,'Figure3_v4_QA_B.png')),png::readPNG(file.path(qa,'Figure3_final_QA_B.png'))))
after<-vapply(protected,hash,character(1));stopifnot(identical(before,after))
caption<-readLines(file.path(src,'Figure3_v4_caption.txt'),encoding='UTF-8')
k<-grep('^A,',caption);stopifnot(length(k)==1)
caption[k]<-sub('Grey is','Study rank is used only for ordering and does not represent chronological or temporal order. Grey is',caption[k],fixed=TRUE)
writeLines(caption,file.path(out,'Figure3_final_caption.txt'),useBytes=TRUE)
write.table(data.frame(file=protected,sha256=before),file.path(qa,'protected_input_sha256.tsv'),sep='\t',quote=FALSE,row.names=FALSE,fileEncoding='UTF-8')
con<-file(png_file,'rb');invisible(readBin(con,'raw',8));dims<-ppm<-NULL
repeat{
 len<-readBin(con,'integer',1,size=4,endian='big');if(!length(len))break
 typ<-rawToChar(readBin(con,'raw',4));raw<-readBin(con,'raw',len);readBin(con,'raw',4)
 if(typ=='IHDR')dims<-readBin(raw[1:8],'integer',2,size=4,endian='big')
 if(typ=='pHYs')ppm<-readBin(raw[1:8],'integer',2,size=4,endian='big')
 if(typ=='IEND')break
};close(con)
stopifnot(identical(dims,c(4322L,5149L)),all(abs(ppm*.0254-600)<.02))
qc<-list(status='DATA_AND_EXPORT_PASS_VISUAL_REVIEW_PENDING',quantitative_data_unchanged=TRUE,all_protected_hashes_unchanged=TRUE,
 B_preview_pixel_identical_to_v4=TRUE,four_panels_unchanged=TRUE,new_fits=0,new_resampling=0,new_quantiles=0,
 A_axis='Studies ranked by observed mean logRR',A_rank_not_chronology_in_caption=TRUE,
 C_legend_spacing_compacted=TRUE,C_zero_line_and_interval_blocks_unchanged=TRUE,
 D_wheel_linear_enlargement=1.12,D_center_text_reduced=TRUE,D_existing_directions_moved_closer=TRUE,
 D_radial_encoding_and_all_169_studies_unchanged=TRUE,width_mm=width_mm,height_mm=height_mm,font=FONT,
 PNG_pixels=as.list(dims),PNG_dpi=as.list(ppm*.0254),pdf_sha256=hash(pdf_file),png_sha256=hash(png_file),
 caption_sha256=hash(file.path(out,'Figure3_final_caption.txt')),template_sha256=hash(template))
jsonlite::write_json(qc,file.path(qa,'Figure3_final_QC.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
capture.output(sessionInfo(),file=file.path(qa,'R_sessionInfo.txt'))
message('Figure 3 final rendered; all protected inputs unchanged; B pixel-identical; 600 dpi.')
