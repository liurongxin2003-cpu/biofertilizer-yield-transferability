#!/usr/bin/env Rscript
# Figure 1: structural counts and frozen validation membership only.
# No models, effects, performance files, cohort modifications or outcome reanalysis.
invisible(Sys.setlocale('LC_ALL', 'Chinese (Simplified)_China.utf8'))
required <- c('ggplot2','patchwork','ragg','digest','jsonlite','png')
missing_packages <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing_packages)) stop('Missing R packages: ',paste(missing_packages,collapse=', '))
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(grid)})
args <- commandArgs(trailingOnly=FALSE)
script_arg <- grep('^--file=',args,value=TRUE)
out <- if (length(script_arg)) dirname(normalizePath(sub('^--file=','',script_arg[1]),winslash='/')) else normalizePath('Schutz_2018_figures/Figure1',winslash='/')
base <- normalizePath(file.path(out,'../..'),winslash='/')
data_dir <- file.path(out,'figure1_plot_data');dir.create(data_dir,recursive=TRUE,showWarnings=FALSE)
source_paths <- c(matrix='Schutz_2018_phase0/analysis_matrix.tsv',
                  splits='Schutz_2018_phase0/split_manifest.tsv',
                  cohort_metadata='Schutz_2018_phase_minus1/phase_minus1_cohort_manifest.tsv',
                  study_metadata='Schutz_2018_phase_minus1/phase_minus1_study_audit.tsv',
                  dependencies='Schutz_2018_phase_minus1/phase_minus1_dependency_audit.tsv')
hash <- function(p) digest::digest(file=p,algo='sha256',serialize=FALSE)
before <- vapply(file.path(base,source_paths),hash,character(1))
names(before) <- names(source_paths)
stopifnot(before[['matrix']]=='66fd6e7ddd65426225167bfa3016f518cad26ff836e0748decdd3afa39209094',
          before[['splits']]=='0c9c0881db2fc520af24aec352f2f9066801cd7cdd0e57ea61c3270b4e0c6583')
rd <- function(name) read.delim(file.path(base,source_paths[[name]]),check.names=FALSE,stringsAsFactors=FALSE,na.strings='NA')
matrix_source <- rd('matrix')
dat <- matrix_source[,c('observation_id','study_id','pair_id','crop','inoculant_function')]
rm(matrix_source)
splits <- rd('splits'); meta <- rd('cohort_metadata'); study_meta <- rd('study_metadata');dep <- rd('dependencies')
stopifnot(nrow(dat)==1723L,length(unique(dat$study_id))==169L,length(unique(dat$crop))==40L,
          length(unique(dat$inoculant_function))==6L,!anyDuplicated(dat$observation_id))
dat$key <- paste(dat$study_id,dat$pair_id,sep=':');meta$key <- paste(meta$study_id,meta$pair_id,sep=':')
stopifnot(!anyDuplicated(meta$key),setequal(dat$key,meta$key))
mm <- match(dat$key,meta$key)
stopifnot(all(dat$crop==meta$crop_group[mm]))
dat$country <- meta$country[mm]
dat$study_label <- study_meta$study_label[match(dat$study_id,study_meta$study_id)]
dat$function_code <- as.character(dat$inoculant_function)
function_labels <- c('N fixing','P solubilizing','N + P traits','Unclear','S/K solubilizing','AMF')
function_names <- c('N_fixing','P_solubilizing','N_and_P_traits','unclear','S_or_K_solubilizer','AMF')
dat$function_group <- function_names[dat$inoculant_function+1L]
stopifnot(all(dat$function_group==meta$function_group[mm]))

write_tsv <- function(d,name) write.table(d,file.path(data_dir,name),sep='\t',row.names=FALSE,quote=FALSE,na='NA',fileEncoding='UTF-8')
unpack <- function(x) if (is.na(x)||!nzchar(x)) character() else strsplit(as.character(x),',',fixed=TRUE)[[1]]
member_digest <- function(ids) digest::digest(paste0(paste(ids,collapse='\n'),'\n'),algo='sha256',serialize=FALSE)
roles <- c('Train','Test','Purge'); membership <- vector('list',nrow(splits)); fold_qc <- vector('list',nrow(splits))
for(i in seq_len(nrow(splits))) {
  s <- splits[i,];ids <- lapply(c('train_ids','test_ids','purge_ids'),function(k)unpack(s[[k]]));names(ids)<-roles
  stopifnot(length(unique(unlist(ids,use.names=FALSE)))==1723L,setequal(unlist(ids,use.names=FALSE),dat$observation_id))
  for(j in seq_along(roles)) {
    low <- tolower(roles[j]);stopifnot(length(ids[[j]])==s[[paste0(low,'_n')]],member_digest(ids[[j]])==s[[paste0(low,'_ids_sha256')]])
  }
  role <- rep(NA_character_,1723L)
  for(j in seq_along(roles)) role[match(ids[[j]],dat$observation_id)] <- roles[j]
  trst <- unique(dat$study_id[role=='Train']);testst <- unique(dat$study_id[role=='Test']);pust <- unique(dat$study_id[role=='Purge'])
  stopifnot(!anyNA(role),!length(intersect(trst,testst)),!length(intersect(trst,pust)))
  if(s$design=='LOSO') stopifnot(all(role[dat$study_id==as.integer(s$held_out_level)]=='Test'),!any(role=='Purge'))
  if(s$design=='LOCO_whole_study_purge') {
    heldstudies<-unique(dat$study_id[dat$crop==s$held_out_level])
    stopifnot(all(role[dat$crop==s$held_out_level]=='Test'),
              all(role[dat$study_id %in% heldstudies & dat$crop!=s$held_out_level]=='Purge'),
              all(role[!dat$study_id %in% heldstudies]=='Train'))
  }
  membership[[i]] <- role
  fold_qc[[i]] <- data.frame(fold_id=s$fold_id,design=s$design,held_out_level=s$held_out_level,
                           train_n=sum(role=='Train'),test_n=sum(role=='Test'),purge_n=sum(role=='Purge'),
                           train_test_study_overlap=0,train_purge_study_overlap=0,complete_partition=TRUE,member_hashes_match=TRUE)
}
names(membership)<-splits$fold_id
write_tsv(do.call(rbind,fold_qc),'all_folds_membership_qc.tsv')
stopifnot(sum(splits$design=='LOSO')==169L,sum(splits$design=='LOCO_whole_study_purge')==40L)

studies <- aggregate(list(comparisons=rep(1L,nrow(dat))),list(study_id=dat$study_id),sum)
studies <- studies[order(-studies$comparisons,studies$study_id),]
studies$rank <- seq_len(nrow(studies));studies$comparison_share<-studies$comparisons/1723
studies$cumulative_comparison_share<-cumsum(studies$comparisons)/1723
studies$study_label <- study_meta$study_label[match(studies$study_id,study_meta$study_id)]
studies$countries <- vapply(studies$study_id,function(id)paste(sort(unique(dat$country[dat$study_id==id])),collapse=' | '),character(1))
studies$crop_groups <- vapply(studies$study_id,function(id)length(unique(dat$crop[dat$study_id==id])),integer(1))
country_support <- aggregate(list(comparisons=rep(1L,nrow(dat))),list(country=dat$country),sum)
country_support$studies <- vapply(country_support$country,function(co)length(unique(dat$study_id[dat$country==co])),integer(1))
country_support <- country_support[order(-country_support$studies,-country_support$comparisons,country_support$country),]
# Ties in country support use country name, not any outcome or model metric.
top_countries <- country_support$country[order(-country_support$studies,country_support$country)][1:5]
studies$country_display <- ifelse(grepl(' | ',studies$countries,fixed=TRUE),'Multiple',ifelse(studies$countries %in% top_countries,studies$countries,'Other'))
country_levels <- c(top_countries,'Other','Multiple')
country_col <- setNames(c('#536F80','#8E939A','#93AFAA','#B59B7A','#91879B','#DFE3E5','#555A60'),country_levels)
write_tsv(country_support,'country_support.tsv');write_tsv(studies,'study_annotations_and_order.tsv')

crops <- aggregate(list(comparisons=rep(1L,nrow(dat))),list(crop=dat$crop),sum)
crops$study_support <- vapply(crops$crop,function(co)length(unique(dat$study_id[dat$crop==co])),integer(1))
crops<-crops[order(-crops$study_support,-crops$comparisons,crops$crop),];crops$rank<-seq_len(nrow(crops))
crops$display_label <- gsub('_',' ',crops$crop,fixed=TRUE)
crops$display_label <- paste0(toupper(substr(crops$display_label,1,1)),substring(crops$display_label,2))
study_order<-studies$study_id;crop_order<-crops$crop
atlas<-as.data.frame(as.table(table(factor(dat$study_id,levels=study_order),factor(dat$crop,levels=crop_order))))
names(atlas)<-c('study_id','crop','comparisons');atlas$study_id<-as.integer(as.character(atlas$study_id));atlas$crop<-as.character(atlas$crop)
atlas$study_rank<-match(atlas$study_id,study_order);atlas$crop_rank<-match(atlas$crop,crop_order)
atlas$country_display<-studies$country_display[atlas$study_rank]
atlas$study_size<-studies$comparisons[atlas$study_rank];atlas$crop_study_support<-crops$study_support[atlas$crop_rank]
cf<-as.data.frame(as.table(table(factor(dat$crop,levels=crop_order),factor(dat$function_code,levels=as.character(0:5)))))
names(cf)<-c('crop','function_code','comparisons');cf$crop<-as.character(cf$crop);cf$function_code<-as.character(cf$function_code)
cf$crop_rank<-match(cf$crop,crop_order);cf$function_label<-function_labels[as.integer(cf$function_code)+1]
stopifnot(nrow(atlas)==169L*40L,nrow(cf)==40L*6L,sum(atlas$comparisons)==1723L,sum(cf$comparisons)==1723L)
write_tsv(atlas,'panel_A_study_crop_counts.tsv');write_tsv(crops,'crop_annotations_and_order.tsv')
write_tsv(cf,'panel_B_crop_function_counts.tsv');write_tsv(studies,'panel_C_study_concentration.tsv')

# Representative folds selected exclusively by structural support/rank.
d_ranks<-c(1L,2L,3L,85L,168L,169L);d_studies<-study_order[d_ranks]
d_folds<-vapply(d_studies,function(id)splits$fold_id[splits$design=='LOSO' & splits$held_out_level==as.character(id)],character(1))
ddata<-do.call(rbind,lapply(seq_along(d_folds),function(i){
  role<-membership[[d_folds[i]]]
  outd<-data.frame(fold_id=d_folds[i],held_out_study=d_studies[i],held_out_study_rank=d_ranks[i],study_id=study_order,study_rank=seq_along(study_order))
  outd$membership<-vapply(study_order,function(id){x<-unique(role[dat$study_id==id]);stopifnot(length(x)==1);x},character(1))
  outd$comparisons<-studies$comparisons;outd
}))
write_tsv(ddata,'panel_D_LOSO_study_membership.tsv')
loco<-splits[splits$design=='LOCO_whole_study_purge',]
less_supported<-loco[!loco$held_out_level %in% crop_order[1:4],]
less_supported<-less_supported[order(-less_supported$purge_n,less_supported$test_n,less_supported$held_out_level),]
e_crops<-c(crop_order[1:4],less_supported$held_out_level[1])
e_folds<-vapply(e_crops,function(co)loco$fold_id[loco$held_out_level==co],character(1))
edata<-do.call(rbind,lapply(seq_along(e_folds),function(i){
  role<-membership[[e_folds[i]]]
  do.call(rbind,lapply(study_order,function(id){
    idx<-which(dat$study_id==id);held<-dat$crop[idx]==e_crops[i]
    do.call(rbind,lapply(c(TRUE,FALSE),function(h){
      ids<-idx[held==h];mem<-if(length(ids))unique(role[ids]) else 'Absent'
      stopifnot(length(mem)==1L)
      data.frame(fold_id=e_folds[i],held_out_crop=e_crops[i],study_id=id,study_rank=match(id,study_order),
                 band=if(h)'Held-out crop' else 'Other crops',membership=mem,comparisons=length(ids))
    }))
  }))
}))
write_tsv(edata,'panel_E_crop_out_study_bands.tsv')
eobs<-do.call(rbind,lapply(seq_along(e_folds),function(i)data.frame(fold_id=unname(e_folds[i]),held_out_crop=unname(e_crops[i]),
   observation_id=dat$observation_id,study_id=dat$study_id,crop=dat$crop,membership=membership[[e_folds[i]]])))
write_tsv(eobs,'panel_E_crop_out_observation_membership.tsv')
detail_fold<-e_folds[which.max(loco$purge_n[match(e_folds,loco$fold_id)])]
detail_role<-membership[[detail_fold]];candidate_ids<-unique(dat$study_id[detail_role=='Purge'])
candidate_counts<-vapply(candidate_ids,function(id)sum(detail_role=='Purge' & dat$study_id==id),integer(1))
detail_study<-candidate_ids[order(-candidate_counts,candidate_ids)][1]
detail_data<-data.frame(dat[dat$study_id==detail_study,c('observation_id','study_id','crop')],membership=detail_role[dat$study_id==detail_study])
stopifnot(!any(detail_data$membership=='Train'),all(c('Test','Purge') %in% detail_data$membership))
write_tsv(detail_data,'panel_E_mixed_crop_study_detail.tsv')

# Existing candidate control signatures: they are NOT proven unique physical controls.
cluster_members<-do.call(rbind,lapply(seq_len(nrow(dep)),function(i){
  pairs<-as.integer(strsplit(dep$pair_ids[i],',',fixed=TRUE)[[1]])
  ii<-which(dat$study_id==dep$study_id[i] & dat$pair_id %in% pairs)
  stopifnot(length(ii)==dep$observations[i])
  data.frame(candidate_cluster=dep$candidate_cluster[i],study_id=dep$study_id[i],pair_id=dat$pair_id[ii],
             observation_id=dat$observation_id[ii],crop=dat$crop[ii],
             status='suspected_shared_control_not_proven',minimum_isolation_unit='study')
}))
stopifnot(nrow(dep)==315L,nrow(cluster_members)==1117L,length(unique(dep$study_id))==90L,!anyDuplicated(cluster_members$observation_id))
small_examples<-dep[dep$observations>=3 & dep$observations<=4,]
small_examples<-small_examples[order(small_examples$study_id,as.integer(sub('.*,','',small_examples$pair_ids))),]
fcluster<-small_examples$candidate_cluster[1]
fexample<-cluster_members[cluster_members$candidate_cluster==fcluster,]
fexample<-fexample[order(fexample$pair_id),]
source_minimal<-rd('matrix')[,c('observation_id','treatment_excel_row','control_excel_row')]
fexample<-merge(fexample,source_minimal,by='observation_id',sort=FALSE);fexample<-fexample[order(fexample$pair_id),]
write_tsv(cluster_members,'panel_F_all_candidate_control_membership.tsv');write_tsv(fexample,'panel_F_representative_candidate_cluster.tsv')
stats<-data.frame(quantity=c('comparisons','studies','crop_groups','function_groups','countries','largest_study_n','largest_study_share',
  'top5_study_share','single_study_crop_groups','suspected_control_clusters','comparisons_in_suspected_control_clusters','studies_with_suspected_controls'),
  value=c(1723,169,40,6,length(unique(dat$country)),max(studies$comparisons),max(studies$comparisons)/1723,
           sum(studies$comparisons[1:5])/1723,sum(crops$study_support==1),315,1117,90))
write_tsv(stats,'figure1_structural_summary.tsv')

FONT<-'Arial';INK<-'#242A30';MID<-'#66737D';RULE<-'#CDD4D9'
ROLECOL<-c(Train='#D8E2E8',Test='#2C607C',Purge='#BB8751',Absent='#FFFFFF')
count_max<-max(atlas$comparisons,cf$comparisons)
count_palette<-grDevices::colorRampPalette(c('#F1F3F4','#DCE8ED','#94B2C2','#527B94','#183F58'))(256)
count_color<-function(n)count_palette[pmin(256L,pmax(1L,1L+as.integer(round(log1p(n)/log1p(count_max)*255))))]
txt<-function(label,x,y,size=7,col=INK,face='plain',just='left',rot=0){
 grid.text(label,x=unit(x,'npc'),y=unit(y,'npc'),just=just,rot=rot,gp=gpar(fontfamily=FONT,fontsize=size,col=col,fontface=face))
}
rect<-function(x,y,w,h,fill,border=NA,lwd=.2){grid.rect(x=unit(x,'npc'),y=unit(y,'npc'),width=unit(w,'npc'),height=unit(h,'npc'),gp=gpar(fill=fill,col=border,lwd=lwd))}
line<-function(x0,y0,x1,y1,col=RULE,lwd=.4,lty=1){grid.segments(x0=unit(x0,'npc'),y0=unit(y0,'npc'),x1=unit(x1,'npc'),y1=unit(y1,'npc'),gp=gpar(col=col,lwd=lwd,lty=lty))}
point<-function(x,y,size=1.0,col=INK,fill=col){grid.points(x=unit(x,'npc'),y=unit(y,'npc'),pch=21,size=unit(size,'mm'),gp=gpar(col=col,fill=fill,lwd=.4))}
header<-function(tag,title,subtitle,compact=FALSE){txt(paste(tag,title),.02,.987,8.5,face='bold');txt(subtitle,.02,if(compact).918 else .952,6.4,col=MID)}
capture_panel<-function(fun,w,h){
  tmp<-tempfile(fileext='.pdf');on.exit(unlink(tmp))
  grid.grabExpr(fun(),width=w,height=h,wrap=TRUE,
                device=function(width,height)grDevices::cairo_pdf(tmp,width=width,height=height,family=FONT))
}

draw_A<-function(){
 header('A','Study × crop evidence atlas','1,723 comparisons · 169 studies · 40 crop groups')
 x0<-.12;x1<-.80;y0<-.19;y1<-.845;dx<-(x1-x0)/40;dy<-(y1-y0)/169
 xx<-x0+(seq_len(40)-.5)*dx;yy<-y1-(seq_len(169)-.5)*dy
 rect((x0+x1)/2,(y0+y1)/2,x1-x0,y1-y0,count_color(0))
 rect(xx[atlas$crop_rank],yy[atlas$study_rank],dx*.97,dy*.94,count_color(atlas$comparisons))
 rect(rep(.102,169),yy,.010,dy,country_col[studies$country_display])
 # Equal-height study-size points, deliberately not a bar annotation.
 line(.84,y0,.84,y1);point(.84+.13*studies$comparisons/max(studies$comparisons),yy,.55,MID)
 txt('Country',.102,.866,5.6,just='center');txt('Study size',.902,.866,6,just='center')
 txt('0',.84,y0-.012,5.5,just='center');txt('80',.97,y0-.012,5.5,just='center')
 txt('n comparisons',.905,y0-.031,5.7,just='center')
 for(r in c(1,25,50,100,150,169)){
   txt(r,.082,yy[r],5.4,just='right');line(.086,yy[r],.091,yy[r],INK,.35)
 }
 txt('Study rank (size ↓)',.031,(y0+y1)/2,6.2,just='center',rot=90)
 # Crop-support annotation, aligned to all 40 matrix columns.
 sy0<-.873;sy1<-.916
 line(x0,sy0,x1,sy0)
 line(xx,rep(sy0,40),xx,sy0+(sy1-sy0)*crops$study_support/max(crops$study_support),MID,.5)
 point(xx,sy0+(sy1-sy0)*crops$study_support/max(crops$study_support),.6,INK)
 txt('Crop support: unique studies',x0,.931,6.1)
 txt(paste0('max ',max(crops$study_support)),.80,.931,5.6,just='right',col=MID)
 for(j in seq_len(40))txt(crops$display_label[j],xx[j],y0-.009,5.2,just='right',rot=60)
 txt('Crops ordered by study support',.48,.037,6.2,just='center')
}

draw_B<-function(){
 header('B','Crop × function evidence','Comparisons; crop order matches A')
 x0<-.36;x1<-.98;y0<-.19;y1<-.914;dx<-(x1-x0)/6;dy<-(y1-y0)/40
 xx<-x0+(seq_len(6)-.5)*dx;yy<-y1-(seq_len(40)-.5)*dy
 rect(xx[as.integer(cf$function_code)+1],yy[cf$crop_rank],dx*.93,dy*.9,count_color(cf$comparisons))
 positive<-cf$comparisons>0
 contrast<-ifelse(log1p(cf$comparisons)/log1p(count_max)>.65,'#FFFFFF',INK)
 for(i in which(positive))txt(cf$comparisons[i],xx[as.integer(cf$function_code[i])+1],yy[cf$crop_rank[i]],5.25,col=contrast[i],just='center')
 for(j in seq_len(40))txt(crops$display_label[j],x0-.018,yy[j],5.5,just='right')
 flabels<-c('N fixing','P solub.','N + P','Unclear','S/K solub.','AMF')
 for(j in seq_len(6))txt(flabels[j],xx[j],y0-.013,5.8,just='right',rot=60)
 txt('Frozen nominal function groups',.67,.036,6.1,just='center')
}

theme_set(theme_classic(base_size=7,base_family=FONT)+theme(
 axis.line=element_line(linewidth=.28,colour=INK),axis.ticks=element_line(linewidth=.25,colour=INK),
 axis.text=element_text(size=6.2,colour=INK),axis.title=element_text(size=6.8),
 plot.title=element_text(size=8.5,face='bold'),plot.subtitle=element_text(size=6.4,colour=MID),
 plot.margin=margin(3,4,3,3),panel.grid.minor=element_blank()))
pC<-ggplot(studies,aes(rank,comparisons))+
 geom_line(colour='#A9B7C0',linewidth=.35)+
 geom_point(colour=MID,size=.85)+
 geom_point(data=studies[1:5,],colour=ROLECOL[['Test']],size=1.2)+
 geom_rug(sides='r',alpha=.18,colour=MID,length=unit(.018,'npc'),linewidth=.3)+
 scale_x_continuous(breaks=c(1,50,100,169),limits=c(1,173),expand=expansion(mult=c(.015,.01)))+
 scale_y_log10(breaks=c(1,2,5,10,20,40,80),labels=c(1,2,5,10,20,40,80),limits=c(.85,95))+
 labs(title='C  Study evidence concentration',
      subtitle=sprintf('Largest study: %.1f%% · top five: %.1f%% of comparisons',100*max(studies$comparison_share),100*sum(studies$comparison_share[1:5])),
      x='Study rank by comparison count',y='Comparisons / study (log scale)')+
 theme(panel.grid.major.y=element_line(colour='#EEF0F2',linewidth=.22))

draw_D<-function(){
 header('D','Whole-study LOSO membership','6 representative folds of 169 · columns = studies',compact=TRUE)
 x0<-.31;x1<-.985;y0<-.25;y1<-.832;dx<-(x1-x0)/169;dy<-(y1-y0)/6
 xx<-x0+(seq_len(169)-.5)*dx;yy<-y1-(seq_len(6)-.5)*dy
 for(i in seq_along(d_folds)) {
   dd<-ddata[ddata$fold_id==d_folds[i],];rect(xx[dd$study_rank],rep(yy[i],169),dx*.95,dy*.68,ROLECOL[dd$membership])
   txt(sprintf('S%03d · n=%d',d_studies[i],studies$comparisons[d_ranks[i]]),x0-.022,yy[i],6.3,just='right')
 }
 for(r in c(1,50,100,169)){txt(r,xx[r],y0-.038,6,just='center');line(xx[r],y0-.012,xx[r],y0-.024,INK,.3)}
 txt('Study rank (same order as A)',.65,.116,6.6,just='center')
 txt('Train–test study overlap = 0 in all 169 folds',.31,.041,6.1)
}

draw_E<-function(){
 header('E','Crop-out with whole-study purge','5 representative folds of 40 · same 169 study columns',compact=TRUE)
 x0<-.225;x1<-.985;y0<-.445;y1<-.867;dx<-(x1-x0)/169;dy<-(y1-y0)/length(e_crops)
 xx<-x0+(seq_len(169)-.5)*dx
 for(i in seq_along(e_crops)) {
   yc<-y1-(i-.5)*dy
   for(b in c('Held-out crop','Other crops')) {
     dd<-edata[edata$fold_id==e_folds[i] & edata$band==b,];yb<-yc+if(b=='Held-out crop')dy*.18 else -dy*.18
     rect(xx[dd$study_rank],rep(yb,169),dx*.95,dy*.32,ROLECOL[dd$membership],border=ifelse(dd$membership=='Absent','#EDF0F2',NA),lwd=.12)
     txt(if(b=='Held-out crop')'H' else 'O',x0-.012,yb,5.1,just='right',col=MID)
   }
   txt(crops$display_label[match(e_crops[i],crop_order)],x0-.045,yc,6.2,just='right')
 }
 for(r in c(1,50,100,169))txt(r,xx[r],y0-.018,5.8,just='center')
 txt('H: held-out crop; O: other crops within each study',x0,.382,5.9)
 # A real mixed-crop study from the maximum-purge displayed fold.
 label<-study_meta$study_label[match(detail_study,study_meta$study_id)]
 held<-loco$held_out_level[match(detail_fold,loco$fold_id)]
 nt<-sum(detail_data$membership=='Test');np<-sum(detail_data$membership=='Purge');total<-nt+np
 line(.02,.343,.985,.343,RULE,.45)
 txt(sprintf('Mixed-crop study S%03d (%s)',detail_study,label),.025,.303,6.7,face='bold')
 txt(sprintf('Holding out %s excludes all %d rows from training',held,total),.025,.260,6.1)
 bx0<-.025;bw<-.945;by<-.165;bh<-.066;frac<-nt/total
 rect(bx0+bw*frac/2,by,bw*frac,bh,ROLECOL[['Test']])
 rect(bx0+bw*frac+bw*(1-frac)/2,by,bw*(1-frac),bh,ROLECOL[['Purge']])
 txt(sprintf('Test: %s (n=%d)',held,nt),bx0,.219,5.8)
 txt(sprintf('Purge: other crops (n=%d)',np),.47,.219,5.8)
 txt(nt,bx0+bw*frac/2,by,5.8,col='white',just='center')
 txt(np,bx0+bw*frac+bw*(1-frac)/2,by,5.8,col='white',just='center')
 txt('Held-out-crop rows are tested; other rows from those studies are purged.',.025,.067,5.9)
}

draw_F<-function(){
 header('F','Potential shared-control dependence','Repeated control signatures; physical sharing unverified',compact=TRUE)
 n<-nrow(fexample);ys<-seq(.81,.46,length.out=n);cx<-.16;cy<-mean(ys);tx<-.58
 # Dashed links deliberately encode a candidate relationship, not a proven replicate identity.
 line(rep(cx+.038,n),rep(cy,n),rep(tx-.021,n),ys,MID,.55,lty=2)
 grid.circle(x=unit(cx,'npc'),y=unit(cy,'npc'),r=unit(2.1,'mm'),gp=gpar(fill='white',col=INK,lwd=.6));txt('C*',cx,cy,6.8,just='center')
 point(rep(tx,n),ys,1.3,ROLECOL[['Test']])
 for(i in seq_len(n)){
   txt(paste0('Pair ',fexample$pair_id[i]),tx+.035,ys[i]+.012,6.1)
   txt(paste0('coded C row ',fexample$control_excel_row[i]),tx+.035,ys[i]-.031,5.5,col=MID)
 }
 txt('Same control',.16,.442,6.1,just='center');txt('signature',.16,.407,6.1,just='center')
 txt(paste0('Example: ',study_meta$study_label[match(fexample$study_id[1],study_meta$study_id)]),.025,.354,6.4)
 line(.02,.311,.985,.311,RULE,.45)
 txt('315 suspected clusters',.025,.261,7.1,face='bold')
 txt(sprintf('1,117 / 1,723 comparisons (%.1f%%)',1117/1723*100),.025,.211,6.6)
 txt('90 / 169 studies contain a candidate cluster',.025,.160,6.3)
 txt('Minimum safe isolation unit: whole study',.025,.080,6.8,face='bold')
 txt('C* is a matching signature, not a confirmed physical control.',.025,.031,5.6,col=MID)
}

draw_key<-function(){
 # A single shared, quantitative support key; country encoding is display-only.
 txt('Comparisons',.013,.78,6.2,face='bold')
 lx<-.16;lw<-.21;ncol<-100;ii<-seq(0,count_max,length.out=ncol)
 # Equal visual steps in log(1+n), with raw-count tick labels.
 vals<-expm1(seq(0,log1p(count_max),length.out=ncol))
 rect(lx+(seq_len(ncol)-.5)*lw/ncol,rep(.77,ncol),lw/ncol,.19,count_color(vals))
 ticks<-unique(c(0,1,5,20,count_max));ticks<-ticks[ticks<=count_max]
 for(v in ticks){x<-lx+lw*log1p(v)/log1p(count_max);txt(v,x,.48,5.4,just='center')}
 txt('Train',.445,.78,6);rect(.427,.78,.015,.20,ROLECOL[['Train']])
 txt('Test',.522,.78,6);rect(.504,.78,.015,.20,ROLECOL[['Test']])
 txt('Purge',.597,.78,6);rect(.579,.78,.015,.20,ROLECOL[['Purge']])
 txt('Absent',.685,.78,6);rect(.667,.78,.015,.20,ROLECOL[['Absent']],RULE)
 txt('H / O membership bands in E',.775,.78,5.6,col=MID)
 txt('Country',.013,.15,6.1,face='bold')
 for(i in seq_along(country_levels)){
   x<-.16+(i-1)*.111;rect(x,.15,.014,.2,country_col[[country_levels[i]]]);txt(country_levels[i],x+.014,.15,5.7)
 }
}

# Capture grid panels using the same R backend and physical aspect ratios as the layout.
gA<-capture_panel(draw_A,4.65,4.8);gB<-capture_panel(draw_B,2.35,4.8)
gD<-capture_panel(draw_D,3.5,1.8);gE<-capture_panel(draw_E,3.5,2.1);gF<-capture_panel(draw_F,3.5,2.1);gkey<-capture_panel(draw_key,7.1,.37)
wA<-wrap_elements(full=gA,clip=FALSE);wB<-wrap_elements(full=gB,clip=FALSE)
wD<-wrap_elements(full=gD,clip=FALSE);wE<-wrap_elements(full=gE,clip=FALSE);wF<-wrap_elements(full=gF,clip=FALSE);wkey<-wrap_elements(full=gkey,clip=FALSE)
top<-(wA|wB)+plot_layout(widths=c(1.9,1))
middle<-(pC|wD)+plot_layout(widths=c(1,1))
bottom<-(wE|wF)+plot_layout(widths=c(1,1))
fig<-(top/middle/bottom/wkey)+plot_layout(heights=c(125,46,57,10))
fig<-fig & theme(plot.background=element_rect(fill='white',colour=NA))
width_mm<-183;height_mm<-238
pdf_file<-file.path(out,'Figure1_draft.pdf');png_file<-file.path(out,'Figure1_draft.png')
grDevices::cairo_pdf(pdf_file,width=width_mm/25.4,height=height_mm/25.4,family=FONT,onefile=TRUE)
print(fig);invisible(dev.off())
ragg::agg_png(png_file,width=width_mm,height=height_mm,units='mm',res=600,background='white')
print(fig);invisible(dev.off())
# A lightweight, identical R rendering for visual inspection; final PNG remains 600 dpi.
ragg::agg_png(file.path(data_dir,'Figure1_QA_preview.png'),width=width_mm,height=height_mm,units='mm',res=180,background='white')
print(fig);invisible(dev.off())
# Panel QA is exported by R as well, never cross-rendered by another language.
ragg::agg_png(file.path(data_dir,'Figure1_QA_panels_A_B.png'),width=width_mm,height=125,units='mm',res=300,background='white')
print(top);invisible(dev.off())
ragg::agg_png(file.path(data_dir,'Figure1_QA_panels_D_E_F.png'),width=width_mm,height=110,units='mm',res=300,background='white')
print((wD|plot_spacer())/(wE|wF)+plot_layout(heights=c(46,57)));invisible(dev.off())
after<-vapply(file.path(base,source_paths),hash,character(1));names(after)<-names(source_paths)
stopifnot(identical(before,after))
# PNG header: dimensions/resolution metadata checked without any graphic conversion.
con<-file(png_file,'rb');signature<-readBin(con,'raw',8);stopifnot(all(as.integer(signature)==c(137,80,78,71,13,10,26,10)))
png_dims<-NULL;ppm<-NULL
repeat{
 len<-readBin(con,'integer',1,size=4,endian='big');if(!length(len))break
 typ<-rawToChar(readBin(con,'raw',4));raw<-readBin(con,'raw',len);readBin(con,'raw',4)
 if(typ=='IHDR')png_dims<-readBin(raw[1:8],'integer',2,size=4,endian='big')
 if(typ=='pHYs')ppm<-readBin(raw[1:8],'integer',2,size=4,endian='big')
 if(typ=='IEND')break
};close(con)
stopifnot(all(abs(png_dims-c(width_mm,height_mm)/25.4*600)<2),all(abs(ppm*.0254-600)<.02))
qc<-list(cohort=c(comparisons=1723,studies=169,crops=40,functions=6),input_unchanged=TRUE,
         source_paths=as.list(source_paths),source_sha256=as.list(before),all_215_folds_validated=TRUE,
         train_test_and_train_purge_study_overlap=0,source_counts_reconciled=TRUE,
         shared_control_status='suspected; not confirmed physical controls',candidate_clusters=315,candidate_comparisons=1117,
         representative_LOSO_held_studies=as.list(d_studies),representative_crop_holdouts=as.list(e_crops),
         crop_detail=list(fold=detail_fold,study=detail_study),shared_control_example=fcluster,
         no_model_training=TRUE,no_performance_files_read=TRUE,no_outcome_reanalysis=TRUE,
         backend='R: ggplot2 + grid + patchwork',font=FONT,width_mm=width_mm,height_mm=height_mm,
         PNG_pixels=as.list(png_dims),PNG_dpi=as.list(ppm*.0254),
         output_sha256=list(pdf=hash(pdf_file),png=hash(png_file)),visual_review='pending inspection')
jsonlite::write_json(qc,file.path(data_dir,'figure1_QC.json'),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(data_dir,'R_sessionInfo.txt'))
caption<-c('Figure 1 | Evidence structure and leakage-resistant validation framework.',
  'A, All 169 studies by 40 frozen crop groups; tile colour encodes comparison count on a log(1+n) display scale. Studies are sorted by decreasing comparison count (ties: numeric study ID); crops by decreasing unique-study support (ties: comparison count, then crop label). The country strip displays the five highest-study-support countries, Other and Multiple; exact countries are retained in source data. Aligned points show study sizes and crop study support.',
  'B, Comparison support for each crop by the six frozen nominal inoculant-function groups. Zero cells are retained; counts are printed in nonzero cells. Crop order and the count colour scale match A. Unclear and combined-function categories retain their original meaning.',
  'C, Ranked study comparison counts, with a right-axis rug showing the size distribution. The five largest studies are highlighted; no outcome, effect estimate, uncertainty interval or significance test is shown.',
  'D, Six illustrative LOSO folds selected by study-size ranks 1, 2, 3, 85, 168 and 169, showing all 169 study columns in A order. Entire test studies are isolated; all 169 LOSO folds were validated against the frozen manifest.',
  'E, Five crop-out folds selected from the four most study-supported crops and the largest-purge remaining crop. Each study column has a held-out-crop band (H) and other-crop band (O); white indicates absence of that class, not missing membership. Held-out-crop comparisons are tested and other comparisons in the same studies are purged. The inset shows an actual mixed-crop study from the largest-purge displayed fold. All 40 crop-out folds satisfy whole-study isolation.',
  'F, A source-identified candidate repeated-control signature linked by dashed lines to its coded pairs. The relationship is suspected, not proof of a single physical control or recoverable plot-level covariance. Overall, 315 candidate clusters contain 1,117 comparisons in 90 studies. The minimum safe separation unit remains the whole study.',
  'Only the frozen 1,723-comparison / 169-study primary cohort is used. No models were trained or refitted, no performance files were consulted and no scientific result/input was altered. Colours for Train/Test/Purge are fixed across D/E. Country is an annotation, not a new predictor. All displayed values and exact membership IDs are available in figure1_plot_data.')
writeLines(caption,file.path(data_dir,'Figure1_draft_caption.txt'),useBytes=TRUE)
message('Figure 1 draft exported; sources unchanged; PNG ',paste(png_dims,collapse=' × '),' pixels; 600 dpi.')
