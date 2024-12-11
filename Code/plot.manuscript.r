

# Nicholas Ducharme-Barth
# 2024/12/04
# Using 2019 skipjack diagnostic model
# Plot models that induce effort creep

#_____________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)
	library(FLR4MFCL)
	library(ggplot2)
	library(ggthemes)
	library(diags4MFCL)
	library(frqit)
	library(sf)
	library(patchwork)
	library(ggchicklet)


#_____________________________________________________________________________________________________________________________
# set working directory
	proj_dir = this.path::this.proj()
	dir_model = paste0(proj_dir,"/Model_Runs/")

#_____________________________________________________________________________________________________________________________
# make plots
	dir_plot = paste0(proj_dir,"/Plot/")
	dir.create(dir_plot,recursive=TRUE,showWarnings=FALSE)

		dir.list = list.dirs(dir_model)
		dir.list = dir.list[grep("tmp",dir.list)]

		rep.list= lapply(c(paste0(proj_dir,"/Model_Runs/Length100_GrowthDiag_Mix1_H0.8/plot-07.par.rep"),paste0(dir.list,"/plot-tmp.par.rep")),read.MFCLRep)
		par.list = lapply(c(paste0(proj_dir,"/Model_Runs/Length100_GrowthDiag_Mix1_H0.8/07.par"),paste0(dir.list,"/tmp.par")),read.MFCLPar,first.yr=1972)
		# frq.list= lapply(paste0(save.dir.vec,"swo.frq"),read.MFCLFrq)
		model_name = c("2019 Diagnostic","MFCL v2084 - 0%/yr","MFCL v2084 - 1.5%/yr","MFCL v2084 - 1%/yr","MFCL v2084 - 2.5%/yr","MFCL v2084 - 2%/yr","MFCL v2084 - 3%/yr","MFCL v2084 - 6%/yr")

		names(par.list) = names(rep.list) = model_name


#_____________________________________________________________________________________________________________________________
# spatial plot
	pb_coast = st_read(paste0(proj_dir,"/background-data/pb_coast.shp"))
	skj_2019 = st_read(paste0(proj_dir,"/background-data/skj_2019.shp"))

 	labels_dt = data.table(x=c(130,180,135,180,125,150,160,190),
						   y=c(40,40,20,20,0,-10,5,-5),
						   label=paste0("Region ",1:8))

	p1 = ggplot() + 
	# ggtitle("2019 WCPFC skipjack tuna") +
	xlab("") +
	ylab("") +
	geom_sf(data=skj_2019,color="#323C46",fill="NA",linewidth=0.5) +
	geom_sf(data=pb_coast,color="#001743",fill="#002364") +
    coord_sf(xlim = c(100, 220), ylim = c(-30, 60), expand = FALSE) +
	geom_label(data=labels_dt,aes(x=x,y=y,label=label),color="#001743", size = 3) +
	theme(panel.grid.major = element_line(color = "white", linetype = "dashed",linewidth=0.5),
		  panel.background = element_rect(fill = "#e9f3f6"),
          panel.border = element_rect(colour = "black", fill=NA))
    # ggsave("skj-regions.png",plot=p1, device = "png", path = "Plot/",
	#   			scale = 0.75, width =8, height = 6, units = c("in"),
	#   			dpi = 300, limitsize = TRUE)  

#_____________________________________________________________________________________________________________________________
# plot cpue


    effort_creep_vec = c(1,1.5,2,2.5,3,6)
	obs_quant = 0.5*(1-95/100)
	frq = readfrq(paste0(proj_dir,"/Model_Runs/Length100_GrowthDiag_Mix1_H0.8/skj.frq"))
	base_cpue = as.data.table(cateffpen(frq)) %>%
		 			   .[,id:=paste0(year,"_",month,"_",fishery)] %>%
					   .[,.(id,year,month,week,fishery,catch,effort,penalty)] %>%
					   .[,ts:=as.numeric(year)+(as.numeric(month)-2)/12] %>%
					   .[,cpue:=catch/effort] %>%
                       .[fishery %in% c(1,4,7,22,24,28,12,19)] %>%
					   .[effort>=0] %>%
					   .[catch>=0] %>%
                       .[penalty>0] %>%
					   .[cpue>0] %>%
                       .[,sd:=1/sqrt(2*penalty)] %>%
                       .[,mean_penalty:=mean(penalty),by=fishery] %>%
                       .[,mean_sd:=mean(sd),by=fishery] %>%
					   .[,cv:=sd/cpue] %>%
					   .[,cpue:=cpue/mean(cpue),by=fishery] %>%
                       .[,upper:=qlnorm(1-obs_quant,meanlog=log(cpue),sdlog=cv)] %>%
                       .[,lower:=qlnorm(obs_quant,meanlog=log(cpue),sdlog=cv)] %>%
					   .[,effort_creep:=factor(0,levels=c(0,effort_creep_vec),labels=paste0(c(0,effort_creep_vec)," %"))]

	# adjust cpue
	effort_creep_dt.list = as.list(rep(NA,length(effort_creep_vec))) 

	for(i in seq_along(effort_creep_vec)){

			catch_dt = as.data.table(cateffpen(frq)) %>%
		 			   .[,id:=paste0(year,"_",month,"_",fishery)] %>%
					   .[,.(id,year,month,week,fishery,catch,effort,penalty)] %>%
					   .[,ts:=as.numeric(year)+(as.numeric(month)-2)/12] %>%
					   .[,e_scalar:=1]
			
			idx_fish = c(1,4,7,22,24,28,12,19)
			for(f in idx_fish)
			{	
				# define scalar to modify effort column and add to catch_dt
					tmp_catch = catch_dt[f==fishery&catch>0&effort>0&ts>1980,]
					tmp_range = range(tmp_catch$ts)
					tmp_sq = seq(from=tmp_range[1],to=tmp_range[2],by=0.25)
					tmp_scalar = rep(NA,length(tmp_sq))
					tmp_scalar[1] = 1
					for(s in 2:length(tmp_scalar))
					{
						tmp_scalar[s] = tmp_scalar[s-1]*(1+(effort_creep_vec[i]/400))
					}
					tmp_catch$e_scalar = tmp_scalar[match(catch_dt[f==fishery&catch>0&effort>0&ts>1980,]$ts,tmp_sq)]
					catch_dt = catch_dt[!(id %in% tmp_catch$id)]
					catch_dt = rbind(catch_dt,tmp_catch) %>% .[order(fishery,year,month)]

				# clean-up
				rm(list=c("tmp_catch","tmp_range","tmp_sq","tmp_scalar"))
			}

						# modify effort with scalar
				catch_dt$effort = catch_dt$effort * catch_dt$e_scalar

			# reformat 
				effort_creep_dt.list[[i]] = catch_dt %>%
							.[order(fishery,year,month)] %>%
							.[,id:=NULL] %>%
							.[,ts:=NULL] %>%
							.[,e_scalar:=NULL] %>%
							.[,effort_creep:=effort_creep_vec[i]]
	}

	effort_creep_dt = rbindlist(effort_creep_dt.list) %>%
		 			   .[,id:=paste0(effort_creep,"_",fishery)] %>%
					   .[,.(effort_creep,id,year,month,week,fishery,catch,effort,penalty)] %>%
					   .[,ts:=as.numeric(year)+(as.numeric(month)-2)/12] %>%
					   .[,cpue:=catch/effort] %>%
                       .[fishery %in% c(1,4,7,22,24,28,12,19)] %>%
					   .[effort>=0] %>%
					   .[catch>=0] %>%
                       .[penalty>0] %>%
					   .[cpue>0] %>%
					   .[,cpue:=cpue/mean(cpue),by=.(id)] %>%
					   .[,.(cpue=mean(cpue)),by=.(id,effort_creep,fishery,year)] %>%
					   .[,effort_creep:=factor(effort_creep,levels=c(0,effort_creep_vec),labels=paste0(c(0,effort_creep_vec)," %"))] %>%
					   .[,fishery:=factor(fishery,levels=c(1,4,7,22,12,19,24,28),labels=c("Region 1: JP PL","Region 2: JP PL","Region 3: JP PL","Region 4: JP PL","Region 5: ID-PH PS","Region 6: PS-A","Region 7: JP PL","Region 8: JP PL"))]
	
	labels_dt = data.table(x=rep(1972,8),
						   y=rep(3,8)) %>%
						  .[,fishery:=factor(c(1,4,7,22,12,19,24,28),levels=c(1,4,7,22,12,19,24,28),labels=c("Region 1: JP PL","Region 2: JP PL","Region 3: JP PL","Region 4: JP PL","Region 5: ID-PH PS","Region 6: PS-A","Region 7: JP PL","Region 8: JP PL"))] %>%
						  .[,label:=levels(fishery)]

	scenario_colors = c("#001743",viridis::viridis(6, begin = 0.2,end = 0.8,direction = 1,option = "H"))

	p2 = base_cpue %>%
	.[,.(cpue=mean(cpue),lower=mean(lower),upper=mean(upper)),by=.(effort_creep,fishery,year)] %>%
	.[,fishery:=factor(fishery,levels=c(1,4,7,22,12,19,24,28),labels=c("Region 1: JP PL","Region 2: JP PL","Region 3: JP PL","Region 4: JP PL","Region 5: ID-PH PS","Region 6: PS-A","Region 7: JP PL","Region 8: JP PL"))] %>%
	ggplot() +
	ylab("") +
	xlab("") +
	scale_x_continuous(breaks=c(1975,1995,2015)) +
    facet_wrap(~fishery) +
    geom_hline(yintercept=1,linetype="dashed") +
	geom_hline(yintercept=0) +
    geom_path(aes(x=year,y=cpue,color=effort_creep),linewidth=1.25) +
	geom_path(data=effort_creep_dt,aes(x=year,y=cpue,color=effort_creep,group=id)) +
	geom_ribbon(aes(x=year,ymin=lower,ymax=upper),alpha=0.2,fill="#002364") +
	geom_label(data=labels_dt,aes(x=x,y=y,label=label),vjust="top",hjust="left",color="#001743", size = 3) +
	scale_color_manual("Effort\ncreep\nscenario",values=scenario_colors,drop=FALSE) +
	scale_fill_manual("Effort\ncreep\nscenario",values=scenario_colors,drop=FALSE) +	
	theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
							panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
							panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
							strip.background =element_rect(fill="white"),
							legend.key = element_rect(fill = "white", color = "black", linetype = "solid")) +
							theme(strip.text = element_blank())
							# theme(text = element_text(size = 20))
    # ggsave("skj-idx.png",plot=p2, device = "png", path = "Plot/",
	#   			scale = 1, width =12, height = 9, units = c("in"),
	#   			dpi = 300, limitsize = TRUE) 


#_____________________________________________________________________________________________________________________________
# extract quantities
    # the relative recruitment (overall rec dev)
	diag_relrec = as.data.table(rel_rec(par.list[[1]])) %>%
				  .[,year:=as.numeric(year)] %>%
                  .[,ts:=as.numeric(year)+(as.numeric(season)-1)/4] %>%
                  .[,type:="relrec"] %>%
                  .[,.(type,year,ts,value)] %>%
				  .[,.(value=mean(value)),by=.(type,year)] %>%
			      .[,effort_creep:=factor(0,levels=c(0,effort_creep_vec),labels=paste0(c(0,effort_creep_vec)," %"))]


	rel_rec_dt.list = as.list(rep(NA,length(par.list)-2))
	for(i in seq_along(rel_rec_dt.list)){
		rel_rec_dt.list[[i]] = as.data.table(rel_rec(par.list[[i+2]])) %>%
				  .[,year:=as.numeric(year)] %>%
                  .[,ts:=as.numeric(year)+(as.numeric(season)-1)/4] %>%
                  .[,type:="relrec"] %>%
                  .[,.(type,year,ts,value)] %>%
				  .[,.(value=mean(value)),by=.(type,year)] %>%
				  .[,model:=names(par.list)[i+2]] %>%
				  .[,.(model,type,year,value)]
	}
	rel_rec_dt = rbindlist(rel_rec_dt.list) %>%
				 .[,effort_creep:=gsub("MFCL v2084 - ","",model)] %>%
				 .[,effort_creep:=as.numeric(gsub("%/yr","",effort_creep))] %>%
				 .[,effort_creep:=factor(effort_creep,levels=effort_creep_vec,labels=paste0(effort_creep_vec," %"))]


    # all catch in mt, exclude LL fleets (3,6,9,17,21,23,27,31)
    diag_catch = as.data.table(catch_pred(rep.list[[1]])) %>%
                 .[!(unit%in%c("3","6","9","17","21","23","27","31"))] %>%
                 .[,.(value=sum(value)),by=.(year,season)] %>%
				 .[,year:=as.numeric(year)] %>%
                 .[,ts:=as.numeric(year)+(as.numeric(season)-1)/4] %>%
                  .[,type:="catch"] %>%
                  .[,.(type,year,value)] %>%
				  .[,.(value=sum(value)/1000),by=.(type,year)]
	
	scale_catch = 1000 # mean(diag_catch$value)

	

	effort_creep_levels = levels(rel_rec_dt$effort_creep)
	pval_vec = slope_vec = rep(NA,length(effort_creep_levels)+1)
	
	slope_vec[1] = round(lm(data=diag_relrec,value~year)$coefficients[2],digits=3)
	pval_vec[1] = summary(lm(data=diag_relrec,value~year))$coefficients[2,4]
	for(i in seq_along(effort_creep_levels)){
		slope_vec[i+1] = round(lm(data=rel_rec_dt[effort_creep==effort_creep_levels[i]],value~year)$coefficients[2],digits=3)
		pval_vec[i+1] = summary(lm(data=rel_rec_dt[effort_creep==effort_creep_levels[i]],value~year))$coefficients[2,4]
	}

	labels_dt = data.table(x=c(1991,2007,1991,2007,1991,2007,1991,2007)-1,
						   y=c(0.55,0.55,0.425,0.425,0.3,0.3,0.175,0.175)-0.025,
						   label=c("Slope",paste0(c(0,effort_creep_vec)," %: ",slope_vec)),
						   bold=c(FALSE,pval_vec>0.05))
	
	segment_dt = data.table(x=c(2004.5,1988.5,2004.5,1988.5,2004.5,1988.5,2004.5)-1.25,
							xend=c(2006.5,1990.5,2006.5,1990.5,2006.5,1990.5,2006.5)-1.25,
							y=c(0.4875,0.3625,0.3625,0.2375,0.2375,0.1125,0.1125),
							effort_creep=factor(c(0,effort_creep_vec),levels=c(0,effort_creep_vec),labels=paste0(c(0,effort_creep_vec)," %")))
	
	p3 = ggplot() +
	xlab("") +
	scale_x_continuous(breaks=c(1975,1995,2015)) +
	scale_y_continuous(name = "Relative recruitment", sec.axis = sec_axis( transform=~.*scale_catch, name="Catch (1000s mt)")) +
	geom_hline(yintercept=1,linetype="dashed") +
	geom_chicklet(data=diag_catch,aes(x=year,y=value/scale_catch),fill="#A5AAAF") +
	geom_hline(yintercept=0) +
	geom_path(data=diag_relrec,aes(x=year,y=value,color=effort_creep),linewidth=1.25) +
	geom_path(data=rel_rec_dt,aes(x=year,y=value,color=effort_creep)) +
	geom_rrect(data=data.table(xmin=1986,xmax=2018,ymin=0.025,ymax=0.575),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="white",color="#001743") +
	geom_text(data=labels_dt[bold==FALSE],aes(x=x,y=y,label=label),hjust="left",vjust="top") +
	geom_text(data=labels_dt[bold==TRUE],aes(x=x,y=y,label=label),hjust="left",vjust="top",fontface="bold") +
	geom_segment(data=segment_dt[effort_creep=="0 %"],aes(x=x,xend=xend,y=y),color=scenario_colors[1],alpha=1,linewidth=1.25) +
	geom_segment(data=segment_dt[effort_creep!="0 %"],aes(x=x,xend=xend,y=y),color=scenario_colors[-1],alpha=1) +
	# geom_smooth(data=diag_relrec,aes(x=year, y=value),color="#001743",method="lm",linewidth=1.15) +
	# geom_smooth(data=rel_rec_dt,aes(x=year, y=value,color=effort_creep),method="lm") +
	scale_color_manual("Effort\ncreep\nscenario",values=scenario_colors,drop=FALSE) +
	scale_fill_manual("Effort\ncreep\nscenario",values=scenario_colors,drop=FALSE) +
	theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
							panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
							panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
							strip.background =element_rect(fill="white"),
							legend.key = element_rect(fill = "white", color = "black", linetype = "solid")) +
							theme(strip.text = element_blank())

#_____________________________________________________________________________________________________________________________
# majuro plot

	diag_sbfo_dt = as.data.table(adultBiomass_nofish(rep.list[[1]])) %>%
	setnames(.,"value","sbfo") %>%
	.[,.(sbfo=sum(sbfo)),by=.(year,season)] %>%
	.[,.(sbfo=mean(sbfo)),by=year]

	diag_sb_dt = as.data.table(adultBiomass(rep.list[[1]])) %>%
	setnames(.,"value","sb") %>%
	.[,.(sb=sum(sb)),by=.(year,season)] %>%
	.[,.(sb=mean(sb)),by=year]

	diag_majuro_dt = data.table(effort_creep=factor(0,levels=c(0,effort_creep_vec),labels=paste0(c(0,effort_creep_vec)," %")))
	diag_majuro_dt$sbsbfo = mean(diag_sb_dt[year%in%as.character(2015:2018)]$sb)/mean(diag_sbfo_dt[year%in%as.character(2008:2017)]$sbfo)
	diag_majuro_dt$ffmsy = mean(rowMeans(FFMSY_ts(rep.list[[1]])@.Data[1,,1,,,1])[c("2014","2015","2016","2017")])

	majuro_dt.list = as.list(rep(NA,length(rep.list)-2))
	for(i in seq_along(majuro_dt.list)){
		diag_sbfo_dt = as.data.table(adultBiomass_nofish(rep.list[[i+2]])) %>%
		setnames(.,"value","sbfo") %>%
		.[,.(sbfo=sum(sbfo)),by=.(year,season)] %>%
		.[,.(sbfo=mean(sbfo)),by=year]

		diag_sb_dt = as.data.table(adultBiomass(rep.list[[i+2]])) %>%
		setnames(.,"value","sb") %>%
		.[,.(sb=sum(sb)),by=.(year,season)] %>%
		.[,.(sb=mean(sb)),by=year]

		majuro_dt.list[[i]] = data.table(effort_creep=names(rep.list)[i+2]) %>%
				 .[,effort_creep:=gsub("MFCL v2084 - ","",effort_creep)] %>%
				 .[,effort_creep:=as.numeric(gsub("%/yr","",effort_creep))] %>%
				 .[,effort_creep:=factor(effort_creep,levels=effort_creep_vec,labels=paste0(effort_creep_vec," %"))]
		majuro_dt.list[[i]]$sbsbfo = mean(diag_sb_dt[year%in%as.character(2015:2018)]$sb)/mean(diag_sbfo_dt[year%in%as.character(2008:2017)]$sbfo)
		majuro_dt.list[[i]]$ffmsy = mean(rowMeans(FFMSY_ts(rep.list[[i+2]])@.Data[1,,1,,,1])[c("2014","2015","2016","2017")])
	}
	majuro_dt = rbind(diag_majuro_dt,rbindlist(majuro_dt.list))[order(effort_creep)] 


	p4 = ggplot() +
		xlab(expression(SB["2015-2018"]/SB["F=0"])) +
		ylab(expression(F["2014-2017"]/F["MSY"])) +
		coord_cartesian(xlim=c(0,1),ylim=c(0,2)) +
        scale_x_continuous(expand = expansion(mult=c(0,0)),breaks=c(0.2,0.4,0.6,0.8)) +
        scale_y_continuous(expand = expansion(mult=c(0,0))) +
        geom_polygon(data=data.table(SB_SBmsy=c(0,0.2,0.2,0,0),F_Fmsy=c(0,0,1,1,0)),aes(x=SB_SBmsy,y=F_Fmsy),fill="red",alpha=0.2) +
        geom_polygon(data=data.table(SB_SBmsy=c(0,0.2,0.2,0,0),F_Fmsy=c(1,1,5,5,1)),aes(x=SB_SBmsy,y=F_Fmsy),fill="red",alpha=0.2) +
        geom_polygon(data=data.table(SB_SBmsy=c(0.2,5,5,0.2,0.2),F_Fmsy=c(1,1,5,5,1)),aes(x=SB_SBmsy,y=F_Fmsy),fill="orange",alpha=0.2) +
        geom_polygon(data=data.table(SB_SBmsy=c(0.2,5,5,0.2,0.2),F_Fmsy=c(0,0,1,1,0)),aes(x=SB_SBmsy,y=F_Fmsy),fill="white",alpha=0.2) +
        geom_hline(yintercept=0,color="black") +
		geom_vline(xintercept=0,color="black") +
        geom_hline(yintercept=1,linewidth=1.15,color="black") +
		geom_vline(xintercept=0.2,linewidth=1.15,color="black") +
		geom_point(data=majuro_dt,aes(x=sbsbfo,y=ffmsy,fill=effort_creep),shape=21,color="black",size=3) +
		scale_color_manual("Effort\ncreep\nscenario",values=scenario_colors,drop=FALSE) +
		scale_fill_manual("Effort\ncreep\nscenario",values=scenario_colors,drop=FALSE) +			
		theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
							panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
							panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
							strip.background =element_rect(fill="white"),
							legend.key = element_rect(fill = "white"))
#_____________________________________________________________________________________________________________________________
# combine plots using patchwork
	p = p1 + p2 + p3 + p4 + 
	plot_layout(guides = 'collect') +
	plot_annotation(tag_levels = 'A',tag_suffix=")")
    ggsave("skj-effort-creep.png",plot=p, device = "png", path = "Plot/",
	  			scale = 1, width =12, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE) 

