

# Nicholas Ducharme-Barth
# 28/03/2022
# 2022 SPC PAW
# Using 2019 skipjack diagnostic model
# Explore effect of inducing effort creep (0%,1%,3%,6%) per year

#_____________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)
	library(FLR4MFCL)
	library(ggplot2)
	library(ggthemes)
	library(diags4MFCL)
	library(frqit)

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

		agg.years = TRUE
		agg.regions=TRUE
		biomass.type = "SSB"
		biomass.units=1000
		recdist.year_range=NULL
		yaxis.free=FALSE
		LRP=0.20
		TRP=NULL
		turbo_pal=function(selected.model.names, all.model.names=selected.model.names)
		{
			n=length(all.model.names)-1
			out = viridis::viridis_pal(alpha = 1, begin = 0.1, end = 0.8, direction = 1, option = "H")(n)
			out = c("black",out)[1:length(all.model.names)]
		    names(out) = all.model.names
		    out = out[selected.model.names]
		    return(out)
		}

        # reorder
        rep.list = rep.list[c(1,2,4,3,6,5,7,8)]
        par.list = par.list[c(1,2,4,3,6,5,7,8)]

	# biomass
			p = plot.biomass(rep.list,model_name,agg.years,agg.regions=TRUE,biomass.type=biomass.type,biomass.units=biomass.units,palette.func=turbo_pal)
			p = p + theme_few(base_size = 20) 
			ggsave(filename=paste0("bio.png"), plot = p, device = "png", path = dir_plot,scale = 1, width = 16, height = 9, units = c("in"),dpi = 300, limitsize = TRUE)

	# biomass-region
			p = plot.biomass(rep.list,model_name,agg.years,agg.regions=FALSE,biomass.type=biomass.type,biomass.units=biomass.units,palette.func=turbo_pal,yaxis.free=TRUE)
			p = p + theme_few(base_size = 20) 
			ggsave(filename=paste0("bio.reg.png"), plot = p, device = "png", path = dir_plot,scale = 1, width = 16, height = 9, units = c("in"),dpi = 300, limitsize = TRUE)

	# depletion
			p = plot.depletion(rep.list,model_name,agg.years,agg.regions=TRUE,biomass.type=biomass.type,LRP=LRP,TRP=TRP,palette.func=turbo_pal)
			p = p + theme_few(base_size = 20) 
			ggsave(filename=paste0("dep.png"), plot = p, device = "png", path = dir_plot,scale = 1, width = 16, height = 9, units = c("in"),dpi = 300, limitsize = TRUE)
			
	# depletion-region
			p = plot.depletion(rep.list,model_name,agg.years,agg.regions=FALSE,biomass.type=biomass.type,LRP=LRP,TRP=TRP,palette.func=turbo_pal)
			p = p + theme_few(base_size = 20) 
			ggsave(filename=paste0("dep.reg.png"), plot = p, device = "png", path = dir_plot,scale = 1, width = 16, height = 9, units = c("in"),dpi = 300, limitsize = TRUE)
				
	# recruitment
		rec_dt = rbindlist(lapply(1:8,function(x)
				{out=as.data.table(rec_region(rep.list[[x]])) %>%
				 .[,.(rec=sum(value)/1000000),by=.(year,season)] %>%
				 .[,.(rec=mean(rec)),by=year] %>%
				 .[,year:=as.numeric(year)] %>%
				 .[,model:=model_name[x]]
				 return(out)}
				 ))
		p = rec_dt %>%
			ggplot() +
			xlab("Year") + ylab("Recruits (millions)") +
			geom_hline(yintercept = 0) +
	  		geom_path( aes(x=year, y=rec, color=model),linewidth=1,alpha=0.6) +
	  		geom_smooth(aes(x=year, y=rec, color=model),method="lm",linewidth=1.15) +
	  		expand_limits(y=c(0)) +
	  		theme_few(base_size = 20) +
	  		scale_color_manual("Model",values = c("black",viridis::viridis_pal(alpha = 1, begin = 0.1, end = 0.8, direction = 1, option = "H")(7)))
			ggsave(filename=paste0("rec.png"), plot = p, device = "png", path = dir_plot,scale = 1, width = 9, height = 9, units = c("in"),dpi = 300, limitsize = TRUE)
	

	# recruitment-region

			rec_dt = rbindlist(lapply(1:8,function(x)
				{out=as.data.table(rec_region(rep.list[[x]])) %>%
				 .[,.(rec=mean(value)/1000000),by=.(year,area)] %>%
				 .[,year:=as.numeric(year)] %>%
				 .[,model:=model_name[x]]
				 return(out)}
				 ))
		p = rec_dt %>%
			ggplot() +
			xlab("Year") + ylab("Recruits (millions)") + facet_wrap(~area,scales="free_y") + 
			geom_hline(yintercept = 0) +
	  		geom_path( aes(x=year, y=rec, color=model),size=1,alpha=0.6) +
	  		geom_smooth(aes(x=year, y=rec, color=model),method="lm",size=1.15) +
	  		expand_limits(y=c(0)) +
	  		theme_few(base_size = 20) +
	  		scale_color_manual("Manual",values = c("black",viridis::viridis_pal(alpha = 1, begin = 0.1, end = 0.8, direction = 1, option = "H")(7)))
			ggsave(filename=paste0("rec.reg.png"), plot = p, device = "png", path = dir_plot,scale = 1, width = 9, height = 9, units = c("in"),dpi = 300, limitsize = TRUE)
	

