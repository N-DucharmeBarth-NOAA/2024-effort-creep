

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
# extract quantities
    # the relative recruitment (overall rec dev)
	diag_relrec = as.data.table(rel_rec(par.list[[1]])) %>%
                  .[,ts:=as.numeric(year)+(as.numeric(season)-1)/4] %>%
                  .[,type:="relrec"] %>%
                  .[,.(type,ts,value)]
    # all catch in mt, exclude LL fleets (3,6,9,17,21,23,27,31)
    diag_catch = as.data.table(catch_pred(rep.list[[1]])) %>%
                 .[!(unit%in%c("3","6","9","17","21","23","27","31"))] %>%
                 .[,.(value=sum(value)),by=.(year,season)] %>%
                 .[,ts:=as.numeric(year)+(as.numeric(season)-1)/4] %>%
                  .[,type:="catch"] %>%
                  .[,.(type,ts,value)]
    # obs cpue
    	frq = readfrq(paste0(proj_dir,"/Model_Runs/Length100_GrowthDiag_Mix1_H0.8/skj.frq"))
		diag_obscpue = as.data.table(cateffpen(frq)) %>%
		 			   .[,id:=paste0(year,"_",month,"_",fishery)] %>%
					   .[,.(id,year,month,week,fishery,catch,effort,penalty)] %>%
					   .[,ts:=as.numeric(year)+(as.numeric(month)-2)/12] %>%
					   .[,cpue:=catch/effort] %>%
                       .[fishery %in% c(1,4,7,22,24,28,12,19)] %>%
                       .[penalty>0] %>%
                       .[,sd:=1/sqrt(2*penalty)] %>%
                       .[,mean_penalty:=mean(penalty),by=fishery] %>%
                       .[,mean_sd:=mean(sd),by=fishery]
               
    # ggplot(diag_obscpue) +
    # facet_wrap(~fleet) +
    # geom_path(aes(x=ts,y=value))


	