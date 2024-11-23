

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
	dir.create(dir_model,recursive=TRUE,showWarnings=FALSE)
	dir_model_source = "stem/"
	model_name_stem = "effort_creep/"
	dir_input_files = paste0(dir_model,dir_model_source,"/")
	dir_mfcl = paste0(proj_dir,"/MFCL/Windows/v2084/")

#_____________________________________________________________________________________________________________________________
# set up effort creep vector
	effort_creep.vec = c(0,1,1.5,2,2.5,3,6)

	for(i in 1:length(effort_creep.vec))
	{	
		#_____________________________________________________________________________________________________________________________
		# define directory structure
			model_name = effort_creep.vec[i]
			dir_run = paste0(dir_model,model_name_stem,model_name,"/") 
			dir.create(dir_run,recursive=TRUE,showWarnings=FALSE)
			dir_plot = paste0(proj_dir,"/Plot/Model_Runs/",model_name_stem,model_name,"/")
			dir.create(dir_plot,recursive=TRUE,showWarnings=FALSE)

		#_____________________________________________________________________________________________________________________________
		# transfer files
			# intial files
				FileList=c("labels.tmp","mfcl.cfg",paste0("PHASE0",1:7,".txt"),"skj.frq","skj.ini","skj.tag")
				file.copy(paste0(dir_input_files,FileList),dir_run,overwrite=TRUE)
			# mfcl files
				FileList=list.files(dir_mfcl)
				file.copy(paste0(dir_mfcl,FileList),dir_run,overwrite=TRUE)


		#_____________________________________________________________________________________________________________________________
		# modify input files
		 	dummy_frq = readfrq(paste0(dir_run,"skj.frq"))
		 	catch_dt = as.data.table(cateffpen(dummy_frq)) %>%
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
						tmp_scalar[s] = tmp_scalar[s-1]*(1+(effort_creep.vec[i]/400))
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
				cateffpen = catch_dt %>%
							.[order(fishery,year,month)] %>%
							.[,id:=NULL] %>%
							.[,ts:=NULL] %>%
							.[,e_scalar:=NULL]

			cateffpen(dummy_frq) = cateffpen
			writefrq(dummy_frq,paste0(dir_run,"skj.frq"))
  			p = plot.frqit(dummy_frq, Frq.names = NULL, fdesc = NULL, save.dir=dir_run, save.name="skj.frq") 

		#_____________________________________________________________________________________________________________________________
		# run MFCL
			setwd(dir_run)
			# build 00.par
			system("powershell ./mfclo64.exe skj.frq skj.ini 00.par -makepar")

			# run phases
			for(p in 1:7)
			{
				system(paste0("powershell ./mfclo64.exe skj.frq 0",p-1,".par 0",p,".par -file PHASE0",p,".txt"))
			}

		rm(list=c("catch_dt","cateffpen","dummy_frq","p"))
	}	

	for(i in 1:length(effort_creep.vec))
	{
		model_name = effort_creep.vec[i]
		dir_run = paste0(dir_model,model_name_stem,model_name,"/")
		dir_run_tmp =  paste0(dir_model,model_name_stem,model_name,"/tmp/")
		dir.create(dir_run_tmp,recursive=TRUE,showWarnings=FALSE)
		FileList=c("labels.tmp","mfcl.cfg","skj.frq","routine_save","skj.tag")
		file.copy(paste0(dir_run,FileList),dir_run_tmp,overwrite=TRUE)
		FileList=list.files(dir_mfcl)
		file.copy(paste0(dir_mfcl,FileList),dir_run_tmp,overwrite=TRUE)
		setwd(dir_run_tmp)
		system("powershell ./mfclo64.exe skj.frq routine_save tmp.par -switch 1 1 1 1")
	}
