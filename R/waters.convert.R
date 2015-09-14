## Functions for converter ####################


no_scan_empty_warn <- function(w) if( any( grepl( "MS1 scans empty. Skipping profile matrix calculation.", w) ) ) 
  invokeRestart( "muffleWarning" )

no_split_error_warn<- function(w) if( any( grepl( "MSn information will be dropped", w) ) ) 
  invokeRestart( "muffleWarning" )


remakeTIC<-function(object){
  object@tic<-numeric(length(object@scanindex))
  for(i in 1:length(object@scanindex)){
    object@tic[i]<-sum(getScan(object, i)[,"intensity"])
  }
  return(object)
}



LockMassFun_rem<-function(object,idx){
  #idx<-AutoLockMass(object)  
  #idx=FindLockMass(object)
  #select=c(0,1)
  
  if(missing(idx)){stop('Parameters are missing')}
  
  select=c(1,1)
  
  #   if (any(diff(idx)==1)){
  #     select=c(0,1)
  #   }else{
  #     select=c(1,1)
  #   }
  
  
  
  idx.TF<-as.logical(rep(select, length(idx)/2))
  idx.TF<-which(idx.TF != 0)
  inx<-rep(TRUE, length(object@scanindex))
  inx[idx[idx.TF]]<-FALSE
  
  
  # separate scans into list
  scans=list()
  for (i in 1:length(object@scanindex)){
    scans[[i]]=getScan(object,i)  
  }
  
  
  # Just removed unwanted scans...
  scans=scans[inx]
  
  
  # Create new scanindex
  scanindex_new=numeric()
  for (i in 1:length(scans)){
    if (i==1){
      scanindex_new[i]=0
    }else{
      scanindex_new[i] = scanindex_new[i-1] + nrow(scans[[(i-1)]])
    }
  }
  
  
  # Put data into new object
  ob<-new("xcmsRaw")
  ob@env = new.env(parent=.GlobalEnv)
  
  scans=do.call(rbind,scans)
  ob@env$mz<-scans[,'mz']
  ob@env$intensity<-scans[,'intensity']
  
  ob@scantime <-object@scantime[inx]
  ob@scanindex <-as.integer(scanindex_new)
  ob@polarity   <- object@polarity[inx]
  ob@acquisitionNum <- 1:length(ob@scanindex)
  ob@profmethod     <- object@profmethod
  ob@mzrange        <- object@mzrange
  ob@filepath       <- object@filepath
  profStep(ob) =   profStep(object)
  ob=remakeTIC(ob)
  
  return(ob)
}







LockMassFun_rep<-function(object,idx){
  #idx<-AutoLockMass(object)
  #idx=FindLockMass(object)
  #select=c(0,1)
  
  
  if(missing(idx)){stop('Parameters are missing')}
  select=c(1,1)
  
  #   if (any(diff(idx)==1)){
  #     select=c(0,1)
  #   }else{
  #     select=c(1,1)
  #   }
  
  idx=idx[idx>1] # We don't allow to replace scan 1 since there is nothing before to replace with
  
  idx.TF<-as.logical(rep(select, length(idx)/2))
  idx.TF<-which(idx.TF != 0)
  inx<-rep(TRUE, length(object@scanindex))
  inx[idx[idx.TF]]<-FALSE
  
  
  # separate scans into list
  scans=list()
  for (i in 1:length(object@scanindex)){
    scans[[i]]=getScan(object,i)  
  }
  
  
  # Here we do the fillin instead
  for (i in idx[idx.TF]){
    scans[[i]]=scans[[i-1]]
  }
  
  
  # Create new scanindex
  scanindex_new=numeric()
  for (i in 1:length(scans)){
    if (i==1){
      scanindex_new[i]=0
    }else{
      scanindex_new[i] = scanindex_new[i-1] + nrow(scans[[(i-1)]])
    }
  }
  
  
  # Put data into new object
  ob<-new("xcmsRaw")
  ob@env = new.env(parent=.GlobalEnv)
  
  scans=do.call(rbind,scans)
  ob@env$mz<-scans[,'mz']
  ob@env$intensity<-scans[,'intensity']
  
  ob@scantime       <-object@scantime
  ob@scanindex      <-as.integer(scanindex_new)
  ob@polarity       <- object@polarity
  ob@acquisitionNum <- object@acquisitionNum
  ob@profmethod     <- object@profmethod
  ob@mzrange        <- object@mzrange
  ob@filepath       <- object@filepath
  profStep(ob)      <- profStep(object)
  ob=remakeTIC(ob)
  
  return(ob)
}



MSnSplit<-function(object,f){
  ob = withCallingHandlers( split(msn2xcms(object),  f   )$'TRUE' , warning = no_split_error_warn )
  
  ob@env$msnIntensity    =    ob@env$intensity
  ob@env$msnMz           =    ob@env$mz
  
  ob@tic                      = numeric()
  ob@msnScanindex             = ob@scanindex
  ob@scanindex                = integer()
  ob@scantime                 = numeric()
  ob@acquisitionNum           = integer()
  ob@mzrange                  = numeric()
  ob@msnAcquisitionNum        = ob@msnAcquisitionNum[f]
  ob@msnPrecursorScan         = ob@msnPrecursorScan[f]
  ob@msnLevel                 = ob@msnLevel[f]
  ob@msnRt                    = ob@msnRt[f]
  ob@msnPrecursorMz           = ob@msnPrecursorMz[f]
  ob@msnPrecursorIntensity    = ob@msnPrecursorIntensity[f]
  ob@msnPrecursorCharge       = ob@msnPrecursorCharge[f]
  ob@msnCollisionEnergy       = ob@msnCollisionEnergy[f]
  
  return(ob)
}




MSnExtractLockmass<-function(object,MSparam=NULL){
  if (missing(MSparam)){stop('MSparam has to be set.')}
  
  lockmass=grep('Lock Mass',MSparam)
  
  if (!(length(lockmass)==0)){
    lockmass=MSparam[lockmass[length(lockmass)-1]]
    lockmass = strsplit(lockmass,"\t")
    lockmass = as.numeric(lockmass[[1]][length(lockmass[[1]])])
    
  }else{
    lockmass=grep('EDC Mass',MSparam)
    lockmass=MSparam[lockmass]
    lockmass = strsplit(lockmass,"\t")
    lockmass = as.numeric(lockmass[[1]][length(lockmass[[1]])])
  }
  
  
  lockhist_freq = as.matrix(table(object@msnPrecursorMz))
  lockhist_masses = as.numeric(rownames(lockhist_freq))
  lockhist_masses = lockhist_masses[   which.min(abs(lockhist_masses-lockmass))   ]
  
  return(lockhist_masses)
}











convert.waters=function(indir,outdir){
  
  ## Get file list
  files = list.files(indir, pattern='.raw', full.names=T)
  
  namelist=character()
  filenames=character()
  
  log_convert = list()
  log_warnings=character()
  log_convert_warning = list()
  log_convert_warning[[1]]=' '
  
  
  
  dir.create(outdir, showWarnings = FALSE,recursive=T)
  
  
  
  for (i in 1:length(files)) {
    
    
    ## Check if lockspray is present. TODO: Scream at waters
    input_method = scan(paste(files[i],'/_extern.inf',sep=""), character(0), sep = "\n",quiet=T)
    lockspray_present = any(grepl('Lock Spray Configuration',input_method))
    
    if (!lockspray_present){
      warning(paste('LockMass seems to be disabled for file ',basename(files[i]),'!', ' Saving data as is.',sep=""),immediate. = T)
    }
    
    ## Get name of sample
    name = scan(paste(files[i],'/_HEADER.TXT',sep=""), character(0), sep = "\n",quiet=T)
    name_where = grepl('Sample Description',name)
    namelist[i] = strsplit(name[name_where],': ')[[1]][2]
    rm(name,name_where)
    
    
    ## Get polarity
    ionmode = scan(paste(files[i],'/_extern.inf',sep=""), character(0), sep = "\n",quiet=T)
    
    
    if (!(any(grepl('ES+',ionmode,fixed = T)) | any(grepl('ES-',ionmode,fixed = T)))){warning('Ionization mode cannot be detected!')
      
    }else{
      if(any(grepl('ES+',ionmode,fixed = T))){ionmode='positive'
      }else{ionmode='negative'}
    }
    
    
    
    
    
    
    ## Get the functions and figure if there is MS1, MSn and/or MSE
    func_start = grep('Function Parameters - Function',input_method,fixed=T)
    func_end   = c(   grep('Function Parameters',input_method,fixed=T)-1   ,   length(input_method) )
    func_end = func_end[   2:length(func_end)   ]
    
    
    MSE=F
    MSn=F
    MS1=F
    
    functions=list()
    for (i2 in 1:(length(func_start))){
      functions[[i2]]=input_method[    func_start[i2]:func_end[i2]    ]  
    }
    
    
    not_ref = unlist(lapply(functions, function(x) !grepl(' - REFERENCE',x[1])))
    MSn = any(unlist(lapply(functions, function(x) grepl('TOF MSMS FUNCTION',x[1],fixed=T))))
    
    # MS1 = any(unlist(lapply(functions, function(x) grepl('TOF MS FUNCTION',x[1],fixed=T))))
    #  MS1 = any(unlist(lapply(functions, function(x) grepl('Trap MS Collision Energy (eV)\t\t\t4.0',x,fixed=T))))
    if (!MSn){MS1=T}
    
    
    if (!MSn){
      for (i2 in which(not_ref)){
        if (any(grepl('Collision Energy Ramp Start',functions[[i2]],fixed=T))){
          MSE=T
        }
      }
    }
    
    
    
    
    
    
    ## Convert the data
    tempfile = tempfile(pattern = basename(gsub('.raw','',files[i])), tmpdir = '/', fileext = ".mzXML")
    
    
    # if there is mobility data we need to make a copy and remove it. Otherwise masswolf doesn't work
    if (file.exists(   paste(  files[i],'/apex3D',sep=''))){
      input_file=  paste(  gsub('.raw','',files[i],fixed=T),'_temp.raw',sep='')
      dir.create(input_file)
      to_copy=list.files(path = files[i], all.files = T,full.names = T, recursive = F,ignore.case = T, include.dirs = T,no.. = T)
      to_copy=to_copy[      !(grepl('/_mob',to_copy,fixed=T) | grepl('/apex3D',to_copy,fixed=T) | grepl('.cdt',to_copy,fixed=T) | grepl('.ind',to_copy,fixed=T))   ]
      temp=file.copy(to_copy, input_file,overwrite=F,recursive=T)
      rm(temp)
    }else{
      input_file= files[i]
    }
    
    
    # Convert with mass wolf
    if (MSE){
      log_convert[[i]]=system(paste('masswolf --MSe --mzXML ',"\"",input_file,"\"",' ',"\"",outdir,tempfile,"\"",sep=""), intern=T)
    }else{
      log_convert[[i]]=system(paste('masswolf --mzXML ',"\"",input_file,"\"",' ',"\"",outdir,tempfile,"\"",sep=""), intern=T)  
    }
    
    # Delete temporary file of created (for files that include mobility data)
    if (file.exists(   paste(  files[i],'/apex3D',sep=''))){
      unlink(input_file,recursive=T)
    }
    
    
    
    # Read the raw Waters data
    raw_file_data = scan(paste(outdir,tempfile,sep=''), character(0), sep = "\n",quiet=T)
    
    scan_start = grep('<scan num=',raw_file_data,fixed=T)
    scan_end = c(    scan_start[2:length(scan_start)]-1    ,  length(raw_file_data)  )
    
    
    raw_file_data_scans=list()
    for (i2 in 1:(length(scan_start))){
      raw_file_data_scans[[i2]]=raw_file_data[    scan_start[i2]:scan_end[i2]    ]
      
    }
    
    
    
    
    
    ## Read retention times and scan numbers in original data
    
    # Get scan retention times
    raw_file_data_ret = lapply(raw_file_data_scans, function(x) as.numeric(str_trim(str_replace(str_replace(x[     grep('retentionTime',x)   ],'retentionTime=\"PT', ''),'S\"',''),side='both'))   )
    raw_file_data_ret = unlist(raw_file_data_ret)
    
    # Get which are lockmass (calibration) scans
    is_lockmass_scan =   lapply(raw_file_data_scans, function(x)   any(grepl('scanType=\"calibration\"',x,fixed=T))   )
    is_lockmass_scan = unlist(is_lockmass_scan)
    
    
    # Get which are MSE or MS2 scans
    #is_MSE_or_MSn_scan =   lapply(raw_file_data_scans, function(x)   any(grepl('msLevel=\"2\"',x,fixed=T))   )
    #is_MSE_or_MSn_scan = unlist(is_MSE_scan)
    
    
    
    
    
    ## Read file
    xr_org = withCallingHandlers( xcmsRaw(paste(outdir,tempfile,sep=""),includeMSn=T), warning = no_scan_empty_warn )
    
    
    
    
    ## if it is an MSE file, write seperate MSE file
    if (MSE){
      xr_MSE =new('xcmsRaw')
      
      #xr_MSE@env = xr_org@env
      copyEnv(xr_org@env,xr_MSE@env)
      xr_MSE@msnScanindex = xr_org@msnScanindex
      xr_MSE@msnAcquisitionNum = xr_org@msnAcquisitionNum
      xr_MSE@msnPrecursorScan = xr_org@msnPrecursorScan
      xr_MSE@msnLevel = xr_org@msnLevel
      xr_MSE@msnRt = xr_org@msnRt
      
      xr_MSE@msnPrecursorMz = xr_org@msnPrecursorMz
      xr_MSE@msnPrecursorIntensity = xr_org@msnPrecursorIntensity
      xr_MSE@msnPrecursorCharge = xr_org@msnPrecursorCharge
      xr_MSE@msnCollisionEnergy = xr_org@msnCollisionEnergy
      xr_MSE@filepath = xr_org@filepath
      
      
      filename=basename(gsub('.raw','_MSE.mzData',files[i]))
      xr_MSE=msn2xcms(xr_MSE)
      write.mzdata(xr_MSE,paste(outdir,'/',filename,sep=""))
    }
    
    
    
    
    ## Write MS1 data
    if (MS1){
      
      
      if (lockspray_present){
        
        idx=match(raw_file_data_ret[is_lockmass_scan],xr_org@scantime)
        idx = idx[!is.na(idx)]
        
        xr_fixed=LockMassFun_rep(xr_org,idx)
        xr_fixed@polarity=factor(rep(ionmode,1,length(xr_fixed@polarity)),levels=c('negative','positive','unknown'))
        
      }else{
        xr_fixed=xr_org
        xr_fixed@polarity=factor(rep(ionmode,1,length(xr_fixed@polarity)),levels=c('negative','positive','unknown'))
      }
      
      
      
      
      # Write fixed file
      filenames[i]=basename(gsub('.raw','.mzData',files[i]))
      write.mzdata(xr_fixed,paste(outdir,'/',filenames[i],sep=""))
      
      temp1 = plotChrom(xr_org)
      temp2 = plotChrom(xr_fixed)
    }
    
    
    
    
    
    
    
    ## Write MSn data
    
    if (MSn){
      
      
      if (lockspray_present){
        
        LockMass = MSnExtractLockmass(xr_org,MSparam=input_method)
        f = xr_org@msnPrecursorMz != LockMass
        
        if(any(f)){
          log_convert_warning[length(log_convert_warning)+1]=print(paste('Lockmass detected at ',formatC(LockMass,digits=7),'Da in range ',format(min(xr_org@msnRt[f])/60,digits=2),'min to ',format(max(xr_org@msnRt[f])/60,digits=2),'min for file: ',basename(files[i]),sep=""))
          xr_fixed = MSnSplit(xr_org,f)
          xr_fixed@polarity=factor(rep(ionmode,1,length(xr_fixed@polarity)),levels=c('negative','positive','unknown'))
        }else{
          log_convert_warning[length(log_convert_warning)+1]=print(paste('Lockmass was enabled but no lockmass scans could be found. Data written as is for file: ',basename(files[i]),sep=""))
          xr_fixed=xr_org
        }
        
        
      }else{
        xr_fixed=xr_org
        xr_fixed@polarity=factor(rep(ionmode,1,length(xr_fixed@polarity)),levels=c('negative','positive','unknown'))
      }
      
      
      # Write fixed file
      filenames[i]=basename(gsub('.raw','.mzData',files[i]))
      write.mzdata(xr_fixed,paste(outdir,'/',filenames[i],sep=""))
      
      temp1=msn2xcms(xr_org)
      temp1=remakeTIC(temp1)    
      temp1=plotTIC(temp1)
      temp2=msn2xcms(xr_fixed)
      temp2=remakeTIC(temp2)    
      temp2=plotTIC(temp2)
    }
    
    
    
    
    ## Debugging plots
    ## debugging, plotting
    dir.create(paste(outdir,'/debug/',sep=""),showWarnings=F)
    png(paste(outdir,'/debug/',basename(gsub('.raw','.png',files[i])),sep=""),width=2000,height=1500,res=300)
    plot(temp1,type='l',col='red',xlim=c(-10,max(temp1[,1])))
    lines(temp2,col='green')
    dev.off()
    rm(temp1,temp2)
    
    
    
    # delete temp file
    unlink(paste(outdir,tempfile,sep=""))
    rm(xr_org,xr_fixed,xr_MSE)
    gc()
    log_warnings[i]=warnings()
  }
  
  
  
  
  write.csv(cbind(filenames,namelist),paste(outdir,'/','filelist.csv',sep=""))
  
  
  fileConn<-file(paste(outdir,'/','log_convert.log',sep=""))
  writeLines(unlist(log_convert), fileConn)
  close(fileConn)
  
  fileConn<-file(paste(outdir,'/','log_warnings.log',sep=""))
  writeLines(unlist(as.character(log_warnings)), fileConn)
  close(fileConn)
  
  fileConn<-file(paste(outdir,'/','log_convert_warning.log',sep=""))
  writeLines(unlist(log_convert_warning), fileConn)
  close(fileConn)
  
  
  
  
}










convert.waters2=function(infiles,outdir,funcs=c(1) ){
  
  
  
  for(i in 1:length(infiles)){
    
    # Figure out which function number each file belongs to
    raw_files <- list.files(infiles[i])
    m <- regexec("^_FUNC00(.*)(\\.)",raw_files)
    m <- regmatches(raw_files, m)
    raw_files <- cbind.data.frame(raw_files = raw_files, func_nr = as.numeric(sapply(m,function(x) x[2])),stringsAsFactors=FALSE)
    
    
    # Do conversion for each function
    for(funcs_sel in funcs){
      
      # make temp .raw folder
      temp_raw <-paste0(outdir,"/temp.raw") 
      dir.create(temp_raw, showWarnings = FALSE,recursive=T)
      
      # copy the required files to a temp dir
      to_copy <- raw_files$raw_files[raw_files$func_nr==funcs_sel | is.na(raw_files$func_nr)]
      to_copy_dest <- sub("00[0-9]","001",to_copy)
      file.symlink(from = paste0(infiles[i],"/",to_copy), to = paste0(temp_raw,"/",to_copy_dest))
      
      # Do the conversion
      system(paste0('masswolf --mzXML ',"\"",temp_raw,"\"",' ',"\"",outdir,"/",basename(sub("\\.[^.]*$", "", infiles[i]) ),"_func",funcs_sel,".mzXML\""), intern=T)  
      unlink(temp_raw,recursive = TRUE)
      
    }
    
  }
  
}