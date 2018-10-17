library(shiny);  

server <- function(input, output){
  
  datannot <- reactiveValues(annot = NULL, outList = NULL, type_ = NULL); 
  result_ <- reactiveValues(result_ = NULL, draw_ = NULL)
  shinyjs::runjs('$("#barplot").hide();$("#col1id").hide();');
  
  observeEvent(input$annotfile,{
    #browser();
    tryCatch({
      
      # file extension
      ext <- toupper(tools::file_ext(input$annotfile$name)); annot = NULL;
      
      if(!(ext %in% c("TXT","CSV","XLSX"))){
        datannot$annot <- NULL; annot = as.data.frame(c("Error"), stringsAsFactors = F); colnames(annot) <- c("Error");
        createAlert(session, "alertannotfile", "alert1", title = "Error", content = "Only .csv, .txt and .xlsx are supported", append = F)
      }
      else{
          if (ext == 'XLSX') annot <- openxlsx::read.xlsx(input$annotfile$datapath)
          else if (ext %in% c("TXT", "CSV")) annot <- read.table(file = input$annotfile$datapath, header = T, sep=",", row.names = NULL, fill = T)
          
          # take care of NA's and empty rows/cols       
          annot <- data.frame(lapply(annot, as.character), stringsAsFactors=F)
          annot <- annot[!apply(is.na(annot) | annot == "", 1, all),] # rows with all NA
          annot <- annot[,!apply(is.na(annot) | annot == "", 2, all)] # cols with all NA
          
          # check for missing columns
          outList.names <- c("PairIndex", "Response", "Drug1", "Drug2", "Conc1", "Conc2", "ConcUnit", "Type")
          mismatch = outList.names[!(colnames(annot) %in% outList.names)];
          
          if(length(mismatch) == 0){
            dataOutput <- data.frame(
              Response = round(as.numeric(annot$Response),2), Drug1 = as.character(annot$Drug1), 
              Drug2 = as.character(annot$Drug2), Conc1 = round(as.numeric(annot$Conc1),2), Conc2 = round(as.numeric(annot$Conc2),2), 
              ConcUnit = as.character(annot$ConcUnit), Type = as.character(annot$Type), stringsAsFactors=F)
            datannot$annot <- dataOutput;
          }
          else{
            mismatch = na.omit(mismatch);
            createAlert(session, "noPDdata", "alertPD", title = "Error", content = paste0("Annotation table doesnot contain " ,paste0(mismatch, collapse = ", "), " column(s). \n See example data."), append = F)
            datannot$annot <- NULL
          }
      }
    }, error = function(e) {
      toastr_error(paste0("Something wrong with your input file that cannot be handled by application"), title = "Unhandled error occurred!", closeButton = !0, progressBar = !0, position = "top-right", preventDuplicates = !0,
                   showDuration = 300, hideDuration = 1000, timeOut = 10000, extendedTimeOut = 1000, showEasing = "swing",
                   hideEasing = "swing", showMethod = "fadeIn", hideMethod = "fadeOut")
    })
    
  })
  
  output$video <- renderUI({
    if(!is.null(input$GetScreenWidth)){
      width = round(input$GetScreenWidth * 0.63); height = round(input$GetScreenHeight * 0.57);
      HTML(paste0('<iframe id="ifr1" src="https://player.vimeo.com/video/290934199" width="',width,'" height="',height,'" frameborder="0" webkitallowfullscreen mozallowfullscreen allowfullscreen></iframe>'))
    }
  })
  
  # once annot file is loaded
  observeEvent(input$analysis_, {
    #browser();
    #browser();
    if(!is.null(datannot$annot)){
      runjs("$('#video').remove(); $('#barplot').show(); $('#col1id').show();");

      withProgress(message = 'Calculations in process', value = 0, {
        library(parallel); library(synergyfinder); library(kriging)
        
  
        #browser()
        #input$conrol_
        # load sample-control
        sample_ = datannot$annot[datannot$annot$Type=="sample", ]; 
  
        # unique combinations
        combis_ = unique(paste0(sample_$Drug1,"|",sample_$Drug2));
        sample_$PairIndex = match(paste0(sample_$Drug1,"|",sample_$Drug2), combis_);
        incProgress(1/3, detail = paste0("Quantifying synergy", paste0(" (",input$state," method)")))
        d_ = do.call("rbind", lapply(1:length(unique(sample_$PairIndex)), function(i) sample_[sample_$PairIndex == i,][1,]))
        
#browser();
        # calculate synergy for all combinations using synergyfinder R package
        
        
        scores <<- synergyfinder::CalculateSynergy(transformInputData(sample_, "inhibition"),input$state)
        synScores <- round(sapply(1:length(combis_), function(i){
          p = PlotSynMatr(scores$scores[[i]], "inhibition", d1N = scores$drug.pairs[i,1], d2N = scores$drug.pairs[i,2]);
          c_ = as.numeric(colnames(scores$scores[[i]])); r_ = as.numeric(rownames(scores$scores[[i]])); Matr = scores$scores[[i]];
          filename_ = file.path(paste0("./plots/",scores$drug.pairs[i,1],"_syn_",scores$drug.pairs[i,2],".png"))
          png(filename = filename_,width=length(c_)*40,height=length(r_)*40, bg = "transparent")
          print(p); dev.off()
          #Vol = sum(unlist(sapply(2:nrow(Matr), function(k){j = 2:(ncol(Matr));v_=((c_[k]-c_[k-1]) * (r_[j]-r_[j-1])) * Matr[j,k]; 
          #if(k==nrow(Matr)) {v_ = v_[-length(v_)]}; v_}))); Vol / ((max(c_) - min(c_)) * (max(r_) - min(r_)))
          # Matr[nrow(Matr), ncol(Matr)] = NA; 
          # mean(Matr[-1,-1], na.rm = T)
          Vol = sum(unlist(sapply(3:nrow(Matr), function(k){j = 3:(ncol(Matr));v_=(log(c_[k]/c_[k-1]) * log(r_[j]/r_[j-1])) * Matr[j,k]; 
          if(k==nrow(Matr)) {v_ = v_[-length(v_)]}; v_}))); Vol / (log(max(c_) / min(c_[-1])) * log(max(r_)/min(r_[-1])))
        }), 3)
        incProgress(1/3, detail = "Quantifying potency")
  
        # calculate normalized volume under each inhibition matrix
        inhScores <- round(unlist(mclapply(1:length(combis_), function(i){
          Matr = scores$dose.response.mats[[i]]; c_ = as.numeric(colnames(Matr)); r_ = as.numeric(rownames(Matr))
          p = PlotRespMatr(Matr, "inhibition", d1N = scores$drug.pairs[i,1], d2N = scores$drug.pairs[i,2])
          filename_ = file.path(paste0("./plots/",scores$drug.pairs[i,1],"_inh_",scores$drug.pairs[i,2],".png"))
          png(filename = filename_,width=length(c_)*40,height=length(r_)*40, bg = "transparent")
          print(p); dev.off()
          # Vol = sum(unlist(sapply(2:nrow(Matr), function(k){j = 2:(ncol(Matr));v_=((c_[k]-c_[k-1]) * (r_[j]-r_[j-1])) * Matr[j,k]; 
          # if(k==nrow(Matr)) {v_ = v_[-length(v_)]}; v_}))); Vol / ((max(c_) - min(c_)) * (max(r_) - min(r_)))
          Vol = sum(unlist(sapply(3:nrow(Matr), function(k){j = 3:(ncol(Matr));v_=(log(c_[k]/c_[k-1]) * log(r_[j]/r_[j-1])) * Matr[j,k]; 
          if(k==nrow(Matr)) {v_ = v_[-length(v_)]}; v_}))); Vol / (log(max(c_) / min(c_[-1])) * log(max(r_)/min(r_[-1])))
        }, mc.cores = 1)), 3)
        incProgress(1/3, detail = "Quantifying toxicity")
  
        # calculate normalized volume under each toxicity matrix
        if(input$conrol_ == "Control"){
          incProgress(0, detail = "Quantifying toxicity from control data")
          control_ = datannot$annot[datannot$annot$Type=="control", ]
          control_$PairIndex = match(paste0(control_$Drug1,"|",control_$Drug2), combis_);
          dt_ <- transformInputData(control_, "toxicity");
          toxScores <- round(sapply(1:length(combis_), function(i){
            Matr = dt_$dose.response.mats[[i]]; c_ = as.numeric(colnames(Matr)); r_ = as.numeric(rownames(Matr))
            p = PlotRespMatr(Matr, "toxicity", d1N = scores$drug.pairs[i,1], d2N = scores$drug.pairs[i,2])
            filename_ = file.path(paste0("./plots/",scores$drug.pairs[i,1],"_tox_",scores$drug.pairs[i,2],".png"))
            png(filename = filename_,width=length(c_)*40,height=length(r_)*40, bg = "transparent")
            print(p); dev.off()
            # Vol = sum(unlist(sapply(2:nrow(Matr), function(k){j = 2:(ncol(Matr));v_=((c_[k]-c_[k-1]) * (r_[j]-r_[j-1])) * Matr[j,k]; 
            # if(k==nrow(Matr)) {v_ = v_[-length(v_)]}; v_}))); Vol / ((max(c_) - min(c_)) * (max(r_) - min(r_)))
            Vol = sum(unlist(sapply(3:nrow(Matr), function(k){j = 3:(ncol(Matr));v_=(log(c_[k]/c_[k-1]) * log(r_[j]/r_[j-1])) * Matr[j,k]; 
            if(k==nrow(Matr)) {v_ = v_[-length(v_)]}; v_}))); Vol / (log(max(c_) / min(c_[-1])) * log(max(r_)/min(r_[-1])))
          }), 3)
        } else if(input$conrol_ == "SynTox"){
          incProgress(0, detail = "Quantifying toxicity based on SynTox equation")
          toxScores <- round(sapply(1:length(combis_), function(i){ 
            Matr = scores$dose.response.mats[[i]]; c_ = as.numeric(colnames(Matr)); r_ = as.numeric(rownames(Matr))
            drug1 = CALC_IC50_EC50_DSS(Matr[-1,1],2);  drug2 = CALC_IC50_EC50_DSS(Matr[1,-1],2);
            yi1 = fitted(drug1$nls_result_ic50); drug1 = drug1$coef_ic50; yi2 = fitted(drug2$nls_result_ic50); drug2 = drug2$coef_ic50
            
            MatrN = round(cbind(Matr[,1]/100,rbind(Matr[1,-1]/100, 1 - sapply(2:nrow(Matr), function(k){j = 2:(ncol(Matr)); 
            ((100-yi1[j-1])/(100-drug1[["MIN"]]))*((100-yi2[k-1])/(100-drug2[["MIN"]]))
            }))),3)
            colnames(MatrN)[1]=colnames(Matr)[1]
            
            # p = PlotRespMatr(MatrN*100, "toxicity", d1N = scores$drug.pairs[i,1], d2N = scores$drug.pairs[i,2])
            # filename_ = file.path(paste0("./plots/",scores$drug.pairs[i,1],"_tox_appox_",scores$drug.pairs[i,2],".png"))
            # png(filename = filename_,width=length(c_)*40,height=length(r_)*40, bg = "transparent")
            # print(p); dev.off()
            
            Vol = sum(unlist(sapply(3:nrow(MatrN), function(k){j = 3:(ncol(MatrN));v_=(log(c_[k]/c_[k-1]) * log(r_[j]/r_[j-1])) * MatrN[j,k]; 
            if(k==nrow(MatrN)) {v_ = v_[-length(v_)]}; v_}))); Vol / (log(max(c_) / min(c_[-1])) * log(max(r_)/min(r_[-1])))
            
          }),3)
  
          # toxScores <- round(sapply(1:20, function(z){ 
          #   Matr = scores$dose.response.mats[[z]]; c_ = as.numeric(colnames(Matr)); r_ = as.numeric(rownames(Matr))
          #   drug1 = CALC_IC50_EC50_DSS(Matr[-1,1],2);  drug2 = CALC_IC50_EC50_DSS(Matr[1,-1],2); WeiSynVol = 0;
          #   yi1 = fitted(drug1$nls_result_ic50); drug1 = drug1$coef_ic50; yi2 = fitted(drug2$nls_result_ic50); drug2 = drug2$coef_ic50
          #   
          #   Vol = sum(unlist(sapply(2:nrow(Matr), function(k){j = 2:(ncol(Matr)); 
          #   Root_=sqrt(((100-yi1[j-1])/(100-drug1[["MIN"]]))*((100-yi2[k-1])/(100-drug2[["MIN"]])));
          #   v_=((c_[k]-c_[k-1]) * (r_[j]-r_[j-1])) * Root_; if(k==nrow(Matr)) {v_ = v_[-length(v_)]}; v_}))); 
          #   Vol / ((max(c_) - min(c_)) * (max(r_) - min(r_)))
          # }), 3)
        } else {
          
          incProgress(0, detail = "Quantifying toxicity (exporting compounds' CIDs from PubChem)")
          
          Dr_Info <- data.frame(compounds_ = unique(c(datannot$annot$Drug1, datannot$annot$Drug2)))
          Dr_Info$CID = mclapply(Dr_Info$compounds_, function(i){
            tryCatch({
              jsonlite::fromJSON(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",gsub(" ","",i),"/cids/JSON"))[[1]][[1]][[1]]
            }, error = function(e) {
              ""
            })
          }, mc.cores = 1)

          missingComp <- mclapply(which(Dr_Info$CID == ""), function(ind){
            tryCatch({
              jsonlite::fromJSON(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/name/",gsub(" ","",Dr_Info[ind,"compounds_"]),"/cids/JSON"))[[1]][[1]]$CID[[1]]
            }, error = function(e) {
              ""
            })
          }, mc.cores = 1)
          Dr_Info$CID[Dr_Info$CID == ""] = missingComp
          
          if(any(Dr_Info$CID == "")) {
     
            toastr_error(title = "Error!", message = paste0("Some compounds were not found in PubChem. <br> Please select another method for toxocity approximation. <br> Missing compounds are: ", paste0(Dr_Info$compounds_[Dr_Info$CID == ""], collapse = ", "), "."), closeButton = !0, progressBar = !0, position = "bottom-left", preventDuplicates = !0,
                           showDuration = 300, hideDuration = 1000, timeOut = 6000, extendedTimeOut = 1000, showEasing = "swing",
                           hideEasing = "swing", showMethod = "fadeIn", hideMethod = "fadeOut")
            
            
            toxScores <- NULL;
          } else {
            
            incProgress(0, detail = "Quantifying toxicity (exporting SMILES from PubChem)")
            Dr_Info$SMILES = mclapply(Dr_Info$CID, function(i){
              tryCatch({
                tmp_ = jsonlite::fromJSON(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",gsub(" ","",i),"/JSON"))$PC_Compounds$props[[1]]
                tmp_$value[tmp_$urn$label == "SMILES" & tmp_$urn$name == "Canonical","sval"]
              }, error = function(e) {
                ""
              })
            }, mc.cores = 1)
            
            library(ChemmineOB); 
            incProgress(0, detail = "Quantifying toxicity (predicting LD50 values for each drug pair)")
            FP4fing <- do.call("rbind", mclapply(Dr_Info$SMILES, function(sm) fingerprint_OB(forEachMol("SMILES",sm[[1]],identity), "FP4")))
            FP4fing = as.data.frame(FP4fing); rownames(FP4fing) = Dr_Info$compounds_
            # saveRDS(FP4fing,"FP4fing.RDS")
            
            library(randomForest);
            rf_models <- list(readRDS("./models/fit1RFfing.RDS"),readRDS("./models/fit2RFfing.RDS"),readRDS("./models/fit3RFfing.RDS"),
                              readRDS("./models/fit4RFfing.RDS"),readRDS("./models/fit5RFfing.RDS"))
            tox <- rowMeans(do.call("cbind", lapply(1:5, function(i) predict(rf_models[[i]], FP4fing))));
            toxScores <- round((as.numeric(tox[match(d_$Drug1,names(tox))]) + as.numeric(tox[match(d_$Drug2,names(tox))]))/2, 3)
          }
          
        }

        # results 
        results_ <- data.frame(synScores=(as.numeric(factor(synScores))) / length(unique(as.numeric(factor(synScores)))), inhScores=(as.numeric(factor(inhScores))) / length(unique(as.numeric(factor(inhScores)))),
                               toxScores=1 - (as.numeric(factor(toxScores))) / length(unique(as.numeric(factor(toxScores)))),
                               synScoresReal = synScores, inhScoresReal = inhScores, toxScoresReal = toxScores)
        results_$score = (results_$synScores + results_$inhScores + results_$toxScores) / 3
        results_$scoreIT = (results_$inhScores + results_$toxScores) / 2
        results_$scoreIS = (results_$synScores + results_$inhScores) / 2
        results_$scoreST = (results_$synScores + results_$toxScores) / 2
        #browser();
        #browser();
        
        # add drugs
                results_$Drug1 <- d_$Drug1; results_$Drug2 <- d_$Drug2;
        #browser();
        #browser();
        
        # get figures
        Inh_base64 <<- unlist(lapply(1:length(combis_), function(i){
         gsub("\r?\n|\r", " ", base64::img(paste0("./plots/",scores$drug.pairs[i,1],"_inh_",scores$drug.pairs[i,2],".png"))) 
        }))
        if(input$conrol_ == "Control"){
          Tox_base64 <<- unlist(lapply(1:length(combis_), function(i){
            gsub("\r?\n|\r", " ", base64::img(paste0("./plots/",scores$drug.pairs[i,1],"_tox_",scores$drug.pairs[i,2],".png"))) 
          }))
        } else {Tox_base64 <<- NULL}
        Syn_base64 <<- unlist(lapply(1:length(combis_), function(i){
          gsub("\r?\n|\r", " ", base64::img(paste0("./plots/",scores$drug.pairs[i,1],"_syn_",scores$drug.pairs[i,2],".png"))) 
        }))
        #browser();    browser();    browser();
        result_$result_ <<- results_;
        result_$result_$barscore <<- result_$result_$score
        #results_ <<- results_
        result_$draw_ <<- !0;
      })
    } else {
      toastr_warning(title = "Warning!", message = "Please upload the data first", closeButton = !0, progressBar = !0, position = "bottom-left", preventDuplicates = !0,
                     showDuration = 300, hideDuration = 1000, timeOut = 6000, extendedTimeOut = 1000, showEasing = "swing",
                     hideEasing = "swing", showMethod = "fadeIn", hideMethod = "fadeOut")
    }
    
    #save.image("image1test.RDATA")
    
  })
  
 output$downlEx <- downloadHandler(
    filename = 'Example_data.xlsx',
    content = function(file) {
      file.copy("./AnnotExample.xlsx",file)
    },
    contentType = NULL
  )
 
 output$exp_excel <- downloadHandler(
   filename = 'Summary_table.xlsx',
   content = function(file) {
     #browser();browser();
     file_ = result_$result_;
     file_ = file_[ ,c("Drug1","Drug2","score","synScoresReal","inhScoresReal","toxScoresReal","synScores","inhScores","toxScores")]
     colnames(file_) = c("Drug1","Drug2","STP score", paste0("Synergy score (",input$state," model)"),"Potency score",paste0("Toxicity score (",
                         names(which(list("From control data" = "Control", "SynTox Equation" = "SynTox", "QSTR based approximation" = "QSTR") =="QSTR")),")"), "Synergy rank", "Potency rank", "Toxicity rank")
     
     openxlsx::write.xlsx(file_, "results_.xlsx", asTable = T);
     file.copy("./results_.xlsx",file)
   },
   contentType = NULL
 )
 
 
    
  output$exp_2d2 <- downloadHandler(
    filename = function() {
      paste('figure_barplot_2D-', Sys.Date(), '.pdf', sep='')
    },
    content = function(file) {
      ggsave("fig2D.pdf", plot2Dbar(), width = 20, height = 20, units = "cm")
      file.copy("./fig2D.pdf",file)
    },
    contentType = NULL
  )
  
  output$exp_2d <- downloadHandler(
    filename = function() {
      paste('figure_2D-', Sys.Date(), '.pdf', sep='')
    },
    content = function(file) {
      #browser();
      ggsave("fig2D.pdf", plot2Dgr(), width = 20, height = 20, units = "cm")
      file.copy("./fig2D.pdf",file)
    },
    contentType = NULL
  )
  
  # results_$vv = sapply(1:20, function(i) mean(results_$ttt.x[i], results_$ttt.y[i]))
  # ggplot(results_, aes(x=toxScoresReal, y=vv)) +
  #   geom_point(size=2) + 
  #   geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + theme_classic()
  # 
  # # calculate synergy for all combinations using synergyfinder R package
  # sample_$PairIndex = match(paste0(sample_$Drug1,"|",sample_$Drug2), combis_);
  # dt_ <- transformInputData(sample_, "inhibition");
  # scores <<- synergyfinder::CalculateSynergy(dt_,"Bliss")
  # 
  # vv = c();
  # for(z in 1:20){
  # 
  #   InteractionCp=Interaction=scores$scores[[z]][-1,-1]; Matr = scores$dose.response.mats[[z]]
  # 
  #   # # drugs
  #   drug1 = CALC_IC50_EC50_DSS(Matr[-1,1],2);  drug2 = CALC_IC50_EC50_DSS(Matr[1,-1],2);
  #   d1Model = drug1$nls_result_ic50; drug1 = drug1$coef_ic50; d2Model = drug2$nls_result_ic50; drug2 = drug2$coef_ic50
  #   #
  #   yi1 = fitted(d1Model)
  #   yi2 = fitted(d2Model)
  # 
  #   #Matr = 100 - Matr;
  #   WeiSynVol = 0;
  #   for(i in 1:7){
  #     for(j in 1:7){
  #       Root_ = sqrt(((100 - yi1[i]) / (100 - drug1[["MIN"]])) * ((100 - yi2[j]) / (100 - drug2[["MIN"]])));
  #       Incr = (as.numeric(colnames(Matr))[i+1]-as.numeric(colnames(Matr))[i]) * (as.numeric(rownames(Matr))[j+1]-as.numeric(rownames(Matr))[j])
  #       WeiSynVol = WeiSynVol +  Incr * Root_
  #     }
  #   }
  #   WeiSynVol = WeiSynVol / ((max(as.numeric(colnames(Matr))) - min(as.numeric(colnames(Matr)))) *
  #                            (max(as.numeric(rownames(Matr))) - min(as.numeric(rownames(Matr)))))
  #   vv = c(vv, WeiSynVol)
  # }


  # Vol = sum(sapply(2:(nrow(Matr)-1), function(i){ j = 2:(ncol(Matr)-1);(c_[i+1]/c_[i] * r_[j+1]/r_[j]) * Matr[i,j]}))
  # Vol / ((max(c_) - min(c_)) * (max(r_) - min(r_)))
  
  
  # Dr_Info <- data.frame(compounds_ = unique(c(datannot$annot$Drug1, datannot$annot$Drug2)))
  # Dr_Info$CID = mclapply(Dr_Info$compounds_, function(i){
  #   tryCatch({
  #     jsonlite::fromJSON(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",gsub(" ","",i),"/cids/JSON"))[[1]][[1]][[1]]
  #   }, error = function(e) {
  #     ""
  #   })
  # }, mc.cores = 1)
  # 
  # missingComp <- mclapply(which(Dr_Info$CID == ""), function(ind){
  #   tryCatch({
  #     jsonlite::fromJSON(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/name/",gsub(" ","",Dr_Info[ind,"compounds_"]),"/cids/JSON"))[[1]][[1]]$CID[[1]]
  #   }, error = function(e) {
  #     ""
  #   })
  # }, mc.cores = 1)
  # Dr_Info$CID[Dr_Info$CID == ""] = missingComp
  # 
  # Dr_Info$SMILES = mclapply(Dr_Info$CID, function(i){
  #   tryCatch({
  #     tmp_ = jsonlite::fromJSON(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",gsub(" ","",i),"/JSON"))$PC_Compounds$props[[1]]
  #     tmp_$value[tmp_$urn$label == "SMILES" & tmp_$urn$name == "Canonical","sval"]
  #   }, error = function(e) {
  #     ""
  #   })
  # }, mc.cores = 1)
  # 
  # library(ChemmineOB)
  # FP4fing <- do.call("rbind", mclapply(Dr_Info$SMILES, function(sm) fingerprint_OB(forEachMol("SMILES",sm[[1]],identity), "FP4")))
  # FP4fing = as.data.frame(FP4fing); rownames(FP4fing) = Dr_Info$compounds_
  # saveRDS(FP4fing,"FP4fing.RDS")
  # 
  # toxpredStruct = readRDS("tox200OurStructure.RDS")
  # toxpredStruct = data.frame(toxpredStruct= toxpredStruct, name= names(toxpredStruct))
  # toxpredStruct = toxpredStruct
  # 
  # dt_ <- openxlsx::read.xlsx("192_Combinations (10).xlsx")
  # View(dt_)
  # dt_$ToxStr = sapply(1:nrow(dt_), function(i) mean(toxpredStruct[toxpredStruct$name == dt_$Drug1[[i]],"toxpredStruct"], 
  #                                                   toxpredStruct[toxpredStruct$name == dt_$Drug2[[i]],"toxpredStruct"]))
  # 
  # 
  # cor(dt_$ToxStr, dt_$toxScores)
  # dtSyn = read.csv2("C:\\Users\\Localadmin_aianevsk\\Documents\\projects\\Anil\\192reall\\SUM_SYN_ANT.csv",sep=",")
  # dtSyn$Drug1 = dtSyn$COMPOUND_A
  # dtSyn$Drug2 = dtSyn$COMPOUND_B
  # dtSyn$ToxStr = sapply(1:nrow(dtSyn), function(i) min(toxpredStruct[toxpredStruct$name == dtSyn$Drug1[[i]],"toxpredStruct"], 
  #                                                      toxpredStruct[toxpredStruct$name == dtSyn$Drug2[[i]],"toxpredStruct"]))
  # dtSyn$PairIndex = gsub("^.*_","",dtSyn$Folder)
  # dttss = dt_[,c("PairIndex","Cell")];  dttss = dttss[!duplicated(dttss), ]
  # dtSyn = merge(dtSyn2, dttss, by = "PairIndex", all.x = T, all.y = F)
  # dtSyn$Comb = paste0(dtSyn$Drug1,"_",dtSyn$Drug2)
  # 
  # for(i in 1:nrow(dt_)){
  #   matr[rownames(matr) == dtSyn[i,"Comb"],colnames(matr) == dtSyn[i,"Cell"]] =  dtSyn[i,"SUM_SYN_ANT_LOEWE"]
  # }
  # dttss2 = dt_[,c("PairIndex","ToxStr")]; dttss2 = dttss2[!duplicated(dttss2), ]
  # dtSyn = merge(dtSyn, dttss2, by = "PairIndex", all.x = T, all.y = F)
  # 
  # matrcp = matr
  # for(i in 1:nrow(dt_)){
  #   matrcp[rownames(matrcp) == dtSyn[i,"Comb"],colnames(matrcp) == dtSyn[i,"Cell"]] =  1/dtSyn[i,"ToxStr.x"]
  # }
  # 
  # corrL = c()
  # for(i in c(3,4,5,9)){
  #   corrL = c(corrL, cor(as.numeric(as.character(matr[1:30,i])),as.numeric(as.character(matrcp[1:30,i]))))
  # }
  # corrB = c()
  # for(i in c(3,4,5,9)){
  #   corrB = c(corrB, cor(as.numeric(as.character(matr[1:30,i])),as.numeric(as.character(matrcp[1:30,i]))))
  # }
  # corrH = c()
  # for(i in c(3,4,5,9)){
  #   corrH = c(corrH, cor(as.numeric(as.character(matr[1:30,i])),as.numeric(as.character(matrcp[1:30,i]))))
  # }
  
  
  
  # matr = matrix(ncol = length(unique(dt_$Cell)), nrow = length(unique(dt_$Comb)))
  # colnames(matr) = unique(dt_$Cell); rownames(matr) = unique(dt_$Comb)
  # for(i in 1:nrow(dt_)){
  #   matr[rownames(matr) == dt_[i,"Comb"],colnames(matr) == dt_[i,"Cell"]] = 1
  # }
  # 
  # 
  # dtTr= read.csv2("C:\\Users\\Localadmin_aianevsk\\Downloads\\trainingset_171130.txt", sep='\t')
  # dtTr = na.omit(dtTr[,c(5,11)])
  # dtVal= read.csv2("C:\\Users\\Localadmin_aianevsk\\Downloads\\validationset.txt", sep='\t')
  # dtVal = na.omit(dtVal[,c(7,11)])
  # dt_ = rbind(dtTr, dtVal)
  # dt_ = dt_[dt_$LD50_mgkg!="",]
  # 
  # 
  # library(ChemmineOB)
  # FP4fing <- c(); sm2 = c();
  # for(sm in 1:length(Dr_Info$SMILES)) {
  #   sm2 = c(sm, sm2)
  #   tmp_ = fingerprint_OB(forEachMol("SMILES",Dr_Info$SMILES[[sm2]][[1]],identity), "FP4")
  #   FP4fing = rbind(FP4fing, tmp_)
  # }
  # 
  # FP4fing = as.data.frame(FP4fing);
  # FP4fing$Tox = dt_$LD50_mgkg
  # 
  # FP4fing = apply(FP4fing,1,paste,collapse="")  
  # FP4fingB = apply(FP4fingB,1,paste,collapse="") 
  # 
  # sapply(1:length(FP4fing), function(i){
  #   which(FP4fing[i] == FP4fingB)
  # }) 
  # 
  
  plot2Dbar <- function(){
    #browser();  browser();

    results_Bar = isolate(result_$result_); results_Bar$combi = paste0(results_Bar$Drug1, " & ", results_Bar$Drug2)

    if(input$icons2 == "synergy, toxicity, potency"){
      
      results_Bar$barscore = result_$result_$score; result_$result_ <<- results_Bar;
      g_ <- labs(title="Top 20 combinations", subtitle="sorted by summary score", caption="Ordered by summary score of potency, toxicity and synergy")
      
    } else if(input$icons2 == "synergy, toxicity") {
      
      results_Bar$barscore = result_$result_$scoreST;  result_$result_ <<- results_Bar;
      g_ <- labs(title="Top 20 combinations", subtitle="sorted by summary score", caption="Ordered by summary score of synergy and toxicity")
      
    } else if(input$icons2 == "synergy, potency") {
      
      results_Bar$barscore = result_$result_$scoreIS;  result_$result_ <<- results_Bar;
      g_ <- labs(title="Top 20 combinations", subtitle="sorted by summary score", caption="Ordered by summary score of synergy and potency")
      
    } else if(input$icons2 == "toxicity, potency") {
      
      results_Bar$barscore = result_$result_$scoreIT;  result_$result_ <<- results_Bar;
      g_ <- labs(title="Top 20 combinations", subtitle="sorted by summary score", caption="Ordered by summary score of toxicity and potency")
      
    } else if(input$icons2 == "synergy") {

      results_Bar$barscore = result_$result_$synScores;  result_$result_ <<- results_Bar;
      g_ <- labs(title="Top 20 combinations", subtitle="sorted by summary score", caption="Ordered by synergy")
      
    } else if(input$icons2 == "potency") {
      
      results_Bar$barscore = result_$result_$inhScores;  result_$result_ <<- results_Bar;
      g_ <- labs(title="Top 20 combinations", subtitle="sorted by summary score", caption="Ordered by potency")
 
    } else if(input$icons2 == "toxicity") {
      
      results_Bar$barscore = result_$result_$toxScores;  result_$result_ <<- results_Bar;
      g_ <- labs(title="Top 20 combinations", subtitle="sorted by summary score", caption="Ordered by toxicity")

    }
    
    results_Bar = results_Bar[order(results_Bar$barscore, decreasing = T),]; results_Bar$Position = 1:nrow(results_Bar)
    results_Bar = head(results_Bar, 20)
    ggplot(data=results_Bar, aes(x=Position, y=barscore)) + theme_classic() +
      geom_bar(stat="identity", width = 0.5, fill="#990746") + scale_x_reverse(breaks=1:20, labels = results_Bar$combi) + 
      ylim(0,1.05) + theme(
        panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), 
        legend.box.background = element_rect(fill = "transparent")
      ) + coord_flip() +
      theme(axis.text.x = element_text(margin=ggplot2::margin(0,0,0,0,"pt"),hjust=1)) + xlab("Combination") + ylab("summary score") + g_
    
  }
  
  plot2Dgr <- function(){
    result_$result_$STP_score = result_$result_$score
    if(input$icons == "synergy vs toxicity"){
      ggplot(result_$result_, aes(x=synScores, y=toxScores)) +
        geom_point(alpha = 0.95, aes(size=STP_score, color=STP_score), position=position_jitter(0.01)) +
        theme_classic() + scale_colour_gradientn(colours = colorRampPalette(c("#F7FBFF", "#990746"))(25)[4:22]) +
        geom_vline(xintercept = 0, linetype="dotted") +  guides(color= guide_legend(reverse = !0), size=guide_legend(reverse = !0)) +
        theme(
          panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), 
          legend.box.background = element_rect(fill = "transparent")
        ) + labs(title="Top 20 combinations", subtitle="synergy and toxicity", 
                 caption="Ordered by summary score of potency, toxicity and synergy") 
    } else if(input$icons == "synergy vs inhibition") {
      ggplot(result_$result_, aes(x=synScores, y=inhScores)) +
        geom_point(alpha = 0.95, aes(size=STP_score, color=STP_score), position=position_jitter(0.01)) +
        theme_classic() + scale_colour_gradientn(colours = colorRampPalette(c("#F7FBFF", "#990746"))(25)[4:22]) +
        geom_vline(xintercept = 0, linetype="dotted") +  guides(color= guide_legend(reverse = !0), size=guide_legend(reverse = !0)) +
        theme(
          panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), 
          legend.box.background = element_rect(fill = "transparent")
        ) + labs(title="Top 20 combinations", subtitle="synergy and inhibition", 
                 caption="Ordered by summary score of potency, toxicity and synergy") 
      
    } else {
      ggplot(result_$result_, aes(x=toxScores, y=inhScores)) +
        geom_point(alpha = 0.95, aes(size=STP_score, color=STP_score), position=position_jitter(0.01)) +
        theme_classic() + scale_colour_gradientn(colours = colorRampPalette(c("#F7FBFF", "#990746"))(25)[4:22]) +
        geom_vline(xintercept = 0, linetype="dotted") +  guides(color= guide_legend(reverse = !0), size=guide_legend(reverse = !0)) +
        theme(
          panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), 
          legend.box.background = element_rect(fill = "transparent")
        ) + labs(title="Top 20 combinations", subtitle="toxicity and inhibition", 
                 caption="Ordered by summary score of potency, toxicity and synergy") 
      
    }
  }
  
  ########
  # first scatter plot  
  output$scatterplot <- renderPlot({
    input$icons
    if(!is.null(isolate(result_$draw_))){
     # browser();
      plot2Dgr()
 
    }
  }, bg="transparent", height = 500)
  
  output$hover_info <- renderUI({
    
   #brower();
    
    if(!is.null(result_$draw_)){
      
      #browser();browser();

      point <- nearPoints(result_$result_, input$plot_hover, threshold = 10, maxpoints = 1, addDist = !1)
      if (nrow(point) == 0) { return(NULL)} 
      
      ind_ = which(result_$result_$Drug1==point$Drug1 & result_$result_$Drug2==point$Drug2)[[1]]

      # actual tooltip created as wellPanel
      wellPanel(
        style = paste0("position:fixed; z-index:100000; background-color: rgba(245, 245, 245, 1); ",
                       "left:2%; top:2%"),
        p(HTML(paste0("Drug1:<b>", paste0(" ",point$Drug1), "</b> ,",
                      "Drug2: <b>", paste0(" ",point$Drug2), "</b><br/>",
                      "<div class='images'>", Inh_base64[ind_], Tox_base64[ind_], "</div><br>", Syn_base64[ind_])))
      )
    }
  })
  
  # ########
  # # bar plot  
  
  output$barplot <- renderPlot({
    if(!is.null(result_$draw_)){ 
      
      plot2Dbar();
      
    }
  }, bg="transparent", height = 500)

  
  output$hover_info2 <- renderUI({

    #browser();
    if(!is.null(result_$draw_)){
      #browser();
      hover <- input$plot_hover2
      
      #browser();browser();

      #point <- nearPoints(result_$result_, hover, threshold = 10, maxpoints = 1, addDist = !1)
      if(is.null(hover$y) || abs(hover$y) > 20 || abs(hover$y) < 0) return(NULL)
      #browser();
      ind_ = -1*round(hover$y); order_ = order(result_$result_$barscore, decreasing = T); point = result_$result_[order_,][ind_,]
      
      #which(result_$result_$Drug1==point$Drug1 & result_$result_$Drug2==point$Drug2)[[1]]

      # actual tooltip created as wellPanel
      wellPanel(
        style = paste0("position:fixed; z-index:100000; background-color: rgba(245, 245, 245, 1); ",
                       "left:2%; top:2%"),
        p(HTML(paste0("Drug1:<b>", paste0(" ",point$Drug1), "</b> ,",
                      "Drug2: <b>", paste0(" ",point$Drug2), "</b><br/>",
                      "<div class='images'>", Inh_base64[order_][ind_], Tox_base64[order_][ind_], "</div><br>", Syn_base64[order_][ind_])))
      )
    }
  })
  
  ########
  # 3d plot
  
  output$plotlcont <- renderUI({
    if(!is.null(result_$draw_)){
      plotlyOutput("plotl", height = "700px")
    }
  })
  
  output$plotl <- renderPlotly({
      plot_ly(x = result_$result_$synScores, y = result_$result_$inhScores, z = result_$result_$toxScores, type = "scatter3d",  color = result_$result_$score,
              text = paste0("Synergy: <b>", result_$result_$synScores, "</b><br>Toxicity: <b>", result_$result_$toxScores, 
                            "</b><br>Potency: <b>",result_$result_$inhScores,"</b><br>Combined score: <b>", round(result_$result_$score,2),"</b>"), hoverinfo = 'text') %>%  
        layout(plot_bgcolor='rgba(254, 247, 234,0)', paper_bgcolor='rgba(254, 247, 234,0)', title = "Combined Synergy, Toxicity, and Potency landscape",
               scene = list(xaxis = list(title = "Synergy", gridcolor = "#696969"), 
                            yaxis = list(title = "Potency", gridcolor = "#696969"), 
                            zaxis = list(title = "Toxicity", gridcolor = "#696969"))) 
  })
  
  output$hover_info3 <- renderUI({
    d <- event_data("plotly_hover");
    if (is.null(d)) return(NULL)

    ind_ = d$pointNumber+1; point = result_$result_[ind_,]
      
      # actual tooltip created as wellPanel
      wellPanel(
        style = paste0("position:fixed; z-index:100000; background-color: rgba(245, 245, 245, 1); ",
                       "left:2%; top:2%"),
        p(HTML(paste0("Drug1:<b>", paste0(" ",point$Drug1), "</b> ,",
                      "Drug2: <b>", paste0(" ",point$Drug2), "</b><br/>",
                      "<div class='images'>", Inh_base64[ind_], Tox_base64[ind_], "</div><br>", Syn_base64[ind_])))
      )
  })
  
  ########
  # Additional funct.
  library(reshape2); library(ggplot2)
  
  # plot dose-response matrix (args: matrix, title, agent names)
  PlotRespMatr <- function(matr, name = "none", d1N = "", d2N = "", col_ = "#b44545"){
    
    data.plot <- melt(matr)
    colnames(data.plot) <- c("y","x","Inhibition")
    data.plot$Inhibition <- round(c(matr))
    data.plot$x <- as.factor(data.plot$x)
    data.plot$y <- as.factor(data.plot$y)
    
    axis.x.text <- round(as.numeric(colnames(matr)), 1)
    axis.y.text <- round(as.numeric(rownames(matr)), 1)
    dose.response.p <- ggplot(data.plot, aes_string(x = "x", y = "y")) + geom_tile(aes_string(fill = "Inhibition"), size=1) + 
      theme(title = element_text(face = "bold", size = 10)) + 
      geom_text(aes_string(fill = "Inhibition", label = "Inhibition"), size = 3.5) + 
      scale_fill_gradient2(low = "grey", high = col_, midpoint = 0, name = paste0("% cell \ninhibition"),na.value="white", limits=c(0, 100)) + 
      scale_x_discrete(labels = axis.x.text) + scale_y_discrete(labels = axis.y.text) + 
      labs(x = d2N, y = d1N) + theme(plot.title = element_text(hjust = 0.5, size = 16)) +
      theme(
        panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent")
      )
    
    dose.response.p <- dose.response.p + theme(axis.text.x = element_text(color = "black", face = "bold", size = 11))
    dose.response.p <- dose.response.p + theme(axis.text.y = element_text(color = "black", face = "bold", size = 11))
    dose.response.p <- dose.response.p + theme(axis.title = element_text(size = 14))
    dose.response.p <- dose.response.p + ggtitle(paste0("\nDose-response matrix (",name,")\n"));
    list(pl = dose.response.p)
  }
  
  # plot synergy matrix (args: matrix, title, agent names)
  PlotSynMatr <- function(matr, name = "none", d1N = "", d2N = "", col_ = "#b44545"){
    
    data.plot <- melt(matr)
    colnames(data.plot) <- c("y","x","z")
    data.plot$z <- round(c(matr),2)
    myPalette <- colorRampPalette(c("green2", "white", "red1"))(100)
    
    dt_ = calcsyn(matr, data.frame(concUnit="nM", drug.row=d1N, drug.col=d2N));
    plot2d = melt(dt_$c)
    names(plot2d) <- c("x","y","z")
    byx = (max(plot2d$x) - min(plot2d$x))/(length(dt_$x.conc) - 1);
    byy = (max(plot2d$y) - min(plot2d$y))/(length(dt_$y.conc) - 1);
    
    syn.plot <- ggplot(plot2d) + aes(x, y, z = z, fill = z)  + geom_raster(interpolate = !0) + 
      geom_contour(color = "white", alpha = 0.5) + 
      scale_fill_gradientn(expression(delta ~ -score), colours = myPalette, limits = c(dt_$start.point, dt_$end.point)) +
    scale_x_continuous(dt_$drug.col, expand = c(0, 0), 
                       breaks = seq(min(plot2d$x), max(plot2d$x), by = (max(plot2d$x) - min(plot2d$x))/(length(dt_$x.conc) - 1)), 
                       labels = round(dt_$x.conc, 2)) + 
    scale_y_continuous(dt_$drug.row, expand = c(0, 0), 
                       breaks = seq(min(plot2d$y), max(plot2d$y), by = (max(plot2d$y) - min(plot2d$y))/(length(dt_$y.conc) - 1)), 
                       labels = round(dt_$y.conc, 2)) +
      theme_bw() + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      theme(axis.text = element_text(size = 10)) + 
      theme(title = element_text(vjust = 12)) + 
      theme(plot.title = element_text(size = 18, margin = ggplot2::margin(b = 25, unit = "pt"))) + 
      ggtitle(paste0("Synergy distribution (",input$state,")"))  + theme(plot.title = element_text(hjust = 0.5)) +
      geom_vline(xintercept = seq(min(plot2d$x), max(plot2d$x), by = byx), linetype = "dotted") + 
      geom_hline(yintercept = seq(min(plot2d$y), max(plot2d$y), by = byy), linetype = "dotted")  +
      theme(
        panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent")
      )

    list(pl = syn.plot)
    }
  
  # Transform input data to appropriate for synergyfinder
  transformInputData <- compiler::cmpfun(function(dataOutput, data_typ = NULL){
    
    if (data_typ == "viability") dataOutput$Response = 100 - dataOutput$Response
    
    dose.response.mats <- list()
    num.pairs <- length(unique(dataOutput$PairIndex)); pairs_ = unique(dataOutput$PairIndex);
    drug.pairs = data.frame(drug.row = 1:num.pairs, drug.col = 1:num.pairs, concUnit = 1:num.pairs, PairIndex = 1:num.pairs)
    error_ <- rep(1, num.pairs); # errors initially 1
    warning_ <- rep("", num.pairs); # warnings initially ""
    
    for(i in 1:num.pairs){
      df2 <- dplyr::arrange(dataOutput[dataOutput$PairIndex==pairs_[i],],Conc1,Conc2)
      Conc1 <- sort(as.numeric(unique(df2$Conc1))); Conc2 <- sort(as.numeric(unique(df2$Conc2)));
      
      tryCatch({
        aaaa = matrix(nrow = length(Conc1), ncol = length(Conc2), byrow = T)
        dimnames(aaaa) = list(Conc1, Conc2) # rowname, colname
        
        for(j in 1:length(Conc1)){
          for(k in 1:length(Conc2)){ 
            resp = as.numeric(df2[((df2$Conc1 == Conc1[j]) & (df2$Conc2 == Conc2[k])),"Response"]);
            if(!identical(resp, numeric(0))) aaaa[j,k] = resp else aaaa[j,k] = NA
          }
        }
        ###dose.response.mats[[i]] = matrix(as.numeric(df2$Response), nrow = length(Conc1), ncol = length(Conc2), byrow = T)
        dose.response.mats[[i]] = aaaa;
        dimnames(dose.response.mats[[i]]) = list(Conc1, Conc2) # rowname, colname
        drug.pairs[i,]$drug.row <- df2$Drug1[1]; drug.pairs[i,]$drug.col<- df2$Drug2[1];
        drug.pairs[i,]$concUnit <- df2$ConcUnit[1]; drug.pairs[i,]$PairIndex <- df2$PairIndex[1];
        error_[i] <- 0;
        if(length(Conc1) * length(Conc2) != length(df2$Response))
          warning_[i] <- paste0("The provided number of responses for ",i," PairIndex is: <b>",length(df2$Response), "</b>. However <b>",length(Conc1) * length(Conc2), "</b> responses were expected. Please reconsider <b>",df2$Drug1[1]," & ",df2$Drug2[1], "</b> drug pair, or continue at your own risk " )
      })
    }
    return(list(dose.response.mats = dose.response.mats, drug.pairs = drug.pairs, error = error_, warning = warning_))
  })
  
  # Dose-response curve-fitting
  CALC_IC50_EC50_DSS <- function(xpr_tbl, DSS_typ, readoutCTX = F, drug_name="drug_name")
  {
    library(drc)
    tryCatch({
      
      mat_tbl <- data.frame(inhibition=as.numeric(xpr_tbl), dose = as.numeric(names(xpr_tbl))); #dose = 10**(1:length(xpr_tbl)));
      mat_tbl$logconc = log10(mat_tbl$dose); mat_tbl$viability = 100 - mat_tbl$inhibition;
      mat_tbl$inhibition2 = mat_tbl$inhibition; mat_tbl$viability2 = mat_tbl$viability;
      mat_tbl <- mat_tbl[order(mat_tbl[,"dose"]),] 
      
      if(any(duplicated(mat_tbl$inhibition))) mat_tbl$inhibition <- seq(from = 0, length.out = length(mat_tbl$inhibition), by = 0.01) + mat_tbl$inhibition; 
      
      
      estimate_param <- tryCatch({drm(inhibition ~ logconc, data = mat_tbl, fct = LL.4(fixed = c(NA, NA, NA,NA),names = c("SLOPE","MIN","MAX","IC50")),logDose=10,control = drmc(errorm = F))}, 
                                 warning=function(w){drm(inhibition ~ logconc, data = mat_tbl, fct = L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),logDose=10)},
                                 error=function(e){drm(inhibition ~ logconc, data = mat_tbl, fct = L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),logDose=10)})
      coef_estim <- coef(estimate_param); names(coef_estim) <- c("SLOPE","MIN","MAX","IC50"); coef_estim["SLOPE"] <- coef_estim["SLOPE"]*-1 
      coef_estim["IC50"] <- ifelse(coef_estim["MAX"]<=coef_estim["MIN"] | coef_estim["IC50"]>max(mat_tbl$dose,na.rm=T), max(mat_tbl$dose,na.rm=T),coef_estim["IC50"])
      coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<0,min(mat_tbl$dose,na.rm=T),coef_estim["IC50"]); coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<0,mean(mat_tbl$dose,na.rm=T),coef_estim["IC50"])
      coef_estim["IC50"] <- log10(coef_estim["IC50"]); coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<min(mat_tbl$logconc),max(mat_tbl$logconc),coef_estim["IC50"])
      coef_estim["IC50"] <- ifelse(all(mat_tbl$inhibition<0),max(mat_tbl$logconc,na.rm=T),coef_estim["IC50"]); coef_estim["MIN"] <- 0; coef_estim["MAX"] <- max(mat_tbl$inhibition,na.rm=T)
      min_lower <- ifelse(min(mat_tbl$inhibition,na.rm=T) > 0,min(mat_tbl$inhibition,na.rm=T),0); min_lower <- ifelse(min_lower >= 100,99,min_lower)
      coef_estim["MAX"] <- ifelse(coef_estim["MAX"]>100,100,coef_estim["MAX"]); coef_estim["MAX"] <- ifelse(coef_estim["MAX"]<0,100,coef_estim["MAX"])
      max_lower <- ifelse(max(mat_tbl$inhibition,na.rm=T)>100,coef_estim["MAX"],max(mat_tbl$inhibition,na.rm=T)); max_lower <- ifelse(max_lower < 0,coef_estim["MAX"],max(mat_tbl$inhibition,na.rm=T)); 
      max_lower <- ifelse(max_lower < 0,0,max_lower); max_lower <- ifelse(max_lower > 100,100,max_lower); run_avg <- caTools::runmean(mat_tbl$inhibition, 10); max_upper <- ifelse(any(run_avg[-nrow(mat_tbl)]>run_avg[nrow(mat_tbl)]),max(mat_tbl$inhibition[run_avg>run_avg[nrow(mat_tbl)]]),coef_estim["MAX"])
      max_upper <- ifelse(any(mat_tbl$inhibition > max_upper),mean(mat_tbl$inhibition[mat_tbl$inhibition > max_upper])+5,max_upper)
      max_upper <- ifelse(max_upper < 0,coef_estim["MAX"],max_upper); max_upper <- ifelse(max_upper > 100,100,max_upper) #coef_estim["MAX"]
      max_upper <- ifelse(max_lower > max_upper,coef_estim["MAX"],max_upper);mean_inh_last = mean(tail(mat_tbl$inhibition,2),na.rm=T)
      if(mean_inh_last < 60) {
        if(mean_inh_last > 25) coef_estim["IC50"] <- mean(mat_tbl$logconc,na.rm=T) else if(mean_inh_last < 25) coef_estim["IC50"] <- max(mat_tbl$logconc,na.rm=T)}
      if(mean(mat_tbl$inhibition[1:3],na.rm=T)<5) coef_estim["IC50"] <- max(mat_tbl$logconc,na.rm=T)
      if(unname(coef_estim["MIN"]) == unname(coef_estim["MAX"])) coef_estim["MAX"] <- coef_estim["MAX"] + 0.001
      
      #adaptive nonlinear Least-Squares algorithm NL2SOL to estimate parameters.
      nls_result_ic50_old <- function(){
        tryCatch({
          nls(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl, algorithm="port", start=list(SLOPE=1,MIN=coef_estim["MIN"][[1]],MAX=coef_estim["MAX"][[1]],IC50=coef_estim["IC50"][[1]]),lower=list(SLOPE=0,MIN=0,MAX=max_lower, IC50=min(mat_tbl$logconc)),upper=list(SLOPE=4,MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)),control=list(warnOnly=T,minFactor = 1/2048))
        }, error = function(e) {minpack.lm::nlsLM(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl,start=list(SLOPE=1, MIN=coef_estim["MIN"][[1]],MAX=coef_estim["MAX"][[1]],IC50=coef_estim["IC50"][[1]]),lower=c(SLOPE=0, MIN=0,MAX=max_lower, IC50=min(mat_tbl$logconc)),upper=c(SLOPE=4, MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)))
        })} 
      # IC50 first
      nls_result_ic50 <- nls_result_ic50_old();
      nls_result_ic50_2 <- tryCatch({
        nls(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl, algorithm="port", start=list(SLOPE=1,MIN=coef_estim["MIN"][[1]],MAX=coef_estim["MAX"][[1]], IC50=median(mat_tbl$logconc)),lower=list(SLOPE=0,MIN=0,MAX=max_lower, IC50=min(mat_tbl$logconc)),upper=list(SLOPE=4,MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)),control=list(warnOnly=T,minFactor = 1/2048))
      },warning = function(w) {nls_result_ic50_old()},error = function(e) {nls_result_ic50_old()})
      
      aaa=tryCatch({summary(nls_result_ic50)},error=function(e){summary(nls_result_ic50_2)})
      bbb=tryCatch({summary(nls_result_ic50_2)},error=function(e){summary(nls_result_ic50)})
      
      sumIC50 = list(aaa,bbb)
      ic50std_resid <- round(sqrt(sum((sumIC50[[1]]$residuals)^2)/(length(sumIC50[[1]]$residuals)-1)),1);
      ic50std_resid2 <- round(sqrt(sum((sumIC50[[2]]$residuals)^2)/(length(sumIC50[[2]]$residuals)-1)),1);
      # continue with the best
      switch_ = which.min(c(ic50std_resid, ic50std_resid2)); nls_result_ic50 = list(nls_result_ic50, nls_result_ic50_2)[[switch_]]
      
      if(coef(nls_result_ic50)["SLOPE"] <= 0.2){if(mean_inh_last > 60) coef_estim["IC50"] <- min(mat_tbl$logconc,na.rm=T);
      nls_result_ic50 <- nls(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl, algorithm="port",start=list(SLOPE=1, MIN=unname(coef_estim["MIN"]),MAX=unname(coef_estim["MAX"]),IC50=unname(coef_estim["IC50"])),lower=list(SLOPE=0.1,MIN=min_lower,MAX=max_lower,IC50=min(mat_tbl$logconc)),upper=list(SLOPE=2.5, MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)),control=list(warnOnly=T,minFactor = 1/2048))
      }
      #Calculate the standard error scores
      sumIC50 = summary(nls_result_ic50); ic50std_Error <- sumIC50$coefficients["IC50","Std. Error"];
      ic50std_resid <- round(sqrt(sum((sumIC50$residuals)^2)/(length(sumIC50$residuals)-1)),1);
      max_signal <- max(mat_tbl$dose,na.rm=T); min_signal <- min(mat_tbl$dose,na.rm=T)
      
      #prepare final data and convert IC50 back from log scale (inverse)
      coef_ic50 <- coef(nls_result_ic50)[c("IC50", "SLOPE","MAX","MIN")]; coef_ic50["IC50"] <- 10^coef_ic50["IC50"]; coef_ic50["IC50"] <- ifelse(coef_ic50["SLOPE"]<0,max_signal,coef_ic50["IC50"])
      coef_ic50["IC50"] <- ifelse(coef_ic50["MAX"]<0,max_signal,coef_ic50["IC50"]);coef_ic50["IC50"] <- ifelse(coef_ic50["MAX"]<10,max_signal,coef_ic50["IC50"])
      coef_ic50["MAX"] <- ifelse(coef_ic50["MAX"]<0,0,coef_ic50["MAX"]);coef_ic50["IC50"] <- ifelse(all(c(max(mat_tbl$inhibition,na.rm=T),min(mat_tbl$inhibition,na.rm=T))>50),min_signal,coef_ic50["IC50"])
      x <- seq(min(mat_tbl$logconc),max(mat_tbl$logconc), length=100); yic <- predict(nls_result_ic50, data.frame(logconc=x))
      perInh <- t(matrix(mat_tbl[,"inhibition"],dimnames=list(paste0(rep("D", length(mat_tbl[,"inhibition"])), 1:length(mat_tbl[,"inhibition"])))))
      coef_tec50 = coef_ic50; 
      coef_tec50["IC50"] <- ifelse(coef_tec50["MAX"] > 25, coef_tec50["IC50"], max(mat_tbl$dose,na.rm=T))
      if(readoutCTX){names(coef_tec50) <- c("TC50","SLOPE","MAX","MIN"); ytec <- yic; perViaTox <- perInh;} else{
        names(coef_tec50) <- c("EC50","SLOPE","MAX","MIN"); coef_tec50["SLOPE"] = -1 * coef_tec50["SLOPE"]; # min - 0, max - 77 in ec50 it is max - 100, min - 23
        tmp = coef_tec50["MAX"]; coef_tec50["MAX"] = 100 - coef_tec50["MIN"]; coef_tec50["MIN"] = 100 - tmp; ytec <- 100 - yic;
        perViaTox <- 100 - perInh;
      }
      ############################# 
      #############    DSS
      coef_ic50 <- c(coef_ic50,Min.Conc.tested=min_signal,Max.Conc.tested=max_signal,IC50_std_error=ic50std_Error,DSS=100)
      coef_tec50 <- c(coef_tec50,Min.Conc.tested=min_signal,Max.Conc.tested=max_signal,TEC50_std_error=ic50std_Error)
      IC50_df <- data.frame(DRUG_NAME="drug_name",ANALYSIS_NAME="IC50", t(as.matrix(coef_ic50)), perInh,GRAPH=NA, DSS = 100, sDSS = "", SE_of_estimate = as.numeric(ic50std_resid))
      
      return (list(coef_ic50=coef_ic50,nls_result_ic50=nls_result_ic50));
    })
  }
  
  
  
  calcsyn <- function(scores, drug.pairs)
  {
    scores.dose <- t(scores)  
    combscores = scores.dose[-1, -1]; combscores[nrow(combscores),ncol(combscores)] <- 'NA'
    summary.score <- round(mean(as.numeric(combscores), na.rm = !0), 3)
    x.conc <- as.numeric(rownames(scores.dose))
    y.conc <- as.numeric(colnames(scores.dose))
    conc.unit <- drug.pairs$concUnit
    unit.text <- paste0("(", conc.unit, ")")
    drug.row <- paste0(drug.pairs$drug.row, " ", unit.text)
    drug.col <- paste0(drug.pairs$drug.col, " ", unit.text)
    color.range <- round(max(abs(max(as.numeric(combscores), na.rm = !0)), abs(min(as.numeric(combscores), na.rm = !0))) + 5, -1)
    start.point <- -color.range; end.point <- color.range  
    pixels.num = 5 * (length(x.conc) - 1) + 2;
    
    # only for visualization  (max BzCl)
    scores.dose[nrow(scores.dose),ncol(scores.dose)] <- max(scores.dose[nrow(scores.dose)-1,ncol(scores.dose)],
                                                            scores.dose[nrow(scores.dose),ncol(scores.dose)-1],
                                                            scores.dose[nrow(scores.dose)-1,ncol(scores.dose)-1])    
    
    kriged =  tryCatch({
      tmp <- cbind(expand.grid(c(0:(length(x.conc) - 1)), c(0:(length(y.conc) - 1))), c(as.matrix(scores.dose)))        
      kriging(tmp[, 1],tmp[, 2], tmp[, 3], lags = ifelse(dim(scores.dose)[1] < 8, 2 ,3), 
              pixels = pixels.num, model = "spherical")
    },error = function(e){
      appro <- function(x, n) approx(x, n=n)$y
      tryCatch({
        m = apply(t(apply(scores.dose, 1, function(x) appro(x, ncol(scores.dose)*2))), 2, function(x) appro(x, nrow(scores.dose)*2))
        tmp2<- cbind(expand.grid(c(0:(nrow(m)-1)),c(0:(ncol(m)-1))), c(as.matrix(m)))
        kriging(tmp2[, 1], tmp2[, 2], tmp2[, 3], lags = ifelse(dim(m)[1] < 8, 2 ,3), pixels = pixels.num, model = "spherical")
      },error = function(e){ 
        m = apply(t(apply(scores.dose, 1, function(x) appro(x, ncol(scores.dose)*3))), 2, function(x) appro(x, nrow(scores.dose)*3))
        tmp2<- cbind(expand.grid(c(0:(nrow(m)-1)),c(0:(ncol(m)-1))), c(as.matrix(m)))
        kriging(tmp2[, 1], tmp2[, 2], tmp2[, 3], lags = ifelse(dim(m)[1] < 8, 2 ,3), pixels = pixels.num, model = "spherical")
      })
    })
    
    xseq <- round(kriged[["map"]]$x/kriged$pixel)
    yseq <- round(kriged[["map"]]$y/kriged$pixel)
    a <- min(xseq):max(xseq); b <- min(yseq):max(yseq)
    na <- length(a); nb <- length(b)
    res1 <- as.double(rep(0, na * nb))
    res2 <- as.integer(rep(0, na * nb))

    z.len = length(kriged[["map"]]$pred); #
    for(idx1 in 1:na) { #
      for(idx2 in 1:nb) { #
        for(idx3 in 1:z.len) { #
          if(xseq[idx3] == a[idx1] && yseq[idx3] == b[idx2]) { #
            indx_ = idx2+(idx1-1)*nb; #
            res1[indx_] <- kriged[["map"]]$pred[idx3] #
            res2[indx_] <- 1 #
            break #
          } #
        } #
      } #
    } #
    
    res1[res2 == 0] <- NA #
    cMat <- matrix(res1, na, nb, byrow = !0)
    
    # most synergystic region
    max_ = r_ = c_ = -999;
    for(i in 1:(ncol(scores.dose)-2)){
      for(j in 1:(nrow(scores.dose)-2)){
        mean_ = mean(scores.dose[j:(j+2),i:(i+2)], na.rm = !0)
        if(mean_ > max_) {
          max_ = mean_; r_ = j; c_ = i;
        }
      }
    }
    
    return(list(c = cMat, conc.unit = conc.unit, drug.row = drug.row, 
                drug.col = drug.col, start.point = start.point, end.point = end.point, 
                summary.score = summary.score, x.conc = x.conc, y.conc = y.conc, pixels.num = pixels.num, r_ = r_, c_ = c_, max_ = round(max_,3)))
  }
  
}




