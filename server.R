options(java.parameters = "-Xss2560k")
library(shiny);  

server <- function(input, output, session){
  
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
            dataOutput$Response[dataOutput$Response>100] = 100;
            datannot$annot <- dataOutput;
            browser();browser()
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

    if(!is.null(datannot$annot)){
      runjs("$('#video').remove(); $('#barplot').show(); $('#col1id').show();");

      withProgress(message = 'Calculations in process', value = 0, {
        library(parallel); library(synergyfinder); library(kriging)
        
  
        # if no toxicity should be calculated
        if(input$conrol_ == "SynTox"){
          updateRadioButtons(session = session, "icons2", choices = list("efficacy & synergy"),
                             selected = "efficacy & synergy");
          updateRadioButtons(session = session, "icons", choices = list("efficacy & synergy"),
                             selected = "efficacy & synergy");
          shinyjs::runjs("$('.tabbable ul li:nth-child(3)').hide()");
        }

        #browser()#input$conrol_# load sample-control
        sample_ = datannot$annot[datannot$annot$Type=="sample", ]; 
        
        # unique combinations
        combis_ = unique(paste0(sample_$Drug1,"|",sample_$Drug2));
        sample_$PairIndex = match(paste0(sample_$Drug1,"|",sample_$Drug2), combis_);
        incProgress(1/3, detail = paste0("Quantifying synergy", paste0(" (",input$state," method)")))
        d_ = do.call("rbind", lapply(1:length(unique(sample_$PairIndex)), function(i) sample_[sample_$PairIndex == i,][1,]))
     
     
        # calculate synergy for all combinations using synergyfinder R package
        scores <<- synergyfinder::CalculateSynergy(transformInputData(sample_, "inhibition"),input$state);
        MostSyn_r <<- c(); MostSyn_c <<- c();
        synScores <- round(sapply(1:length(combis_), function(i){
          p = PlotSynMatr(scores$scores[[i]], "inhibition", d1N = scores$drug.pairs[i,1], d2N = scores$drug.pairs[i,2]);
          MostSyn_r <<- c(MostSyn_r, p$r_); MostSyn_c <<- c(MostSyn_c, p$c_); p$r_ = p$c_ = NULL;
          c_ = as.numeric(colnames(scores$scores[[i]])); r_ = as.numeric(rownames(scores$scores[[i]])); Matr = scores$scores[[i]];
          filename_ = file.path(paste0("./plots/",scores$drug.pairs[i,1],"_syn_",scores$drug.pairs[i,2],".png"))
          png(filename = filename_,width=320,height=320, bg = "transparent")
          print(p); dev.off()
          Vol = sum(unlist(sapply(3:nrow(Matr), function(k){j = 3:(ncol(Matr));v_=(log(c_[j]/c_[j-1]) * log(r_[k]/r_[k-1])) * Matr[k,j]; 
          if(k==nrow(Matr)) {v_ = v_[-length(v_)]}; v_}))); Vol / (log(max(c_) / min(c_[-1])) * log(max(r_)/min(r_[-1])))
          
        }), 3)
        # calculate normalized volume under each synergy matrix (At The Most Synergistic Area)
        synScoresMSA <- round(unlist(mclapply(1:length(combis_), function(i){
          Matr = scores$scores[[i]][c(1,MostSyn_c[i]:(MostSyn_c[i]+2)), c(1,MostSyn_r[i]:(MostSyn_r[i]+2))];
          c_ = as.numeric(colnames(Matr)); r_ = as.numeric(rownames(Matr))
          Vol = sum(unlist(sapply(3:nrow(Matr), function(k){j = 3:(ncol(Matr));v_=(log(c_[j]/c_[j-1]) * log(r_[k]/r_[k-1])) * Matr[k,j]; 
          if(k==nrow(Matr)) {v_ = v_[-length(v_)]}; v_}))); Vol / (log(max(c_) / min(c_[-1])) * log(max(r_)/min(r_[-1])))
        }, mc.cores = 1)), 3)
        
        
        incProgress(1/3, detail = "Quantifying efficacy")
  
        # calculate normalized volume under each inhibition matrix
        inhScores <- round(unlist(mclapply(1:length(combis_), function(i){
          Matr = scores$dose.response.mats[[i]]; c_ = as.numeric(colnames(Matr)); r_ = as.numeric(rownames(Matr))
          p = PlotRespMatr(Matr, "inhibition", d1N = scores$drug.pairs[i,1], d2N = scores$drug.pairs[i,2])
          filename_ = file.path(paste0("./plots/",scores$drug.pairs[i,1],"_inh_",scores$drug.pairs[i,2],".png"))
          png(filename = filename_,width=320,height=320, bg = "transparent")
          print(p); dev.off()
          Vol = sum(unlist(sapply(3:nrow(Matr), function(k){j = 3:(ncol(Matr));v_=(log(c_[j]/c_[j-1]) * log(r_[k]/r_[k-1])) * Matr[k,j]; 
          if(k==nrow(Matr)) {v_ = v_[-length(v_)]}; v_}))); Vol / (log(max(c_) / min(c_[-1])) * log(max(r_)/min(r_[-1])))
        }, mc.cores = 1)), 3)
        incProgress(1/3, detail = "Quantifying toxicity")
        # calculate normalized volume under each inhibition matrix (At The Most Synergistic Area)
        inhScoresMSA <- round(unlist(mclapply(1:length(combis_), function(i){
          Matr = scores$dose.response.mats[[i]][c(1,MostSyn_c[i]:(MostSyn_c[i]+2)), c(1,MostSyn_r[i]:(MostSyn_r[i]+2))];
          c_ = as.numeric(colnames(Matr)); r_ = as.numeric(rownames(Matr))
          Vol = sum(unlist(sapply(3:nrow(Matr), function(k){j = 3:(ncol(Matr));v_=(log(c_[j]/c_[j-1]) * log(r_[k]/r_[k-1])) * Matr[k,j]; 
          if(k==nrow(Matr)) {v_ = v_[-length(v_)]}; v_}))); Vol / (log(max(c_) / min(c_[-1])) * log(max(r_)/min(r_[-1])))
        }, mc.cores = 1)), 3)
        
  
        

        incProgress(1/3, detail = "Quantifying toxicity")
  
        # calculate normalized volume under each toxicity matrix
        if(input$conrol_ == "Control"){
          incProgress(0, detail = "Quantifying toxicity from control data")
          control_ = datannot$annot[datannot$annot$Type=="control", ]
          control_$PairIndex = match(paste0(control_$Drug1,"|",control_$Drug2), combis_);
          dt_ <- transformInputData(control_, "toxicity");
          toxScores = round(sapply(1:length(combis_), function(i){
            Matr = dt_$dose.response.mats[[i]]; c_ = as.numeric(colnames(Matr)); r_ = as.numeric(rownames(Matr))
            p = PlotRespMatr(Matr, "toxicity", d1N = scores$drug.pairs[i,1], d2N = scores$drug.pairs[i,2])
            filename_ = file.path(paste0("./plots/",scores$drug.pairs[i,1],"_tox_",scores$drug.pairs[i,2],".png"))
            png(filename = filename_,width=320,height=320, bg = "transparent")
            print(p); dev.off()
            Vol = sum(unlist(sapply(3:nrow(Matr), function(k){j = 3:(ncol(Matr));v_=(log(c_[j]/c_[j-1]) * log(r_[k]/r_[k-1])) * Matr[k,j]; 
            if(k==nrow(Matr)) {v_ = v_[-length(v_)]}; v_}))); Vol / (log(max(c_) / min(c_[-1])) * log(max(r_)/min(r_[-1])))
          }), 3)
          toxScoresMSA <- round(unlist(mclapply(1:length(combis_), function(i){
            Matr = dt_$dose.response.mats[[i]][c(1,MostSyn_c[i]:(MostSyn_c[i]+2)), c(1,MostSyn_r[i]:(MostSyn_r[i]+2))];
            c_ = as.numeric(colnames(Matr)); r_ = as.numeric(rownames(Matr))
            Vol = sum(unlist(sapply(3:nrow(Matr), function(k){j = 3:(ncol(Matr));v_=(log(c_[j]/c_[j-1]) * log(r_[k]/r_[k-1])) * Matr[k,j]; 
            if(k==nrow(Matr)) {v_ = v_[-length(v_)]}; v_}))); Vol / (log(max(c_) / min(c_[-1])) * log(max(r_)/min(r_[-1])))
          }, mc.cores = 1)), 3)
        } else if(input$conrol_ == "SynTox"){
          incProgress(0, detail = "No toxicity assessment")
          toxScores <- toxScoresMSA <-  rep(0, length(combis_));
        } 

        # results 
        results_ <- data.frame(synScores = rank(synScores) / length(combis_), inhScores = rank(inhScores) / length(combis_), 
                               s_inhScores = rank(c(inhScores - toxScores)) / length(combis_),
                               synScoresReal = synScores, inhScoresReal = inhScores, toxScoresReal = toxScores, s_inhReal = inhScores - toxScores,
                               synScoresMSA = rank(synScoresMSA) / length(combis_), inhScoresMSA = rank(inhScoresMSA) / length(combis_), 
                               s_inhScoresMSA = rank(c(inhScoresMSA - toxScoresMSA)) / length(combis_),
                               synScoresRealMSA = synScoresMSA, inhScoresRealMSA = inhScoresMSA, toxScoresRealMSA = toxScoresMSA, s_inhRealMSA = inhScoresMSA - toxScoresMSA)
        results_$score = round((results_$s_inhScores + results_$synScores) / 2, 3)
        results_$scoreIS = round((results_$inhScores + results_$synScores) / 2, 3)
        results_$scoreMSA = round((results_$s_inhScoresMSA + results_$synScoresMSA) / 2, 3)
        results_$scoreISMSA = round((results_$inhScoresMSA + results_$synScoresMSA) / 2, 3)

    
        browser();browser();browser();
        
        # add drugs
        results_$Drug1 <- d_$Drug1; results_$Drug2 <- d_$Drug2;

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
        
        save(result_$result_, result_$draw_, result_, Inh_base64, Tox_base64, Syn_base64, file = "results_.RDATA")
        
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
     browser(); browser();browser();browser();
     file_ = result_$result_; 
     file_2 = result_$result_; 
     
       if(input$conrol_ == "SynTox"){
         file_ = file_[ ,c("Drug1","Drug2","scoreIS","synScoresReal","inhScoresReal","toxScoresReal", "scoreIS")]
         colnames(file_) = c("Drug1","Drug2","STE score (selective efficacy & synergy)", paste0("Synergy volume score (",input$state," model)"),"Efficacy volume score",
                             "Toxicity volume score", "STE score (efficacy & synergy)")
         file_$`Toxicity volume score` = NA; file_$`STE score (selective efficacy & synergy)` = NA
         
       } else {
         file_ = file_[ ,c("Drug1","Drug2","score","synScoresReal","inhScoresReal","toxScoresReal","scoreIS")]
         colnames(file_) = c("Drug1","Drug2","STE score (selective efficacy & synergy)", paste0("Synergy volume score (",input$state," model)"),"Efficacy volume score",
                             "Toxicity volume score", "STE score (efficacy & synergy)")
       }

       if(input$conrol_ == "SynTox"){
         file_2 = file_2[ ,c("Drug1","Drug2","scoreISMSA","synScoresRealMSA","inhScoresRealMSA","toxScoresRealMSA","scoreISMSA")]
         colnames(file_2) = c("Drug1","Drug2","STE score (selective efficacy & synergy)", paste0("Synergy volume score (",input$state," model)"),"Efficacy volume score",
                              "Toxicity volume score", "STE score (efficacy & synergy)")
         file_2$`Toxicity volume score` = NA; file_2$`STE score (selective efficacy & synergy)` = NA
         
         
       } else {
         file_2 = file_2[ ,c("Drug1","Drug2","scoreMSA","synScoresRealMSA","inhScoresRealMSA","toxScoresRealMSA","scoreISMSA")]
         colnames(file_2) = c("Drug1","Drug2","STE score (selective efficacy & synergy)", paste0("Synergy volume score (",input$state," model)"),"Efficacy volume score",
                             "Toxicity volume score", "STE score (efficacy & synergy)")
       }
      
     wb = openxlsx::createWorkbook();
     openxlsx::addWorksheet(wb = wb, sheetName = "Full matrix"); openxlsx::addWorksheet(wb = wb, sheetName = "Most synergistic area"); 
     
     file_[["Dose-response matrix"]] = ""; file_[["Synergy plot"]] = ""; file_[["Dose-response matrix (Control)"]] = "";
     file_ = file_[order(file_$`STE score (selective efficacy & synergy)`, decreasing = !0), ]
     
     openxlsx::writeDataTable(wb, sheet = "Full matrix", file_, colNames =T)
     invisible(lapply(1:nrow(file_), function (i) {
       img_ <- paste0("./plots/",file_$Drug1[i],"_inh_", file_$Drug2[i], ".png");  
       openxlsx::insertImage(wb, sheet = "Full matrix", file = file.path(img_),
                             width = 700, height = 700, startRow = i+1, startCol = 8, units = "px", dpi = 360)}))
     invisible(lapply(1:nrow(file_), function (i) {
       img_ <- paste0("./plots/",file_$Drug1[i],"_syn_", file_$Drug2[i], ".png");  
       openxlsx::insertImage(wb, sheet = "Full matrix", file = file.path(img_),
                             width = 700, height = 700, startRow = i+1, startCol = 9, units = "px", dpi = 360)}))
     invisible(lapply(1:nrow(file_), function (i) {
       img_ <- paste0("./plots/",file_$Drug1[i],"_tox_", file_$Drug2[i], ".png");  
       openxlsx::insertImage(wb, sheet = "Full matrix", file = file.path(img_),
                             width = 700, height = 700, startRow = i+1, startCol = 10, units = "px", dpi = 360)}))

     #change style 
     openxlsx::setColWidths(wb, sheet = "Full matrix", cols = 1:10, widths=c(12, 12, 38, 9, 9, 9, 29, 30, 30, 30), ignoreMergedCells = F)
     openxlsx::setRowHeights(wb, sheet = "Full matrix", rows = 1:nrow(file_)+1, heights = 142.5)
     style_ = openxlsx::createStyle(halign = "center", valign = "center",  wrapText = T)
     style_2 = openxlsx::createStyle( fgFill  = "#bcdfeb", halign = "center", valign = "center",  wrapText = T)
     openxlsx::addStyle(wb, sheet = "Full matrix", style_, rows = 2:(nrow(file_)+1), cols = 1:10, gridExpand = T)
     if(input$conrol_ != "SynTox"){ openxlsx::addStyle(wb, sheet = "Full matrix", style_2, rows = 2:(nrow(file_)+1), cols = 3, gridExpand = T) }
     else  { openxlsx::addStyle(wb, sheet = "Full matrix", style_2, rows = 2:(nrow(file_)+1), cols = 7, gridExpand = T) }
    
     
     
     
     file_2[["Dose-response matrix"]] = ""; file_2[["Synergy plot"]] = ""; file_2[["Dose-response matrix (Control)"]] = "";
     file_2 = file_2[order(file_2$`STE score (selective efficacy & synergy)`, decreasing = !0), ]
     
     openxlsx::writeDataTable(wb, sheet = "Most synergistic area", file_2, colNames =T)
     invisible(lapply(1:nrow(file_2), function (i) {
       img_ <- paste0("./plots/",file_2$Drug1[i],"_inh_", file_2$Drug2[i], ".png");  
       openxlsx::insertImage(wb, sheet = "Most synergistic area", file = file.path(img_),
                             width = 700, height = 700, startRow = i+1, startCol = 8, units = "px", dpi = 360)}))
     invisible(lapply(1:nrow(file_2), function (i) {
       img_ <- paste0("./plots/",file_2$Drug1[i],"_syn_", file_2$Drug2[i], ".png");  
       openxlsx::insertImage(wb, sheet = "Most synergistic area", file = file.path(img_),
                             width = 700, height = 700, startRow = i+1, startCol = 9, units = "px", dpi = 360)}))
     invisible(lapply(1:nrow(file_2), function (i) {
       img_ <- paste0("./plots/",file_2$Drug1[i],"_tox_", file_2$Drug2[i], ".png");  
       openxlsx::insertImage(wb, sheet = "Most synergistic area", file = file.path(img_),
                             width = 700, height = 700, startRow = i+1, startCol = 10, units = "px", dpi = 360)}))
     
     #change style 
     openxlsx::setColWidths(wb, sheet = "Most synergistic area", cols = 1:10, widths=c(12, 12, 38, 9, 9, 9, 29, 30, 30, 30), ignoreMergedCells = F)
     openxlsx::setRowHeights(wb, sheet = "Most synergistic area", rows = 1:nrow(file_2)+1, heights = 142.5)
     style_ = openxlsx::createStyle(halign = "center", valign = "center",  wrapText = T)
     style_2 = openxlsx::createStyle( fgFill  = "#bcdfeb", halign = "center", valign = "center",  wrapText = T)
     openxlsx::addStyle(wb, sheet = "Most synergistic area", style_, rows = 2:(nrow(file_2)+1), cols = 1:10, gridExpand = T)
     if(input$conrol_ != "SynTox"){ openxlsx::addStyle(wb, sheet = "Most synergistic area", style_2, rows = 2:(nrow(file_2)+1), cols = 3, gridExpand = T)}
     else {openxlsx::addStyle(wb, sheet = "Most synergistic area", style_2, rows = 2:(nrow(file_2)+1), cols = 7, gridExpand = T)}
     
     
     
     
     
     
     openxlsx::saveWorkbook(wb, file = "results_.xlsx", overwrite = T)
     file.copy("./results_.xlsx",file)
   },
   contentType = NULL
 )
 
 
    
  output$exp_2d2 <- downloadHandler(
    filename = function() {
      paste('figure_barplot_2D-', Sys.Date(), '.pdf', sep='')
    },
    content = function(file) {
      ggsave("fig2D.pdf", plot2Dbar(), width = 25, height = 20, units = "cm")
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
      ggsave("fig2D.pdf", plot2Dgr(), width = 25, height = 20, units = "cm")
      file.copy("./fig2D.pdf",file)
    },
    contentType = NULL
  )
  
################################################################################################################################################
  
  plot2Dbar <- function(){

    results_Bar = isolate(result_$result_); results_Bar$combi = paste0(results_Bar$Drug1, " & ", results_Bar$Drug2)

    if(input$icons2 == "selective efficacy & synergy"){
      
      if(input$rankbasedon == "Full matrix"){
        results_Bar$barscore = result_$result_$score; result_$result_ <<- results_Bar;
        g_ <- labs(title=paste0("Top ",min(20,length(results_Bar$combi))," combinations"), subtitle="sorted by summary score", caption="Ordered by summary score of selective efficacy & synergy")
      } else {
        results_Bar$barscore = result_$result_$scoreMSA; result_$result_ <<- results_Bar;
        g_ <- labs(title=paste0("Top ",min(20,length(results_Bar$combi))," combinations"), subtitle="sorted by summary score", caption="Ordered by summary score of selective efficacy & synergy calculated for Most synergistic area (3x3 dose window)")
      }
      
    } else if(input$icons2 == "efficacy & synergy") {
      
      if(input$rankbasedon == "Full matrix"){
        results_Bar$barscore = result_$result_$scoreIS;  result_$result_ <<- results_Bar;
        g_ <- labs(title=paste0("Top ",min(20,length(results_Bar$combi))," combinations"), subtitle="sorted by summary score", caption="Ordered by summary score of efficacy and synergy")
      } else {
        results_Bar$barscore = result_$result_$scoreISMSA;  result_$result_ <<- results_Bar;
        g_ <- labs(title=paste0("Top ",min(20,length(results_Bar$combi))," combinations"), subtitle="sorted by summary score", caption="Ordered by summary score of efficacy and synergy for Most synergistic area (3x3 dose window)")
        
      }
    }
    
    results_Bar = results_Bar[order(results_Bar$barscore, decreasing = T),]; results_Bar$Position = 1:nrow(results_Bar)
    results_Bar = head(results_Bar, 20)
    ggplot(data=results_Bar, aes(x=Position, y=barscore)) + theme_classic() +
      geom_bar(stat="identity", width = 0.5, fill="#990746") + scale_x_reverse(breaks=1:min(20,length(results_Bar$combi)), labels = results_Bar$combi) + 
      ylim(0,1.05) + theme(
        panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), 
        legend.box.background = element_rect(fill = "transparent")
      ) + coord_flip() +
      theme(axis.text.x = element_text(margin=ggplot2::margin(0,0,0,0,"pt"),hjust=1)) + xlab("Drug combination") + ylab("Summary score") + g_
  }
  
  plot2Dgr <- function(){
    if(input$rankbasedon2 == "Full matrix"){
      result_$result_$STE_score = result_$result_$score; 
      if(input$icons == "selective efficacy & synergy"){ 
        result_$result_$synScores_xaxis = result_$result_$synScoresReal; result_$result_$synScores_yaxis = result_$result_$s_inhReal; 
      } else {
        result_$result_$synScores_xaxis = result_$result_$synScoresReal; result_$result_$synScores_yaxis = result_$result_$inhScoresReal;
      }
    } else {
      result_$result_$STE_score = result_$result_$scoreMSA;
      if(input$icons == "selective efficacy & synergy"){ 
        result_$result_$synScores_xaxis = result_$result_$synScoresRealMSA; result_$result_$synScores_yaxis = result_$result_$s_inhRealMSA;
      } else {
        result_$result_$synScores_xaxis = result_$result_$synScoresRealMSA; result_$result_$synScores_yaxis = result_$result_$inhScoresRealMSA;
      }
    }
    if(input$icons == "selective efficacy & synergy"){
      g_ = ggplot(result_$result_, aes(x=synScores_xaxis, y=synScores_yaxis)) +
        geom_point(alpha = 0.95, aes(size=STE_score, color=STE_score), position=position_jitter(0.01)) +
        theme_classic() + scale_colour_gradientn(colours = colorRampPalette(c("#F7FBFF", "#990746"))(25)[4:22]) +
        geom_vline(xintercept = 0, linetype="dotted") +  guides(color= guide_legend(reverse = !0), size=guide_legend(reverse = !0)) +
        theme(
          panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), 
          legend.box.background = element_rect(fill = "transparent")
        )  
      if(input$rankbasedon2 == "Full matrix"){
        g_ + xlab(paste0("Synergy score (",input$state,")")) + ylab("Selective efficacy volume score") + labs(title="Top 20 combinations", subtitle="selective efficacy & synergy", 
                                                                                                              caption="Ordered by summary score of selective efficacy & synergy") 
      } else {
        g_ + xlab(paste0("Most synergistic area score (",input$state,")")) + ylab("Selective efficacy volume score calculated for most synergistic area") + labs(title="Top 20 combinations", subtitle="selective efficacy & synergy", 
                                                                                                                            caption="Ordered by summary score of selective efficacy & synergy calculated for Most synergistic area (3x3 dose window)") 
      }
        
    } else if(input$icons == "efficacy & synergy") {
      g_ = ggplot(result_$result_, aes(x=synScores_xaxis, y=synScores_yaxis)) +
        geom_point(alpha = 0.95, aes(size=STE_score, color=STE_score), position=position_jitter(0.01)) +
        theme_classic() + scale_colour_gradientn(colours = colorRampPalette(c("#F7FBFF", "#990746"))(25)[4:22]) +
        geom_vline(xintercept = 0, linetype="dotted") +  guides(color= guide_legend(reverse = !0), size=guide_legend(reverse = !0)) +
        theme(
          panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), 
          legend.box.background = element_rect(fill = "transparent")
        ) 
      if(input$rankbasedon2 == "Full matrix"){
        g_ + xlab(paste0("Synergy score (",input$state,")")) + ylab("Efficacy volume score") + labs(title="Top 20 combinations", subtitle="efficacy & synergy", 
                                                                                                   caption="Ordered by summary score of efficacy & synergy")
      } else {
        g_ + xlab(paste0("Most synergistic area score (",input$state,")")) + ylab("Efficacy volume score calculated for most synergistic area") + labs(title="Top 20 combinations", subtitle="efficacy & synergy", 
                                                                                                                 caption="Ordered by summary score of efficacy & synergy calculated for Most synergistic area (3x3 dose window)")
      }
    }
  }
  
  ########
  # first scatter plot  
  output$scatterplot <- renderPlot({
    list(input$icons, input$rankbasedon2)
    if(!is.null(isolate(result_$draw_))){

      plot2Dgr()
 
    }
  }, bg="transparent", height = 500)
  
  output$hover_info <- renderUI({

    if(!is.null(result_$draw_)){

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
    input$rankbasedon
    
    if(!is.null(result_$draw_)){ 
      plot2Dbar();
    }
  }, bg="transparent", height = 500)

  
  output$hover_info2 <- renderUI({

    if(!is.null(result_$draw_)){

      hover <- input$plot_hover2
      #point <- nearPoints(result_$result_, hover, threshold = 10, maxpoints = 1, addDist = !1)
      if(is.null(hover$y) || abs(hover$y) > 20 || abs(hover$y) < 0) return(NULL)
      ind_ = round(hover$y); order_ = order(result_$result_$barscore, decreasing = T); point = result_$result_[order_,][ind_,]

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
    
    if(input$rankbasedon2 == "Full matrix"){
    
      col_ = result_$result_$score;
      text_ = paste0("<b>STE score: ", round(col_,2),"</b><br>Selective efficacy rank: <b>", round(result_$result_$s_inhScores,2), "</b><br>Synergy rank: <b>",  round(result_$result_$synScores,2), "</b><br>Efficacy rank: <b>", 
                     round(result_$result_$inhScores,2), "</b><br>Synergy volume score (",input$state,"): <b>",  round(result_$result_$synScoresReal,2), 
                     "</b><br>Efficacy volume score(",input$conrol_,"): <b>",  round(result_$result_$inhScoresReal,2),  "</b><br>Toxicity volume score: <b>",  round(result_$result_$toxScoresReal,2), "</b>");
      title_ = "Combined Synergy, Toxicity and Efficacy landscape";
      plot_ly(x = result_$result_$synScoresReal, y = result_$result_$inhScoresReal, z = result_$result_$toxScoresReal, type = "scatter3d",  color = col_,
              text = text_, hoverinfo = 'text') %>%  
        layout(plot_bgcolor='rgba(254, 247, 234,0)', paper_bgcolor='rgba(254, 247, 234,0)', title = title_,
               scene = list(xaxis = list(title = "Synergy volume score", gridcolor = "#696969"), 
                            yaxis = list(title = "Efficacy volume score", gridcolor = "#696969"), 
                            zaxis = list(title = "Toxicity volume score", gridcolor = "#696969"),  aspectmode='cube')) 
    }
    else {
      
      col_ = result_$result_$scoreMSA;
      text_ = paste0("<b>STE score (Most synergistic area): ", round(col_,2),"</b><br>Selective efficacy rank for most synergistic area: <b>", round(result_$result_$s_inhScoresMSA,2), "</b><br>Synergy rank for most synergistic area: <b>",  round(result_$result_$synScoresMSA,2), "</b><br>Efficacy rank for most synergistic area: <b>", 
                     round(result_$result_$inhScoresMSA,2), "</b><br>Synergy volume score (",input$state,") for most synergistic area: <b>",  round(result_$result_$synScoresRealMSA,2), 
                     "</b><br>Efficacy volume score(",input$conrol_,") for most synergistic area: <b>",  round(result_$result_$inhScoresRealMSA,2),  "</b><br>Toxicity volume score for most synergistic area: <b>",  round(result_$result_$toxScoresRealMSA,2), "</b>");
      title_ = "Combined Synergy, Toxicity and Efficacy landscape calculated for most synergistic area";
      plot_ly(x = result_$result_$synScoresRealMSA, y = result_$result_$inhScoresRealMSA, z = result_$result_$toxScoresRealMSA, type = "scatter3d",  color = col_,
              text = text_, hoverinfo = 'text') %>%  
        layout(plot_bgcolor='rgba(254, 247, 234,0)', paper_bgcolor='rgba(254, 247, 234,0)', title = title_,
               scene = list(xaxis = list(title = "Most synergistic area volume score", gridcolor = "#696969"), 
                            yaxis = list(title = "Efficacy volume score", gridcolor = "#696969"), 
                            zaxis = list(title = "Toxicity volume score", gridcolor = "#696969"),  aspectmode='cube')) 
      
    }
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
  
  
  observeEvent(input$rankbasedon2, {
    updateSelectInput(session = session, "rankbasedon", choices = list("Full matrix", "Most synergistic area (3x3 dose window)"),
                      selected = input$rankbasedon2);
  })
  observeEvent(input$rankbasedon, {
    updateSelectInput(session = session, "rankbasedon2", choices = list("Full matrix", "Most synergistic area (3x3 dose window)"),
                      selected = input$rankbasedon);
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
    
    rect <- data.frame(xmin = byx*(dt_$r_-1)+1, xmax = byx*((dt_$r_-1)+2)+1, ymin = byy*(dt_$c_-1)+1, ymax = byy*((dt_$c_-1)+2)+1)
    syn.plot <- syn.plot + 
      geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="grey90", alpha=0.2, inherit.aes = !1)
    

    list(pl = syn.plot, r_=dt_$r_, c_=dt_$c_)
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
    for(i in 2:(ncol(scores.dose)-2)){
      for(j in 2:(nrow(scores.dose)-2)){
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




