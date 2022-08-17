# server functions' file

#source("scripts/read_data.R")
# source("scripts/compat_fxs.R")
# source("scripts/PT_fxs.R")
# source("scripts/ET_fxs.R")
# source("scripts/Lima_fxs.R")
# source("scripts/UK_fxs.R")

library(DT)
library(tidyverse)
library(openxlsx)
library(gtsummary)

library(histoc)

#############################################################
# function for multiple donors
res_mult <- function(df.donors = donors,
                     df.candidates = candidates,
                     df.abs = abs.d,
                     algorithm = pts,
                     n = 0,
                     check.validity = FALSE,
                     ...){
  
all_pairs <- donor_recipient_pairs(df.donors = df.donors,
                                   df.candidates = df.candidates,
                                   df.abs = df.abs,
                                   algorithm = algorithm,
                                   n = 0,
                                   check.validity = FALSE,
                                   ...)

# for loop to select 2 available candidates for a pool of donors
used.candidates <- NULL
result <- NULL
for(i in 1: length(all_pairs)){
  tmp <- all_pairs[[i]][!ID %in% used.candidates][1:2,]
  
  result <- data.table::rbindlist(list(result, tmp))
  used.candidates <- c(used.candidates, tmp$ID)
}

result[!is.na(ID),] %>% 
  rowwise() %>% 
  mutate(txScore = histoc::txscore(recipient.age = age
                                   , recipient.dialysis = dialysis
                                   , donor.age = donor_age
                                   , mmHLA_A = mmA
                                   , mmHLA_B = mmB
                                   , mmHLA_DR = mmDR)$prob5y
  ) %>% ungroup()

}
############################################

function(input, output, session) {
  
  output$ex.cands<- renderDataTable({
    
    datatable(histoc::candidates,
              rownames = FALSE)
  })
  output$ex.abs<- renderDataTable({
    datatable(histoc::cabs, #ex.abs
              rownames = FALSE)
  })
  output$ex.donors<- renderDataTable({
    datatable(histoc::donors, # ex.donors
              rownames = FALSE)
  })
  
  # reactive to the input UK files
  ukf<-reactive(input$ukfiles)
  
  # read candidates' file
  datasetCands <- reactive({
    
    file_cands <- input$file_cand
    
    if (is.null(file_cands))
      return(NULL)
    
    if (input$fileSepDF == 1) {
      data<-read.csv(file_cands$datapath)
    } else if (input$fileSepDF == 2) {
      read.delim(file_cands$datapath)
    } else if (input$fileSepDF == 3) {
      data<-read.csv2(file_cands$datapath)
    } else {data<-read.table(file_cands$datapath)}
    
     
    validate(
      if(ukf() == 1){need(identical(colnames(data),c("ID","bg","A1","A2","B1","B2","DR1","DR2","age","dialysis","cPRA", "Tier", "MS", "RRI")), 
                             "Candidates column names are not the necessary for UK algorithm!")} else {need(identical(colnames(data),c("ID","bg","A1","A2","B1","B2","DR1","DR2","age","dialysis","cPRA")), 
           "Candidates column names are not identical to example data!")}
      )
    
    data %>% 
      mutate_at(vars(ID, A1,A2,B1,B2,DR1,DR2),as.character) 
    
    
  })
  
  # read donors' file 
  datasetDonors <- reactive({
    
    file_donors <- input$file_donor
    
    if (is.null(file_donors))
      return(NULL)
    
    if (input$fileSepDF == 1) {
      data<-read.csv(file_donors$datapath)
    } else if (input$fileSepDF == 2) {
      read.delim(file_donors$datapath)
    } else if (input$fileSepDF == 3) {
      data<-read.csv2(file_donors$datapath)
    } else {data<-read.table(file_donors$datapath)}
    
    
    validate(
      if(ukf() == 1){
        need(identical(colnames(data),c(colnames(histoc::donors),"DRI")), 
             "Donors column names are not the necessary for UK algorithm!")
      } else {
      need(identical(colnames(data),colnames(histoc::donors)), 
           "Donors column names are not identical to example data!")}
    )
    
    data %>% 
      mutate_at(vars(ID,A1,A2,B1,B2,DR1,DR2),as.character) %>% 
      mutate_at(vars(age), as.numeric)
    
  })
  
  # read candidates' antibodies' file
  datasetAbs <- reactive({
    
    file_abss <- input$file_abs
    
    if (is.null(file_abss))
      return(NULL)
    
    if (input$fileSepDF == 1) {
      data<-read.csv(file_abss$datapath)
    } else if (input$fileSepDF == 2) {
      read.delim(file_abss$datapath)
    } else if (input$fileSepDF == 3) {
      data<-read.csv2(file_abss$datapath)
    } else {data<-read.table(file_abss$datapath)}
    
    
    validate(
      need(identical(colnames(data),colnames(histoc::cabs)), 
           "HLA antibodies column names are not identical to example data!")
    )
    
    data %>% 
      mutate_at(vars(ID), as.character)
    
  })
  
  output$sel.cands<- renderDataTable({
    datasetCands() %>% datatable(rownames = FALSE)
  })
  
  output$sel.abs<- renderDataTable({
    datasetAbs() %>% datatable(rownames = FALSE)
  })
  
  output$sel.donors<- renderDataTable({
    datasetDonors() %>% datatable(rownames = FALSE)
  })
  

  ############################
  ### PT algorithm
  ############################
 
  ## for 10 best candidates selected according to unique donor 
  
  output$res1 <- renderDataTable({
    
    if (input$dataInput == 1) {candidates<-histoc::candidates} else {candidates<-datasetCands()}
    if (input$dataInput == 1) {abs.d<-histoc::cabs} else {abs.d<-datasetAbs()}
    
    validate(
      need(candidates != "", "Please select a candidates data set!")
    )
    
    validate(
      need(abs.d != "", "Please select candidates' HLA antibodies data set!")
    )
    
      dt<-histoc::pts(iso = input$iso, # isogroup compatibility
                  dABO = input$dabo, # donor's blood group
                  dA = c(input$a1,input$a2),
                  dB = c(input$b1,input$b2),
                  dDR = c(input$dr1,input$dr2),
                  donor.age = input$dage, # donor's age
                  data = candidates, # data file with candidates
                  points.80 = as.numeric(input$pra8), # points for a PRA equal or higher than 80%
                  points.50 = as.numeric(input$pra5), # points for a PRA equal or higher than 50%
                  points.dialysis = input$dialysis, # points for each month on dialysis
                  points.age = input$age_dif, # points for age difference in PT punctuation table
                  itemA = as.numeric(input$a), # points for A) on PT points table
                  itemB = as.numeric(input$b), # points for B) on PT points table
                  itemC = as.numeric(input$c), # points for C) on PT points table
                  itemD = as.numeric(input$d), # points for D) on PT points table
                  itemE = as.numeric(input$e), # points for E) on PT points table
                  df.abs = abs.d, # candidates' HLA antibodies
                  n = 10,
                  check.validity = FALSE)
      
      dt <- dt %>%
        rowwise() %>% 
        mutate(txScore = histoc::txscore(recipient.age = age
                                 , recipient.dialysis = dialysis
                                 , donor.age = donor_age
                                 , mmHLA_A = mmA
                                 , mmHLA_B = mmB
                                 , mmHLA_DR = mmDR)$prob5y
               ) %>% ungroup()

    datatable(dt, options = list(pageLength = 5, dom = 'tip'))
  })

  
  #### to reset PT sidebarpanel
  observeEvent(input$reset_inputPT, {
    shinyjs::reset("side-panelPT")
  })
  
  N <- 2
  
  compute_resm <- reactiveVal()
 
  observeEvent(input$Go, {
    
    compute_resm(NULL)
    
    withProgress(message = 'Calculation in progress, be patient!', {
      for(i in 1:N){
        # Long Running Task
        Sys.sleep(1)
        # Update progress
        incProgress(1/N)
        }

    if (input$dataInput == 1) {candidates<-histoc::candidates} else {candidates<-datasetCands() %>% 
      select(ID, bg,A1,A2,B1,B2,DR1,DR2,age,dialysis,cPRA)}
    if (input$dataInput == 1) {abs.d<-histoc::cabs} else {abs.d<-datasetAbs()}
    if (input$dataInput == 1) {donors<-histoc::donors} else {donors<-datasetDonors()}
    
    validate(
      need(candidates != "", "Please select a candidates data set!")
    )
    
    validate(
      need(abs.d != "", "Please select candidates' HLA antibodies data set!")
    )
    
    validate(
      need(donors != "", "Please select donors' data set!")
    )
    
    dt <- res_mult(df.donors = donors,
                   df.candidates = candidates,
                   df.abs = abs.d,
                   algorithm = pts,
                   n = 0,
                   check.validity = FALSE,
                   iso = input$iso,
                   points.80 = as.numeric(input$pra8),
                   points.50 = as.numeric(input$pra5),
                   points.dialysis = input$dialysis,
                   points.age = input$age_dif,
                   itemA = as.numeric(input$a),
                   itemB = as.numeric(input$b),
                   itemC = as.numeric(input$c),
                   itemD = as.numeric(input$d),
                   itemE = as.numeric(input$e)
                   )

    compute_resm(dt)
    })
    
    
  })

  output$resm <- renderDataTable({
    compute_resm()
  })

  # Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("PT_results", ".csv", sep = "")
    },
    content = function(file) {
      write.csv2(compute_resm(), file, row.names = FALSE, fileEncoding="latin1")
    }
  )
  
  ## Resume dataset results from PT algorithm
  output$resumePT <-
    render_gt({
      
      validate(
        need(compute_resm() != "", "Results will be presented after the run!")
      )
      
      tabsum<-compute_resm() %>% 
        select(bg, age, dialysis, cPRA, HI, mmHLA, txScore) %>% 
        rename(`Blood group` = bg,
               `receptores' age (years)` = age,
               `time on dialysis (months)` = dialysis,
               `Hiper Immunized` = HI,
               `HLA miss matchs` = mmHLA,
               TxScore = txScore)
      
      tbl_summary(tabsum) %>% as_gt()
    })

  
  ############################
  ### ET algorithm
  ############################

  ## for 10 first candidates selected according to unique donor 
  output$res1ET <- renderDataTable({
    
    if (input$hlafreqs == 1) {
      hlaA <- histoc::hlaApt
      hlaB <- histoc::hlaBpt
      hlaDR <- histoc::hlaDRpt
      ABOfreq <- histoc::ABOpt
      } else {
        hlaA <- histoc::hlaAet
        hlaB <- histoc::hlaBet
        hlaDR <- histoc::hlaDRet
        ABOfreq <- histoc::ABOpt
        }

    if (input$dataInput == 1) {candidates<-histoc::candidates} else {candidates<-datasetCands()}
    if (input$dataInput == 1) {abs.d<-histoc::cabs} else {abs.d<-datasetAbs()}

    validate(
      need(candidates != "", "Please select a candidates data set!")
    )

    validate(
      need(abs.d != "", "Please select candidates' HLA antibodies data set!")
    )

    dt<-et(iso = input$isoET, 
           dABO = input$daboET, 
           dA = c(input$a1ET,input$a2ET),
           dB = c(input$b1ET,input$b2ET),
           dDR = c(input$dr1ET,input$dr2ET),
           donor.age = input$dageET,
           data = candidates,
           month = as.numeric(input$tdET),
           mm0 = as.numeric(input$mm0),
           mm1 = as.numeric(input$mm1),
           mm2 = as.numeric(input$mm2),
           mm3 = as.numeric(input$mm3),
           mm4 = as.numeric(input$mm4),
           mm5 = as.numeric(input$mm5),
           mm6 = as.numeric(input$mm6),
           df.abs = abs.d,
           hlaA = hlaA,
           hlaB = hlaB,
           hlaDR = hlaDR,
           abo.freq = ABOfreq,
           n = 10)

    dt <- dt %>%
      rowwise() %>% 
      mutate(txScore = histoc::txscore(recipient.age = age
                                       , recipient.dialysis = dialysis
                                       , donor.age = donor_age
                                       , mmHLA_A = mmA
                                       , mmHLA_B = mmB
                                       , mmHLA_DR = mmDR)$prob5y
      ) %>% ungroup()
    
    datatable(dt, options = list(pageLength = 5, dom = 'tip'))
  }) 
  
  observeEvent(input$reset_inputET, {
    shinyjs::reset("side-panelET")
  })
 
  ## compute nultiple results for ET algorithm
  compute_resmET <- reactiveVal()
  
  observeEvent(input$GoET, {
    
    compute_resmET(NULL)

    withProgress(message = 'Calculation in progress, be patient!', {
      for(i in 1:N){
        # Long Running Task
        Sys.sleep(0.01)
        # Update progress
        incProgress(1/N)
      }

      if (input$hlafreqs == 1) {
        hlaA <- histoc::hlaApt
        hlaB <- histoc::hlaBpt
        hlaDR <- histoc::hlaDRpt
        ABOfreq <- histoc::ABOpt
      } else {
        hlaA <- histoc::hlaAet
        hlaB <- histoc::hlaBet
        hlaDR <- histoc::hlaDRet
        ABOfreq <- histoc::ABOpt
      }
     
      if (input$dataInput == 1) {candidates<-histoc::candidates} else {candidates<-datasetCands()}
      if (input$dataInput == 1) {abs.d<-histoc::cabs} else {abs.d<-datasetAbs()}
      if (input$dataInput == 1) {donors<-histoc::donors} else {donors<-datasetDonors()}
      
      validate(
        need(candidates != "", "Please select a candidates data set!")
      )
      
      validate(
        need(abs.d != "", "Please select candidates' HLA antibodies data set!")
      )
      
      validate(
        need(donors != "", "Please select donors' data set!")
      )
      
      dt <- res_mult(df.donors = donors,
                     df.candidates = candidates,
                     df.abs = abs.d,
                     algorithm = et,
                     n = 0,
                     check.validity = FALSE,
                     iso = input$isoET,
                     month = as.numeric(input$tdET),
                     mm0 = as.numeric(input$mm0),
                     mm1 = as.numeric(input$mm1),
                     mm2 = as.numeric(input$mm2),
                     mm3 = as.numeric(input$mm3),
                     mm4 = as.numeric(input$mm4),
                     mm5 = as.numeric(input$mm5),
                     mm6 = as.numeric(input$mm6),
                     hlaA = hlaA,
                     hlaB = hlaB,
                     hlaDR = hlaDR,
                     abo.freq = ABOfreq,
                     )

      compute_resmET(dt)
      
      })
    
  })
  
  output$resmET <- renderDataTable({
    compute_resmET()
  })
  
  # Downloadable csv of selected dataset ----
  output$downloadDataET <- downloadHandler(
    filename = function() {
      paste("ET_results", ".csv", sep = "")
    },
    content = function(file) {
      write.csv2(compute_resmET(), file, row.names = FALSE, fileEncoding="latin1")
    }
  )
  
  
  ## Resume dataset results from ET algorithm
  output$resumeET <-
    render_gt({
      
      validate(
        need(compute_resmET() != "", "Results will be presented after the run!")
      )
      
      tabsum<-compute_resmET() %>% 
        select(bg, age, dialysis, cPRA, HI, mmHLA, txScore) %>% 
        rename(`Blood group` = bg,
               `receptores' age (years)` = age,
               `time on dialysis (months)` = dialysis,
               `Hiper Immunized` = HI,
               `HLA miss matchs` = mmHLA,
               TxScore = txScore)
      
      tbl_summary(tabsum) %>% as_gt()
    })
  
  
  ############################
  ### Lima algorithm
  ############################
  
  #### to reset LIMA sidebarpanel
  observeEvent(input$reset_inputLIMA, {
    shinyjs::reset("side-panelLima")
  })
  
  
  output$res1LIMA <- renderDataTable({
    
    if (input$dataInput == 1) {candidates<-histoc::candidates} else {candidates<-datasetCands()}
    if (input$dataInput == 1) {abs.d<-histoc::cabs} else {abs.d<-datasetAbs()}
    
    validate(
      need(candidates != "", "Please select a candidates data set!")
    )
    
    validate(
      need(abs.d != "", "Please select candidates' HLA antibodies data set!")
    )

    dt<-histoc::lima(iso = input$isoLIMA, # isogroup compatibility
                   dABO = input$daboLIMA, # donor's blood group
                   dA = c(input$a1LIMA,input$a2LIMA),
                   dB = c(input$b1LIMA,input$b2LIMA),
                   dDR = c(input$dr1LIMA,input$dr2LIMA),
                   donor.age = input$dageLIMA, # donor's age
                   data = candidates, # data file with candidates
                   df.abs = abs.d, # candidates' HLA antibodies
                   n = 10,
                   q2 = input$td2q,
                   q3 = input$td3q,
                   cPRA1 = 50,
                   cPRA2 = 85,
                   check.validity = FALSE)
    
    dt <- dt %>%
      rowwise() %>% 
      mutate(txScore = histoc::txscore(recipient.age = age
                                       , recipient.dialysis = dialysis
                                       , donor.age = donor_age
                                       , mmHLA_A = mmA
                                       , mmHLA_B = mmB
                                       , mmHLA_DR = mmDR)$prob5y
      ) %>% ungroup()
    
    datatable(dt, options = list(pageLength = 5, dom = 'tip'))
    
    })
  
  ## compute multiple results for LIMA algorithm
  compute_resmLIMA <- reactiveVal()
  
  observeEvent(input$GoLIMA, {
    
    compute_resmLIMA(NULL)
    
    withProgress(message = 'Calculation in progress, be patient!', {
      for(i in 1:N){
        # Long Running Task
        Sys.sleep(1)
        # Update progress
        incProgress(1/N)
      }
      
      if (input$dataInput == 1) {candidates<-histoc::candidates} else {candidates<-datasetCands()}
      if (input$dataInput == 1) {abs.d<-histoc::cabs} else {abs.d<-datasetAbs()}
      if (input$dataInput == 1) {donors<-histoc::donors} else {donors<-datasetDonors()}
      
      validate(
        need(candidates != "", "Please select a candidates data set!")
      )
      
      validate(
        need(abs.d != "", "Please select candidates' HLA antibodies data set!")
      )
      
      validate(
        need(donors != "", "Please select donors' data set!")
      )
      
      dt <- res_mult(df.donors = donors,
                     df.candidates = candidates,
                     df.abs = abs.d,
                     algorithm = lima,
                     n = 0,
                     check.validity = FALSE,
                     iso = input$isoLIMA,
                     q2 = input$td2q,
                     q3 = input$td3q,
                     cPRA1 = 50,
                     cPRA2 = 85
                     )
 
      compute_resmLIMA(dt)

    })
    })
  
  output$resmLIMA <- renderDataTable({
    compute_resmLIMA()
  })
  
  # Downloadable csv of selected dataset ----
  output$downloadDataLIMA <- downloadHandler(
    filename = function() {
      paste("LIMA_results", ".csv", sep = "")
    },
    content = function(file) {
      write.csv2(compute_resmLIMA(), file, row.names = FALSE, fileEncoding="latin1")
    }
  )
  
  ## Resume dataset results from LIMA algorithm
  output$resumeLIMA <-
    render_gt({
      
      validate(
        need(compute_resmLIMA() != "", "Results will be presented after the run!")
      )
      
      tabsum<-compute_resmLIMA() %>% 
        select(bg, age, dialysis, cPRA, HI, mmHLA, txScore) %>% 
        rename(`Blood group` = bg,
               `receptores' age (years)` = age,
               `time on dialysis (months)` = dialysis,
               `Hiper Immunized` = HI,
               `HLA miss matchs` = mmHLA,
               TxScore = txScore)
      
      tbl_summary(tabsum) %>% as_gt()
    })
  
  ############################
  ### UK Transplant
  ############################

  # Donor-recipient risk index combinations
  # multiplyer<-reactive({input$multipleUK})
  
  output$tableDRriskUK <- renderDT({
    dt<-data.frame(R1=c(1000,700,350,0),
                   R2=c(700,1000,500,350),
                   R3=c(350,500,1000,700),
                   R4=c(0,350,700,1000))
    rownames(dt)<-c("D1","D2","D3","D4")
    
    datatable(dt * 1, # multiplyer()
              selection = 'single', escape=FALSE,
              options = list(searching = FALSE, dom = 't'))
  })
  
  # Matchability - illustrative plot
  pointsM<-reactive({
    input$mUK * (1 + ((1:9) / input$nUK)^input$oUK)
  })
  output$matchability<-renderPlot({
    ggplot(data.frame(Points = pointsM(), MatchScore = 1:9)) +
      geom_line(aes(MatchScore,Points)) + 
      ggtitle("Points scored illustration")
  })
  
  #### to reset UK sidebarpanel
  observeEvent(input$reset_inputUK, {
    shinyjs::reset("side-panelUK")
  })
  
  
  ############################ UK algorithm ###################

  ## compute DRI for one donor
  driv<-reactive({
    exp(0.023 * (input$dageUK-50) +
          -0.152 * ((input$dheightUK - 170) / 10) +
          0.149 * ifelse(input$dhtUK == 'Yes', 1, 0) +
          -0.184 * ifelse(input$dsexUK == 'Female', 1, 0) +
          0.19 * ifelse(input$dcmvUK == 'Yes', 1, 0) +
          -0.023 * (input$dgfrUK-90)/10 +
          0.015 * input$dhospUK
    )

  })

  output$dri<-renderText({
   paste("DRI is:", round(driv(),2), "; the donor belong to",
         ifelse(driv() <= 0.79, "D1",
                ifelse(driv() <= 1.12,"D2",
                       ifelse(driv() <= 1.5,"D3","D4"))))
    })


  output$res1UK <- renderDataTable({

    if (input$dataInput == 1) {candidates<-histoc::candidates.uk} else {candidates<-datasetCands()}
    if (input$dataInput == 1) {abs.d<-histoc::cabs} else {abs.d<-datasetAbs()}

    validate(
      need(candidates != "", "Please select a candidates data set!")
    )

    validate(
      need(abs.d != "", "Please select candidates' HLA antibodies data set!")
    )

    dt<-histoc::uk(DRI = ifelse(driv() <= 0.79, "D1",
                               ifelse(driv() <= 1.12,"D2",
                                      ifelse(driv() <= 1.5,"D3","D4"))), # Donor RisK Index group
                  dA = c(input$a1UK,input$a2UK), # donor's HLA typing'
                  dB = c(input$b1UK,input$b2UK),
                  dDR = c(input$dr1UK,input$dr2UK),
                  dABO = input$daboUK, # donors' blood group
                  donor.age = input$dageUK, # donors' age
                  data = candidates, # data file with candidates
                  D1R1 = 1000,
                  D1R2 = 700,
                  D1R3 = 350,
                  D1R4 = 0,
                  D2R1 = 700,
                  D2R2 = 1000,
                  D2R3 = 500,
                  D2R4 = 350,
                  D3R1 = 350,
                  D3R2 = 500,
                  D3R3 = 1000,
                  D3R4 = 700,
                  D4R1 = 0,
                  D4R2 = 350,
                  D4R3 = 700,
                  D4R4 = 1000,
                  ptsDial = input$tdUK,
                  a1 = input$aa1UK, 
                  a2 = input$aa2UK, 
                  b1 = input$bb1UK, 
                  b2 = input$bb2UK, 
                  b3 = input$bb3UK, 
                  m = input$mUK, 
                  nn = input$nUK, 
                  o = input$oUK, 
                  mm1 = as.numeric(input$mm1UK), 
                  mm23 = as.numeric(input$mm23UK), 
                  mm46 = as.numeric(input$mm46UK), 
                  pts = input$bloodUK, 
                  df.abs = abs.d, 
                  n = 10,
                  check.validity = FALSE
                  )

    dt <- dt %>%
      rowwise() %>% 
      mutate(txScore = histoc::txscore(recipient.age = age
                                       , recipient.dialysis = dialysis
                                       , donor.age = donor_age
                                       , mmHLA_A = mmA
                                       , mmHLA_B = mmB
                                       , mmHLA_DR = mmDR)$prob5y
      ) %>% ungroup()

    datatable(dt, options = list(pageLength = 5, dom = 'tip'))

  })

  ## compute nultiple results for UK algorithm
  compute_resmUK <- reactiveVal()

  observeEvent(input$GoUK, {

    compute_resmUK(NULL)

    withProgress(message = 'Calculation in progress, be patient!', {
      for(i in 1:N){
        # Long Running Task
        Sys.sleep(1)
        # Update progress
        incProgress(1/N)
      }


      if (input$dataInput == 1) {candidates<-histoc::candidates.uk} else {candidates<-datasetCands()}
      if (input$dataInput == 1) {abs.d<-histoc::cabs} else {abs.d<-datasetAbs()}
      if (input$dataInput == 1) {donors<-histoc::donors.uk %>%
        select(ID, bg, A1, A2, B1, B2, DR1, DR2, age, DRI)} else {donors<-datasetDonors()}


      validate(
        need(candidates != "", "Please select a candidates data set!")
      )

      validate(
        need(abs.d != "", "Please select candidates' HLA antibodies data set!")
      )

      validate(
        need(donors != "", "Please select donors' data set!")
      )
      
      dt <- res_mult(df.donors = donors,
                     df.candidates = candidates,
                     df.abs = abs.d,
                     algorithm = uk,
                     n = 0,
                     check.validity = FALSE,
                     DRI = ifelse(driv() <= 0.79, "D1",
                                  ifelse(driv() <= 1.12,"D2",
                                         ifelse(driv() <= 1.5,"D3","D4"))),
                     D1R1 = 1000,
                     D1R2 = 700,
                     D1R3 = 350,
                     D1R4 = 0,
                     D2R1 = 700,
                     D2R2 = 1000,
                     D2R3 = 500,
                     D2R4 = 350,
                     D3R1 = 350,
                     D3R2 = 500,
                     D3R3 = 1000,
                     D3R4 = 700,
                     D4R1 = 0,
                     D4R2 = 350,
                     D4R3 = 700,
                     D4R4 = 1000,
                     ptsDial = input$tdUK,
                     a1 = input$aa1UK, 
                     a2 = input$aa2UK, 
                     b1 = input$bb1UK, 
                     b2 = input$bb2UK, 
                     b3 = input$bb3UK, 
                     m = input$mUK, 
                     nn = input$nUK, 
                     o = input$oUK, 
                     mm1 = as.numeric(input$mm1UK), 
                     mm23 = as.numeric(input$mm23UK), 
                     mm46 = as.numeric(input$mm46UK), 
                     pts = input$bloodUK
                     )

      # add to reactiveval
      compute_resmUK(dt)

    })
  })

  output$resmUK <- renderDataTable({
    compute_resmUK()
  })

  # Downloadable csv of selected dataset ----
  output$downloadDataUK <- downloadHandler(
    filename = function() {
      paste("UK_results", ".csv", sep = "")
    },
    content = function(file) {
      write.csv2(compute_resmUK(), file, row.names = FALSE, fileEncoding="latin1")
    }
  )

  ## Resume dataset results from UK algorithm
  output$resumeUK <-
    render_gt({

      validate(
        need(compute_resmUK() != "", "Results will be presented after the run!")
      )

      tabsum<-compute_resmUK() %>%
        select(bg, age, dialysis, cPRA, Tier, mmHLA, txScore) %>%
        rename(`Blood group` = bg,
               `receptores' age (years)` = age,
               `time on dialysis (months)` = dialysis,
               `Tier` = Tier,
               `HLA miss matchs` = mmHLA,
               TxScore = txScore)

      tbl_summary(tabsum) %>% as_gt()
    })
  
#   ######################################################################################
#   ############ ligação entre inputs ############
#   ############   
#   ############ donor's age
#   observeEvent(input$dage,{
#     new <- input$dage
#     updateSliderInput(session, "dageET", value = new) 
#   })
#   
#   observeEvent(input$dageET,{ 
#     new <- input$dageET
#     updateSliderInput(session, "dage", value = new) 
#   })
#   
#   observeEvent(input$dage,{
#     new <- input$dage
#     updateSliderInput(session, "dageLIMA", value = new) 
#   })
#   
#   observeEvent(input$dageLIMA,{ 
#     new <- input$dageLIMA
#     updateSliderInput(session, "dage", value = new) 
#   })
#   
#   observeEvent(input$dage,{
#     new <- input$dage
#     updateSliderInput(session, "dageUK", value = new) 
#   })
#   
#   observeEvent(input$dageUK,{ 
#     new <- input$dageUK
#     updateSliderInput(session, "dage", value = new) 
#   })
#   
#   ############ donor's ABO
#   observeEvent(input$dabo,{
#     new <- input$dabo
#     updateRadioButtons(session, "daboET", selected = new)
#   })
#   
#   observeEvent(input$daboET,{
#     new <- input$daboET
#     updateRadioButtons(session, "dabo", selected = new)
#   })
#   
#   observeEvent(input$dabo,{
#     new <- input$dabo
#     updateRadioButtons(session, "daboLIMA", selected = new)
#   })
#   
#   observeEvent(input$daboLIMA,{
#     new <- input$daboLIMA
#     updateRadioButtons(session, "dabo", selected = new)
#   })
#   
#   observeEvent(input$dabo,{
#     new <- input$dabo
#     updateRadioButtons(session, "daboUK", selected = new)
#   })
#   
#   observeEvent(input$daboUK,{
#     new <- input$daboUK
#     updateRadioButtons(session, "dabo", selected = new)
#   })
#   
# ############ donor's typing HLA-A a1
# observeEvent(input$a1,{
#   new <- input$a1
#   updateTextAreaInput(session, "a1ET", value  = new)
# })
# 
# observeEvent(input$a1ET,{
#   new <- input$a1ET
#   updateTextAreaInput(session, "a1", value  = new)
# })
# 
# observeEvent(input$a1,{
#   new <- input$a1
#   updateTextAreaInput(session, "a1LIMA", value  = new)
# })
# 
# observeEvent(input$a1LIMA,{
#   new <- input$a1LIMA
#   updateTextAreaInput(session, "a1", value  = new)
# })
# 
# observeEvent(input$a1,{
#   new <- input$a1
#   updateTextAreaInput(session, "a1UK", value  = new)
# })
# 
# observeEvent(input$a1UK,{
#   new <- input$a1UK
#   updateTextAreaInput(session, "a1", value  = new)
# })
# 
# ############ donor's typing HLA-A a2
# observeEvent(input$a2,{
#   new <- input$a2
#   updateTextAreaInput(session, "a2ET", value  = new)
# })
# 
# observeEvent(input$a2ET,{
#   new <- input$a2ET
#   updateTextAreaInput(session, "a2", value  = new)
# })
# 
# observeEvent(input$a2,{
#   new <- input$a2
#   updateTextAreaInput(session, "a2LIMA", value  = new)
# })
# 
# observeEvent(input$a2LIMA,{
#   new <- input$a2LIMA
#   updateTextAreaInput(session, "a2", value  = new)
# })
# 
# observeEvent(input$a2,{
#   new <- input$a2
#   updateTextAreaInput(session, "a2UK", value  = new)
# })
# 
# observeEvent(input$a2UK,{
#   new <- input$a2UK
#   updateTextAreaInput(session, "a2", value  = new)
# })
# 
# ############ donor's typing HLA-B b1
# observeEvent(input$b1,{
#   new <- input$b1
#   updateTextAreaInput(session, "b1ET", value  = new)
# })
# 
# observeEvent(input$b1ET,{
#   new <- input$b1ET
#   updateTextAreaInput(session, "b1", value  = new)
# })
# 
# observeEvent(input$b1,{
#   new <- input$b1
#   updateTextAreaInput(session, "b1LIMA", value  = new)
# })
# 
# observeEvent(input$b1LIMA,{
#   new <- input$b1LIMA
#   updateTextAreaInput(session, "b1", value  = new)
# })
# 
# observeEvent(input$b1,{
#   new <- input$b1
#   updateTextAreaInput(session, "b1UK", value  = new)
# })
# 
# observeEvent(input$b1UK,{
#   new <- input$b1UK
#   updateTextAreaInput(session, "b1", value  = new)
# })
# 
# ############ donor's typing HLA-B b2
# observeEvent(input$b2,{
#   new <- input$b2
#   updateTextAreaInput(session, "b2ET", value  = new)
# })
# 
# observeEvent(input$b2ET,{
#   new <- input$b2ET
#   updateTextAreaInput(session, "b2", value  = new)
# })
# 
# observeEvent(input$b2,{
#   new <- input$b2
#   updateTextAreaInput(session, "b2LIMA", value  = new)
# })
# 
# observeEvent(input$b2LIMA,{
#   new <- input$b2LIMA
#   updateTextAreaInput(session, "b2", value  = new)
# })
# 
# observeEvent(input$b2,{
#   new <- input$b2
#   updateTextAreaInput(session, "b2UK", value  = new)
# })
# 
# observeEvent(input$b2UK,{
#   new <- input$b2UK
#   updateTextAreaInput(session, "b2", value  = new)
# })
# 
# ############ donor's typing HLA-DR dr1
# observeEvent(input$dr1,{
#   new <- input$dr1
#   updateTextAreaInput(session, "dr1ET", value  = new)
# })
# 
# observeEvent(input$dr1ET,{
#   new <- input$dr1ET
#   updateTextAreaInput(session, "dr1", value  = new)
# })
# 
# observeEvent(input$dr1,{
#   new <- input$dr1
#   updateTextAreaInput(session, "dr1LIMA", value  = new)
# })
# 
# observeEvent(input$dr1LIMA,{
#   new <- input$dr1LIMA
#   updateTextAreaInput(session, "dr1", value  = new)
# })
# 
# observeEvent(input$dr1,{
#   new <- input$dr1
#   updateTextAreaInput(session, "dr1UK", value  = new)
# })
# 
# observeEvent(input$dr1UK,{
#   new <- input$dr1UK
#   updateTextAreaInput(session, "dr1", value  = new)
# })
# 
# ############ donor's typing HLA-DR dr2
# observeEvent(input$dr2,{
#   new <- input$dr2
#   updateTextAreaInput(session, "dr2ET", value  = new)
# })
# 
# observeEvent(input$dr2ET,{
#   new <- input$dr2ET
#   updateTextAreaInput(session, "dr2", value  = new)
# })
# 
# observeEvent(input$dr2,{
#   new <- input$dr2
#   updateTextAreaInput(session, "dr2LIMA", value  = new)
# })
# 
# observeEvent(input$dr2LIMA,{
#   new <- input$dr2LIMA
#   updateTextAreaInput(session, "dr2", value  = new)
# })
# 
# observeEvent(input$dr2,{
#   new <- input$dr2
#   updateTextAreaInput(session, "dr2UK", value  = new)
# })
# 
# observeEvent(input$dr2UK,{
#   new <- input$dr2UK
#   updateTextAreaInput(session, "dr2", value  = new)
# })
  
}

