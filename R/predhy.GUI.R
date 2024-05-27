#' @title Graphical User Interface for R package predhy
#' @description Graphical User Interface for cross validation, genotype conversion and hybrid performance prediction.
#' @return No return value, called for Graphical User Interface
#' @examples
#' {
#' predhy.GUI()}
#' @export

predhy.GUI <- function(){
  if(interactive()){
    ui <- fluidPage(navbarPage(title = h3('Predhy'),
          tabPanel(h4('cv'),
          navlistPanel(widths = c(3,9),
          tabPanel(h2('Evaluate Trait Predictability via Cross Validation'),title = 'Description',
             helpText('The cv function evaluates trait predictability based on eight GS methods via k-fold cross validation. The trait predictability is defined as the squared Pearson correlation coefficient between the
                      observed and the predicted trait values.')),
          tabPanel(h2('Input genotype'),title='Genotype Input',
             fluidRow(column(width=4,fileInput('inbred_gen1',label = h4('inbred_gen'))),
                      column(width = 8,helpText('A matrix for genotypes of parental lines in numeric format, coded as 1, 0 and
                                                -1. The row.names of inbred_gen must be provied. It can be obtained from the original genotype using',code('convertgen'),'function')))),
          tabPanel(h2('Input hybrid phenotype'),title = 'Hybrid Phenotype',
             fluidRow(column(width = 4,fileInput('hybrid_phe1',label = h4('hybrid_phe'))),
                      column(width=8,helpText('A data frame with three columns. The first column and the second column are the
                                              names of male and female parents of the corresponding hybrids, respectively;the third column is the phenotypic values of hybrids. The names of male and
                                              female parents must match the rownames of inbred_gen. Missing (NA) valuesare not allowed.')))),
          tabPanel(h2('Input parent phenotype'),title='Parent Phenotype',
                   fluidRow(column(width=4,radioButtons('cvinbredpheinput',label = h4('A matrix of a phenotypic values of parent (Optional)'),
                                                        c('Not included'='NULL','Input parent phenotype'='input'))),
                            column(width=8,conditionalPanel(condition = 'input.cvinbredpheinput=="input"',
                                                            fileInput('inbred_phe',label=h4('Parent Phenotype')))))),
          tabPanel(h2('Input design matrix of the fixed effects & domiance genotypes'),title='Optional Input',
             fluidRow(column(width=4,radioButtons('cvfixinput',label = h4('design matrix of the fixed effects(Optional)'),
                                         c('Not included'='NULL','Input a design matrix'='input'))),
                      column(width=8,conditionalPanel(condition = 'input.cvfixinput=="input"',fileInput('fix',label=h4('fixed effects'))))),
             fluidRow(column(width=4,radioButtons('cvgendinput',label=h4('domiance genotypes(Optional)'),
                                         c('Not included'='NULL','Include domiance genotypes'='input'))))),
          tabPanel(h2('Select models & other parameters'),title = 'Parameters',
             fluidRow(column(width=6,selectInput('cvmethod',label=h4('method,eight GS methods'),
                                         choices = list('GBLUP'='GBLUP','BayesB'='BayesB',
                                                   'RKHS'='RKHS','PLS'='PLS',
                                                   'LASSO'='LASSO','EN'='EN',
                                                   'XGBoost'='XGBoost','LightGBM'='LightGBM','ALL'='ALL'),
                                                   selected = 'GBLUP'))),
             fluidRow(column(width=6,numericInput('nfold',label=h4('the number of folds'), value = 5))),
             fluidRow(column(width=4,numericInput('ntimes',label=h4('replicates'), value = 1)),
                                  column(width=4,numericInput('cvseed',label=h4('the random number'),value = 133))),
             fluidRow(column(width=4,numericInput('cpu',label=h4('the number of CPU'), value = 1)))),
                      tabPanel(h2('Trait predictability (R^2)'),title = 'CV Results',
             fluidRow(column(width=4,h3(textOutput('cv1')))),
             fluidRow(column(width=12,plotOutput('cvp'))))
                                    )),
          tabPanel(h4('predhy.predict'),
          navlistPanel(widths = c(3,9),
             tabPanel(h2('Predict the Performance of Hybrids'),title = 'Description',
                      helpText('Predict all potential crosses of a given set of parents using a subset of crosses as the training sample.')),
             tabPanel(h2('Input genotype & phenotype'),title = 'Input files',
                fluidRow(column(width=4,fileInput('inbred_gen',label = h4('inbred_gen'))),
                         column(width=8,helpText('A matrix for genotypes of parental lines in numeric format, coded as 1, 0 and
                                      -1. The row.names of inbred_gen must be provied. It can be obtained from the original genotype using',code('convertgen'),'function'))),
                fluidRow(column(width = 4,fileInput('hybrid_phe',label = h4('hybrid_phe'))),
                         column(width=8,helpText('A data frame with three columns. The first column and the second column are the
                                      names of male and female parents of the corresponding hybrids, respectively;the third column is the phenotypic values of hybrids. The names of male and
                                      female parents must match the rownames of inbred_gen. Missing (NA) valuesare not allowed.'))),
                fluidRow(column(width=4,radioButtons('inbredpheinput1',label = h4('A matrix of a phenotypic values of parent (Optional)'),
                                                     c('Not included'='NULL','Input parent phenotype'='input'))),
                         column(width=8,conditionalPanel(condition = 'input.inbredpheinput1=="input"',
                                                         fileInput('inbred_phe1',label=h4('Parent Phenotype')))))),
             tabPanel(h2('Select methods & models'),title = 'Methods & Models',
                fluidRow(column(width=6,selectInput('method',label=h4('method,eight GS methods'),
                                         choices = list('GBLUP'='GBLUP','BayesB'='BayesB',
                                         'RKHS'='RKHS','PLS'='PLS',
                                         'LASSO'='LASSO','EN'='EN',
                                         'XGBoost'='XGBoost','LightGBM'='LightGBM'),
                                          selected = 'GBLUP'))),
                fluidRow(column(width = 6,selectInput('model',label=h4('the prediction model'),
                                          choices = list('the additive model'='A',
                                                         'the additive-dominance model'='AD',
                                                         'the additive-phenotypic model'='A-P',
                                                         'the additive-dominance-phenotypic model'='AD-P'),
                                          selected = 'the additive model')))),
             tabPanel(h2('Select hybrids'),title = 'Selection',
                fluidRow(column(width=6,selectInput('select',label = h4('the selection of hybrids based on the prediction results'),
                                        choices = list('all potential crosses'='all',
                                        'the top n crosses'='top',
                                        'the bottom n crosses'='bottom')))),
                fluidRow(column(width=6,numericInput('number',label=h4('the number of selected top or bottom hybrids,only when select = "top" or select = "bottom".'),
                                                     value = 100)))),
             tabPanel(h2('Phenotypic values of the predicted hybrids'),title='Phenotypic values',
                fluidRow(column(width = 6,downloadLink('predres',label = h4('Predict & Download Results')))),
                fluidRow(column(width=12,DT::dataTableOutput('predhyres1')))))),
          tabPanel(h4('predhy.predict_NCII'),
          navlistPanel(widths = c(3,9),
             tabPanel(h2('Predict the Performance of Hybrids'),title = 'Description',
                      helpText('Predict all potential crosses of a given set of parents using a subset of crosses as the training sample.')),
             tabPanel(h2('Input genotype & phenotype'),title = 'Input files',
                fluidRow(column(width=4,fileInput('inbred_gen_NCII',label = h4('inbred_gen'))),
                         column(width=8,helpText('A matrix for genotypes of parental lines in numeric format, coded as 1, 0 and
                                      -1. The row.names of inbred_gen must be provied. It can be obtained from the original genotype using',code('convertgen'),'function'))),
                fluidRow(column(width = 4,fileInput('hybrid_phe_NCII',label = h4('hybrid_phe'))),
                                column(width=8,helpText('A data frame with three columns. The first column and the second column are the
                                      names of male and female parents of the corresponding hybrids, respectively;the third column is the phenotypic values of hybrids. The names of male and
                                      female parents must match the rownames of inbred_gen. Missing (NA) valuesare not allowed.'))),
                fluidRow(column(width=4,radioButtons('inbredpheinput2',label = h4('A matrix of a phenotypic values of parent (Optional)'),
                                                     c('Not included'='NULL','Input parent phenotype'='input'))),
                         column(width=8,conditionalPanel(condition = 'input.inbredpheinput2=="input"',
                                                         fileInput('inbred_phe2',label=h4('Parent Phenotype')))))),
             tabPanel(h2('Input names of parents'),title = 'Parent names',
                fluidPage(column(width = 6,fileInput('male_name',label = h4('male_name'))),
                          column(width = 6,helpText('a vector of the names of male parents.'))),
                fluidPage(column(width = 6,fileInput('female_name',label = h4('female_name'))),
                          column(width = 6,helpText('a vector of the names of female parents.')))),
             tabPanel(h2('Select methods & models'),title = 'Methods & Models',
                fluidRow(column(width=6,selectInput('method_NCII',label=h4('method,eight GS methods'),
                                        choices = list('GBLUP'='GBLUP','BayesB'='BayesB',
                                                       'RKHS'='RKHS','PLS'='PLS','LASSO'='LASSO','EN'='EN',
                                                       'XGBoost'='XGBoost','LightGBM'='LightGBM'),
                                                  selected = 'GBLUP'))),
                fluidRow(column(width = 6,selectInput('model_NCII',label=h4('the prediction model'),
                                choices = list('the additive model'='A',
                                               'the additive-dominance model'='AD',
                                               'the additive-phenotypic model'='A-P',
                                               'the additive-dominance-phenotypic model'='AD-P'),
                                selected = 'the additive model')))),
             tabPanel(h2('Select hybrids'),title = 'Selection',
                fluidRow(column(width=6,selectInput('select_NCII',label = h4('the selection of hybrids based on the prediction results'),
                                        choices = list('all potential crosses'='all',
                                                       'the top n crosses'='top',
                                                       'the bottom n crosses'='bottom')))),
                fluidRow(column(width=6,numericInput('number_NCII',label=h4('the number of selected top or bottom hybrids,only when select = "top" or select = "bottom".'),
                                                     value = 100)))),
             tabPanel(h2('Phenotypic values of the predicted hybrids'),title='Phenotypic values',
                fluidRow(column(width = 6,downloadLink('predres_NCII',label = h4('Predict & Download Results')))),
                fluidRow(column(width=12,DT::dataTableOutput('predhyres_NCII'))))
             )),
           tabPanel(h4('convertgen'),
           navlistPanel(widths = c(3,9),
              tabPanel(h2('Convert Genotype'),title = 'Description',
                       helpText('Convert genotypes in HapMap format or in numeric format for hybrid package.')),
              tabPanel(h2('Input genotype'),title = 'Input files',
                 fluidRow(column(width = 5,selectInput('type',label=h4('file type'),
                                 choices = list('HapMap format with single bit'='hHapMap format with single bit',
                                                'HapMap format with double bit'='HapMap format with double bit',
                                                'numeric format'='numeric format'))),
                          column(width = 7,helpText('the type of genotype. There are three options: "hmp1" for genotypes in HapMap format with single bit, 
                                                    "hmp2" for genotypes in HapMap format with double bit, and "num" for genotypes in numeric format.'))),
                 fluidRow(column(width=5,fileInput('geno',label = h4('input_geno'))),
                          column(width=7,helpText('genotype in HapMap format or in numeric format. 
                                                  The names of individuals should be provided. Missing (NA) values are allowed.')))),
              tabPanel(h2('Parameters for SNPs'),title = 'Parameters',
                 fluidRow(column(width = 4,sliderInput('missingrate',label=h4('missingrate'),
                                                       min=0,max=0.3,value=0.2)),
                          column(width = 8,helpText('max missing percentage for each SNP, default is 0.2.'))),
                 fluidRow(column(width = 4,sliderInput('maf',label=h4('maf'),
                                                       min=0,max=0.2,value=0.05)),
                          column(width = 8,helpText('minor allele frequency for each SNP, default is 0.05.'))),
                 fluidRow(column(width = 6,checkboxInput('impute',label='Impute',value = TRUE)),
                          column(width = 8,helpText('logical. If TRUE, imputation. Default is TRUE.')))),
              tabPanel(h4('Converted genotype'),title = 'Results',
                 fluidRow(column(width=12,DT::dataTableOutput('convertview'))),
                 fluidRow(column(width=6,downloadLink('convered',label=h4('Download Genotype'))))))),
          
          tabPanel(h4('crodesign'),
            navlistPanel(widths = c(3,9),
              tabPanel(h2('Generate Mating Design'),title = 'Description',
                      helpText('Generate a mating design for a subset of crosses based on a balanced random partial rectangle cross-design (BRPRCD) (Xu et al. 2016).')),
              tabPanel(h2('Input parent names '),title = 'Parent names',
                       fluidRow(column(width = 4,fileInput('male_parents_name',label = h4('male parent name'))),
                                column(width = 8,helpText('a table containing names of male parents'))),
                       fluidRow(column(width = 4,fileInput('female_parents_name',label = h4('female parent name'))),
                                column(width = 8,helpText('a table containing names of female parents')))
                       ),
              tabPanel(h2('Parameters'),title = 'Input parameters',
                      fluidRow(column(width = 3,numericInput('d',label=h4('percentage'),value = 50)),
                               column(width = 9,helpText('an integer denoting 1/d percentage of crosses to be evaluated in the field.'))
                               ),
                      fluidRow(column(width = 3,offset=0,numericInput('seed_cd',label =h4('seed'),value = 123)),
                               column(width = 9,helpText('the random number')))
                      ),
              
              tabPanel(h2('Results'),title = 'Results',
                       fluidRow(column(width =4,downloadLink('crodesign_download',label = h4('Download crodesign')))),
                       fluidRow(column(width = 12,DT::dataTableOutput('crodesign'))))))
    ))
    
###################################server#######################################    
  server <- function(input, output) {
    options(shiny.maxRequestSize=10*1024^3) #max upload = 10G
####################################cv##########################################    
  fix <- reactive({
    req(input$fix)
    fix <-as.matrix(fread(input$fix$datapath,header = T,stringsAsFactors=F))
  })
    
  inbred <- reactive({
    req(input$inbred_gen1)
    inbred <-as.data.frame(fread(input$inbred_gen1$datapath,header = T,stringsAsFactors=F))
    row.names(inbred) <- inbred[,1]
    inbred <- inbred[,-1]
  })
  
  hybrid <- reactive({
    req(input$hybrid_phe1)
    hybrid <-as.data.frame(fread(input$hybrid_phe1$datapath,header = T,stringsAsFactors=F))
  })
  
  inbred_phe <- reactive({
    req(input$inbred_phe)
    inbred_phe <-as.data.frame(fread(input$inbred_phe$datapath,header = T,stringsAsFactors=F))
    row.names(inbred_phe) <- inbred_phe[,1]
    inbred_phe <- as.matrix(inbred_phe[,2,drop=FALSE])
  })
  
  infered <- reactive({
    infergen(inbred(),hybrid())
  })
  
  gena <- reactive({
    infered()$add
  })

  gend <- reactive({
    infered()$dom
  })

  cvmethod <- reactive({switch(input$cvmethod,
                               'GBLUP'='GBLUP','BayesB'='BayesB',
                               'RKHS'='RKHS','PLS'='PLS',
                               'LASSO'='LASSO','EN'='EN',
                               'XGBoost'='XGBoost','LightGBM'='LightGBM',
                               'ALL'='ALL')})

  nfold <- reactive({
    input$nfold
  })

  ntimes <- reactive({
    input$ntimes
  })

  cvseed <- reactive({
    input$cvseed
  })

  cpu <- reactive({
    input$cpu
  })

  cvres <- reactive({
    if(input$cvfixinput=='input'){
      if(input$cvgendinput=='input'){
        if(input$cvinbredpheinput=='input'){
          cv(fix=fix(),gena = gena(),gend = gend(),parent_phe = inbred_phe(),hybrid_phe =hybrid(),method=cvmethod(),nfold=nfold(),nTimes=ntimes(),seed=cvseed(),CPU=cpu(),drawplot = F)
        }else{
          cv(fix=fix(),gena = gena(),gend = gend(),parent_phe = NULL,hybrid_phe=hybrid(),method=cvmethod(),nfold=nfold(),nTimes=ntimes(),seed=cvseed(),CPU=cpu(),drawplot = F)
        }
      }else{
        if(input$cvinbredpheinput=='input'){
          cv(fix=fix(),gena = gena(),gend = NULL,parent_phe = inbred_phe(),hybrid_phe=hybrid(),method=cvmethod(),nfold=nfold(),nTimes=ntimes(),seed=cvseed(),CPU=cpu(),drawplot = F)
        }else{
          cv(fix=fix(),gena = gena(),gend = NULL,parent_phe = NULL,hybrid_phe=hybrid(),method=cvmethod(),nfold=nfold(),nTimes=ntimes(),seed=cvseed(),CPU=cpu(),drawplot = F)
        }

      }
    }else{
      if(input$cvgendinput=='input'){
        if(input$cvinbredpheinput=='input'){
          cv(fix=NULL,gena = gena(),gend = gend(),parent_phe = inbred_phe(),hybrid_phe=hybrid(),method=cvmethod(),nfold=nfold(),nTimes=ntimes(),seed=cvseed(),CPU=cpu(),drawplot = F)
        }else{
          cv(fix=NULL,gena = gena(),gend = gend(),parent_phe = NULL,hybrid_phe=hybrid(),method=cvmethod(),nfold=nfold(),nTimes=ntimes(),seed=cvseed(),CPU=cpu(),drawplot = F)
        }
      }else{
        if(input$cvinbredpheinput=='input'){
          cv(fix=NULL,gena = gena(),gend = NULL,parent_phe = inbred_phe(),hybrid_phe=hybrid(),method=cvmethod(),nfold=nfold(),nTimes=ntimes(),seed=cvseed(),CPU=cpu(),drawplot = F)
        }else{
          cv(fix=NULL,gena = gena(),gend = NULL,parent_phe = NULL,hybrid_phe=hybrid(),method=cvmethod(),nfold=nfold(),nTimes=ntimes(),seed=cvseed(),CPU=cpu(),drawplot = F)
        }
      }
    }
  })

  output$cv1 <- renderText({
    if(cvmethod()!='ALL'){
      paste(cvmethod(),'R^2=',mean(cvres()))}
  })

  output$cvp <-renderPlot({
    if(cvmethod()=='ALL'){
      mycolor <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
      plotres <- function(d){
        names <- attr(d,'dimnames')[[2]]
        attr(d,'dimnames') <- NULL
        d <- as.matrix(d)
        d <- apply(d,2,mean)
        barplot(height = d,names.arg = names,col=mycolor,las=2,ylim = c(0,max(d)+0.05),xlim=c(0,10),
                main = 'Trait predictability of 8 methods',ylab = "R2")
        locat <- seq(0.75,9,length.out=8)
        text(locat,d+0.015,round(d,3),cex=0.9)
      }
      plotres(cvres())
    }
  })
  #####################################predhy.predict#############################  
  inbred_gen <-reactive({
    req(input$inbred_gen)
    fix <-as.data.frame(fread(input$inbred_gen$datapath,header = T,stringsAsFactors=F))
    row.names(fix) <- fix[,1]
    fix <- as.matrix(fix[,-1])
  })
  
  hybrid_phe <-reactive({
    req(input$hybrid_phe)
    fix <-as.data.frame(fread(input$hybrid_phe$datapath,header = T,stringsAsFactors=F))
  })
  
  inbred_phe1 <- reactive({
    req(input$inbred_phe1)
    inbred_phe1 <-as.data.frame(fread(input$inbred_phe1$datapath,header = T,stringsAsFactors=F))
    row.names(inbred_phe1) <- inbred_phe1[,1]
    inbred_phe1 <- as.matrix(inbred_phe1[,2,drop=FALSE])
  })
  
  method <- reactive({switch(input$method,
                             'GBLUP'='GBLUP','BayesB'='BayesB',
                             'RKHS'='RKHS','PLS'='PLS',
                             'LASSO'='LASSO','EN'='EN',
                             'XGBoost'='XGBoost','LightGBM'='LightGBM')
    })
  
  model <- reactive({switch(input$model,
                            'A'='A',
                            'AD'='AD',
                            'A-P'='A-P',
                            'AD-P'='AD-P')
    })
  
  select <- reactive({switch(input$select,
                             'all'='all',
                             'top'='top',
                             'bottom'='bottom')
    })
  
  number <- reactive({ifelse(select()=='all',NULL,input$number)})
  
  pred <- reactive({
    if(input$inbredpheinput1=='input'){
      predhy.predict(inbred_gen(),hybrid_phe = hybrid_phe(),parent_phe = inbred_phe1(),method(),model(),select(),number())
    }else{
      predhy.predict(inbred_gen(),hybrid_phe = hybrid_phe(),parent_phe = NULL,method(),model(),select(),number())
    }
  })
  
  output$predhyres1 <- DT::renderDataTable({
    pred()
  })
  
  output$predres <- downloadHandler(
    filename = function(){
      paste('predhy results','-',method(),Sys.Date(),'.csv', sep = "")
    },
    content = function(file) {
      write.csv(cbind(rownames(pred()),pred()), file, row.names = FALSE)
    })
  
  #################################predhy.predict_NCII############################
  inbred_gen_NCII <-reactive({
    req(input$inbred_gen_NCII)
    fix <-as.data.frame(fread(input$inbred_gen_NCII$datapath,header = T,stringsAsFactors=F))
    row.names(fix) <- fix[,1]
    fix <- fix[,-1]
  })
  
  hybrid_phe_NCII <-reactive({
    req(input$hybrid_phe_NCII)
    fix <-as.data.frame(fread(input$hybrid_phe_NCII$datapath,header = T,stringsAsFactors=F))
    })
  
  inbred_phe2 <- reactive({
    req(input$inbred_phe2)
    inbred_phe2 <-as.data.frame(fread(input$inbred_phe2$datapath,header = T,stringsAsFactors=F))
    row.names(inbred_phe2) <- inbred_phe2[,1]
    inbred_phe2 <- as.matrix(inbred_phe2[,2,drop=FALSE])
  })
  
  male_name <- reactive({
    req(input$male_name)
    fix <-c(fread(input$male_name$datapath,header = T,stringsAsFactors=F))[[1]]
    })
  
  female_name <- reactive({
    req(input$female_name)
    fix <-c(fread(input$female_name$datapath,header = T,stringsAsFactors=F))[[1]]
    })
  
  method_NCII <- reactive({switch(input$method_NCII,
                                  'GBLUP'='GBLUP','BayesB'='BayesB',
                                  'RKHS'='RKHS','PLS'='PLS',
                                  'LASSO'='LASSO','EN'='EN',
                                  'XGBoost'='XGBoost','LightGBM'='LightGBM')
    })
  
  model_NCII <- reactive({switch(input$model_NCII,
                                 'A'='A',
                                 'AD'='AD',
                                 'A-P'='A-P',
                                 'AD-P'='AD-P')
    })
  
  select_NCII <- reactive({switch(input$select_NCII,
                                  'all'='all',
                                  'top'='top',
                                  'bottom'='bottom')
    })
  
  number_NCII <- reactive({ifelse(select_NCII() == 'all',NULL,input$number_NCII)})
  
  pred_NCII <- reactive({
    if(input$inbredpheinput2=='input'){
    predhy.predict_NCII(inbred_gen = inbred_gen_NCII(),hybrid_phe = hybrid_phe_NCII(),
                        parent_phe = inbred_phe2(),
                        male_name = male_name(),female_name = female_name(),
                        method = method_NCII(),model = model_NCII(),
                        select = select_NCII(),number = number_NCII())
    }else{
    predhy.predict_NCII(inbred_gen = inbred_gen_NCII(),hybrid_phe = hybrid_phe_NCII(),
                        parent_phe = NULL,
                        male_name = male_name(),female_name = female_name(),
                        method = method_NCII(),model = model_NCII(),
                        select = select_NCII(),number = number_NCII())
    }
  })
  
  output$predhyres_NCII <- DT::renderDataTable({
    pred_NCII()
  })
  
  output$predres_NCII <- downloadHandler(
    filename = function() {
      paste('predhy results','-',method(),Sys.Date(),".csv", sep = "")
    },
    content = function(file) {
      write.csv(cbind(rownames(pred_NCII()),pred_NCII()), file, row.names = FALSE)
    })
  #################################convertgen#####################################
  filetype <- reactive({
    switch(input$type,
           'HapMap format with single bit'='hmp1',
           'HapMap format with double bit'='hmp2',
           'numeric format'='num')
  }) 
  
  rawgene<-reactive({
    if(filetype() == 'num'){
      req(input$geno)
      geno_raw <- as.data.frame(fread(input$geno$datapath,header = T,stringsAsFactors = F))
      row.names(geno_raw) <- geno_raw[,1]
      geno_raw <- geno_raw[,-1]
    }else{
      req(input$geno)
      geno_raw <- as.data.frame(fread(input$geno$datapath,header = T,stringsAsFactors = F))
    }
  })
  
  missingrate <- reactive({
    input$missingrate
  })
  
  maf <- reactive({
    input$maf
  })
  
  impute <- reactive({
    input$impute
  })
  
  convert <- reactive({
    convertgen(rawgene(),filetype(),missingrate(),maf(),impute())
  })
  
  output$convered <- downloadHandler(
    filename = function() {
      paste('convertgen results',Sys.Date(),".csv", sep = "")
    },
    content = function(file) {
      write.csv(convert(), file, row.names = T)
    })
  
  output$convertview <- DT::renderDataTable({
    convert()
  })
##########################crodesign################################################
  seed_cd <- reactive({input$seed_cd})
  d <- reactive({input$d})
  
  male_names <- reactive({
    req(input$male_parents_name)
    fix <-c(fread(input$male_parents_name$datapath,header = T,stringsAsFactors=F))[[1]]
  })
  
  female_names <- reactive({
    req(input$female_parents_name)
    fix <-c(fread(input$female_parents_name$datapath,header = T,stringsAsFactors=F))[[1]]
  })
  
  crodesignres <- reactive({
    crodesign(d(),male_names(),female_names(),seed = seed_cd())
  })
  
  output$crodesign <- DT::renderDataTable({
    return(crodesignres())
  })
  
  output$crodesign_download <- downloadHandler(
    filename = function() {
      paste('crodesign_results',Sys.Date(),".csv", sep = "")
    },
    content = function(file) {
      write.csv(crodesignres(), file, row.names = FALSE)
    })
  
}

shinyApp(ui = ui,server = server)}}
