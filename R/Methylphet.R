Methylphet <-
function(traindata.mat =NA, traindata.methyl1=NA,traindata.methyl2=NA, 
                      goldstandard=NA, ChIPseqPeaks = NA, 
                      OtherGenomicFeatures.train=NA,
                      testdata.mat =NA, testdata.methyl1=NA,testdata.methyl2=NA,
                      OtherGenomicFeatures.test=NA)
{
  require(Biostrings)
  require(GenomicRanges)
  require(IRanges)
  #require(BSgenome.Hsapiens.UCSC.hg18)
  #require(BSgenome.Mmusculus.UCSC.mm9)
  require(randomForest)
  #require(GenomicFeatures)
  
  

    
  
  
  ##### train_esti function estimates the parameters.
  
  est.beta =
    function(k, n) {
      ix = n >= 5
      k = k[ix]
      n = n[ix]
      N=length(k)
      mu=sum(k)/sum(n)
      theta=k/n
      s2=(N*sum(n*(theta-mu)^2))/((N-1)*sum(n))
      tmp=mu*(1-mu)/N*sum(1/n)
      
      if(s2<tmp) tmp2=0.001
      else tmp2=s2-tmp
      
      M=(mu*(1-mu)-s2) / tmp2
      alpha=mu*M
      beta=(1-mu)*M
      c(alpha, beta)
    }
  
  
  ######
  
  
  train_esti =
    function(data_5mc, binding_info,win.no, win.len)
    {all_data = data.frame(data_5mc,binding = binding_info)
     dat = t(rbind(as.matrix(all_data[all_data$binding==1,grep("^methy",names(all_data))]),
                   as.matrix(all_data[all_data$binding==1,grep("^total",names(all_data))])))
     ncol = ncol(dat)/2
     est.peak = as.data.frame(apply(dat, 1, function(y) est.beta(y[1:ncol], y[ncol+(1:ncol)])))
     names(est.peak) = paste("windows",1:(win.no*2+1),sep="")
     
     row.names(est.peak) = c("alpha","beta")
     
     est.nopeak = est.beta(as.vector(all_data[all_data$binding==0, grep("^methy",names(all_data))]),
                           as.vector(all_data[all_data$binding==0, grep("^total",names(all_data))]))
     alpha_peak = est.peak[1,]
     beta_peak = est.peak[2,]
     alpha_nopeak = est.nopeak[1]
     beta_nopeak = est.nopeak[2]
     paras = list(alpha_peak = alpha_peak , beta_peak = beta_peak,
                  alpha_nopeak = alpha_nopeak, beta_nopeak = beta_nopeak)
     
    }
  
  #####
  
  
  Train_machine =
    function(data_5mc,binding_flag = NA , win.no = 10,win.len = 30,peak.GR.data = FALSE,peak.GR=NA)
    {
      require(Biostrings)
      require(GenomicRanges)
      
      if (peak.GR.data == TRUE)
      {
        #peak1 = read.table(chip_datapath)
        #peak.GR = GRanges(seqnames = Rle(peak1[,1]),ranges=IRanges(start=peak1[,2]+1, end=peak1[,3]))
        sites = GRanges(seqnames=Rle(data_5mc$seqnames),
                        ranges = IRanges(start = data_5mc$start -300, end = data_5mc$end + 300),
                        motif = data_5mc$motif)
        binding_flag = as.numeric(sites %over% peak.GR)
      }
      paras = train_esti(data_5mc, binding_flag, win.no, win.len)
    }
  
  #####
  
  ##### predict function will predict which part is the binding sites.
  
  llratio =
    function(k,n,paras) {
      llr = rep(0, length(n))
      ix = (n > 5)
      if (length(which(ix == TRUE)) > 0)
      {
        k = as.matrix(k[ix])
        n = as.matrix(n[ix])
        alpha1 = as.matrix(paras$alpha_peak)[ix]
        beta1 = as.matrix(paras$beta_peak)[ix]
        ll1 = log(beta(k + alpha1, n - k + beta1)/beta(alpha1, beta1))
        ll2 = log(beta(k + paras$alpha_nopeak, n - k + paras$beta_nopeak)/beta(paras$alpha_nopeak, paras$beta_nopeak))
        llr[ix] = ll1 -ll2
      }
      llr
    }
  
  
  #####
  
  Predict =
    function(paras, data_5mc)
    {
      require(Biostrings)
      require(GenomicRanges)
      dat = as.matrix(data_5mc[,c(grep("^methy",names(data_5mc)),grep("^total",names(data_5mc)))])
      ncol = ncol(dat)/2
      llr.21 = apply(dat, 1, function(y) llratio(y[1:ncol], y[ncol+(1:ncol)],paras))     ##### 21 windows
      sites = GRanges(seqnames=Rle(data_5mc$seqnames),
                      ranges = IRanges(start = data_5mc$start, end = data_5mc$end),
                      motif = data_5mc$motif)
      elementMetadata(sites)[["predict"]] = colSums(llr.21)
      sites
    }
  
  
  #####
  
  Test_machine =
    function(data_5mc,paras,binding_flag = NA,peak.GR.data = FALSE,peak.GR = NA)
    {
      require(Biostrings)
      require(GenomicRanges)
      
      pred_results = Predict(paras, data_5mc)
      
      if (peak.GR.data ==TRUE)
      {
        sites = GRanges(seqnames=Rle(data_5mc$seqnames),
                        ranges = IRanges(start = data_5mc$start -300, end = data_5mc$end +300 ),
                        motif = data_5mc$motif)
        if (class(peak.GR) =="logical")
          binding_flag = NA   else
            binding_flag = as.numeric(sites %over% peak.GR)
      }
      
      Pred_all = cbind(as.data.frame(pred_results),binding_flag)
      
    }
  
  ChIPseqpeaks2goldstandard = 
    function(traindata.mat, ChIPseqPeaks )
    { require(Biostrings)
      require(GenomicRanges)
      sites.tmp = GRanges(seqnames=Rle(traindata.mat$seqnames),
                          ranges = IRanges(start = traindata.mat$start-300, end = traindata.mat$end+300),
                          motif = traindata.mat$motif)
      binding_flag = as.numeric(sites.tmp %over% ChIPseqPeaks)
      
    }
  
##################################################################################################################  
  
  if   (any(is.na(traindata.mat)) | any(is.na(traindata.methyl1)  ))
    stop("Complete training data is needed")
  
  
  if   (any(is.na(goldstandard)) & any(is.na(ChIPseqPeaks)))
    stop("Golden Standard or Chip-seq Peaks for training data is needed")
  
  
  if   (any(is.na(goldstandard)) ==1)
    goldstandard =   ChIPseqpeaks2goldstandard(traindata.mat, ChIPseqPeaks)
  
  if (any(is.na(traindata.methyl2))==1)
    mc.type = 1 else 
      mc.type = 2
  
  
  if (mc.type==1)
  {
    traindata.all = cbind(traindata.mat,traindata.methyl1)
    testdata.all = cbind(testdata.mat,testdata.methyl1)
    
    paras =
      Train_machine(traindata.all, binding_flag = goldstandard)
    
    
    Train.score =
      Test_machine(traindata.all,paras, binding_flag = goldstandard)
    
    Pred.score =
      Test_machine(testdata.all,paras)
    
    

    
    X.train.rf=cbind(Train.score$motif, Train.score$predict, OtherGenomicFeatures.train)
    colnames(X.train.rf)[c(1,2)]=c('motif.score','methylation.score')
    idx=which(is.na(X.train.rf))
    X.train.rf[idx]=0
    rf.train=randomForest(X.train.rf, factor(Train.score$binding_flag), na.action=na.omit, ntree=500)
    
    
    X.test.rf=cbind(Pred.score$motif, Pred.score$predict, OtherGenomicFeatures.test)
    colnames(X.test.rf)[c(1,2)]=c('motif.score','methylation.score' )
    idx=which(is.na(X.test.rf))
    X.test.rf[idx]=0
    rf.predict.result=predict(rf.train, X.test.rf, na.action=na.omit, type='prob')
    
    res=list(rfmodel=rf.train,
             predict = cbind(testdata.all[1:3],Methyl.score1 = Pred.score$predict,rf.predict.result),
             parameters = paras)
  }
  

  if (mc.type==2)
  {
    traindata.all1 = cbind(traindata.mat,traindata.methyl1)
    testdata.all1 = cbind(testdata.mat,testdata.methyl1)
    
    traindata.all2 = cbind(traindata.mat,traindata.methyl2)
    testdata.all2 = cbind(testdata.mat,testdata.methyl2)
    
    paras1 =
      Train_machine(traindata.all1, binding_flag = goldstandard)
    
    
    Train.score =
      Test_machine(traindata.all1,paras1, binding_flag = goldstandard)
    
    Pred.score =
      Test_machine(testdata.all1,paras1)
    
    
    paras2 =
      Train_machine(traindata.all2, binding_flag = goldstandard)
    
    
    Train.score2 =
      Test_machine(traindata.all2,paras2, binding_flag = goldstandard)
    
    Pred.score2 =
      Test_machine(testdata.all2,paras2)
    
    
    
    X.train.rf=cbind(Train.score$motif, Train.score$predict,Train.score2$predict, OtherGenomicFeatures.train)
    colnames(X.train.rf)[c(1,2,3)]=c('motif.score','methylation.score1','methylation.score2')
    idx=which(is.na(X.train.rf))
    X.train.rf[idx]=0
    rf.train=randomForest(X.train.rf, factor(Train.score$binding_flag), na.action=na.omit, ntree=500)
    
    
    X.test.rf=cbind(Pred.score$motif, Pred.score$predict,Pred.score2$predict, OtherGenomicFeatures.test)
    colnames(X.test.rf)[c(1,2,3)]=c('motif.score','methylation.score1','methylation.score2' )
    idx=which(is.na(X.test.rf))
    X.test.rf[idx]=0
    rf.predict.result=predict(rf.train, X.test.rf, na.action=na.omit, type='prob')
    
    res=list(rfmodel=rf.train,
             predict = cbind(testdata.all1[1:3],
                             Methyl.score1 = Pred.score$predict,Methyl.score2 = Pred.score2$predict,rf.predict.result),
             parameters1 = paras1,parameters2 = paras2)
  }
  
  
  res
  
}
