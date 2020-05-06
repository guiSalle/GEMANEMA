require(ggplot2)
require(lme4)
require(lmerTest)
require(lattice)
require(nlme)
require(Hmisc)
require(MASS)
require(reshape2)
require(stringr)
require(car)
require(scales)
require(ggrepel)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

theme_set(theme_bw())
colsRS = c('#e08214','#8073ac',
         '#d73027','#4575b4')
colEnv = c('#bf812d','#35978f',
           '#c51b7d','#7fbc41')

draft = '~/Documents/INRA/GEMANEMA/Paper/'
setwd('~/Documents/INRA/GEMANEMA/Paper/')

###===============================================================================================
## Create datasets for both lines
###===============================================================================================

####------=========== Create G1 dataset ======----------#####

##--- Load needed files
g1 = read.csv(file='./data/parasito_infestation1et2.csv',header=T,sep=",")

# retrieve G1 litter size
gen1 = read.csv(file='./data/G1sol0acendants.csv',header=T,sep=';')

# Weight data
pds = read.table("./data/gemanema.pesees.txt",header=T)
colnames(gen1)[1] ='IPG'

##--- Merge genetic information
g1 = merge(g1,gen1[,c('IPG','ipg_pe','mod_al','sol0')],by='IPG')
g1 = g1[which(!is.na(g1$n1)),] #-- suppr. individuals with missing data n = 6
g1 = g1[which(!is.na(g1$ht141)),]
#[1] 86 21
g1$lot2 = substr(g1$lot,5,length(g1$lot))
g1$lot = g1$lot2
g1$lot2 = NULL
g1$sexe = factor(g1$sexe)

####---- Redefine blocks
g1$gen = factor(substr(g1$lot,1,1))
g1$cpt = factor(substr(g1$lot,3,5))
table(g1$gen,g1$sexe,g1$cpt)
g1$n1 = NULL
g1$n2 = NULL
g1$n3 = NULL
g1$n4 = NULL

###-- Pds
pds$IPG = factor(pds$agneau+20000e6)
pds$agneau = NULL
colnames(pds) = c('pds011','pds141','pds241','pds301',
                  'pds022','pds142','pds242','pds302',
                  'IPG')
g1 = merge(g1,pds,by = 'IPG')
g1$gen = paste0(g1$gen,'-G1')

####----  G1 data - without stress conditions

wk = melt(g1[grep('-nostress',g1$lot),],c('IPG','lot','sexe','gen','cpt','mod_al'))
wk$lot = factor(wk$lot)

##-- FEC dataset G1
fec=wk[wk$variable %in% c("opg241","opg301","opg242","opg302"),]
fec$day=sapply(str_split(fec$variable,"opg"),function(x) x[2])
fec$inf = 1
fec$inf[fec$day %in% c("02","242","302")]=2
fec$inf=factor(fec$inf)
fec$day=substr(fec$day,1,2)
fec$day[fec$day=='01']='0'
fec$day[fec$day=='02']='0'
fec$day=factor(fec$day)

##-- Weight G1
wgt = wk[grep('pds',wk$variable),]
wgt$day = sapply(str_split(wgt$variable,"pds"),function(x) x[2])
wgt$inf = 1
wgt$inf[wgt$day %in% c("022","142","242","302")] = 2
wgt$inf = factor(wgt$inf)
wgt$day = substr(wgt$day,1,2)
wgt$day[wgt$day=='01'] = '0'
wgt$day[wgt$day=='02'] = '0'
wgt$day = factor(wgt$day)

####------=========== Create G2 dataset ======----------#####

### Files needed for G2
g2 = read.csv(file = './data/gemanema2017.csv',header=T,sep=',')
colnames(g2)[10] = 'IPG'
g2$idxa = NULL
g2$lignee = NULL
g2$modal = NULL

# retrieve G2 litter size
gen2 = read.csv(file = './data/G2_1sol0acendants.csv',header=T,sep=';')
colnames(gen2)[2] = 'IPG'
g2 = merge(g2,gen2[,c('IPG','mod_al','sol0')],by='IPG')

# Remove useless variables
g2$parasite = NULL
g2$ht01 = NULL
vec = grep('cocc',colnames(g2))
vec = c(vec,grep('disco',colnames(g2)))
g2 = g2[,which(!(seq(1:dim(g2)[2]) %in% vec))]

# Remove discarded individuals in the course of experiment
g2 = g2[g2$inoc!='retire ',]
g2$inoc = factor(g2$inoc)

## Add individual weights
pds17 = read.csv(file = './data/poids_2017.csv',sep = ';',header = T)
pds17$date = factor(as.Date(pds17$D_PESEE,"%d/%m/%Y"))
pds17 = pds17[pds17$ID_ANIMAL %in% g2$agneau,]

pds011 = pds17[pds17$D_PESEE=='04/04/2017',c('ID_ANIMAL','POIDS')]
pds011$pds011=pds011$POIDS/100
pds011$POIDS=NULL
pds141 = pds17[pds17$D_PESEE=='18/04/2017',c('ID_ANIMAL','POIDS')]
pds141$pds141=pds141$POIDS/100
pds141$POIDS=NULL
pds301 = pds17[pds17$D_PESEE=='04/05/2017',c('ID_ANIMAL','POIDS')]
pds301$pds301=pds301$POIDS/100
pds301$POIDS=NULL
pds022 = pds17[pds17$D_PESEE=='18/05/2017',c('ID_ANIMAL','POIDS')]
pds022$pds022=pds022$POIDS/100
pds022$POIDS=NULL
pds142 = pds17[pds17$D_PESEE=='01/06/2017',c('ID_ANIMAL','POIDS')]
pds142$pds142=pds142$POIDS/100
pds142$POIDS=NULL
pds302 = pds17[pds17$D_PESEE=='19/06/2017',c('ID_ANIMAL','POIDS')]
pds302$pds302=pds302$POIDS/100
pds302$POIDS=NULL

g2$ID_ANIMAL = g2$agneau
g2 = merge(g2,pds011,by='ID_ANIMAL')
g2 = merge(g2,pds141,by='ID_ANIMAL')
g2 = merge(g2,pds301,by='ID_ANIMAL')
g2 = merge(g2,pds022,by='ID_ANIMAL')
g2 = merge(g2,pds142,by='ID_ANIMAL')
g2 = merge(g2,pds302,by='ID_ANIMAL')

##-- Reformat variable names for further comparison with g1
g2$lot = paste0(g2$gp,'-',g2$inoc)
g2$gen = g2$gp

g2$opg242 = g2$opg252 ## 1-day offset between g1 and g2 for Day24 of 2nd infection
g2$opg302 = g2$opg322 ## 2-day offset between g1 and g2 for Day30 of 2nd infection

g2$gen = factor(paste0(g2$gen, '-G2'))
g2$sexe = factor(g2$sexe)

g2 = g2[,c('IPG','pere','lot','sexe','opg241','opg301','opg302','gen','inoc','mod_al','sol0',
           'pds011', 'pds141', 'pds301', 'pds022', 'pds142', 'pds302')]

####---- Output g1 and g2 group information for G x E variance partitioning
gpinfo = rbind(g1[,c('IPG','lot')], g2[,c('IPG','lot')])
write.table(gpinfo,file='./data/group_info.tsv', quote=F,row.names=F)

####----  G2 data - H. contortus only 
wk2 = melt(g2[grep('hc',g2$lot),],c('IPG','pere','lot','sexe','gen','inoc','mod_al'))
wk2$lot = factor(wk2$lot)

##-- FEC dataset G2
fec2 = wk2[wk2$variable %in% c("opg241","opg301","opg302"),]
fec2$day = sapply(str_split(fec2$variable,"opg"),function(x) x[2])
fec2$inf = 1
fec2$inf[fec2$day %in% c("02","242","302")]=2
fec2$inf=factor(fec2$inf)
fec2$day=substr(fec2$day,1,2)
fec2$day[fec2$day=='01']='0'
fec2$day[fec2$day=='02']='0'
fec2$day=factor(fec2$day)

##-- Weight G2
wgt2 = wk2[grep('pds',wk2$variable),]
wgt2$day = sapply(str_split(wgt2$variable,"pds"),function(x) x[2])
wgt2$inf = 1
wgt2$inf[wgt2$day %in% c("022","142","242","302")] = 2
wgt2$inf = factor(wgt2$inf)
wgt2$day = substr(wgt2$day,1,2)
wgt2$day[wgt2$day=='01'] = '0'
wgt2$day[wgt2$day=='02'] = '0'
wgt2$day = factor(wgt2$day)

###---- Inbreeding within each line
inbreeding = read.csv(file = './data/inbreeding_lines.csv',header=T,sep=';')
inbreeding = inbreeding[inbreeding$Gene %in% c('_G1_','G2_1'),] ## retain G1 and G2 data

###===============================================================================================
## Achieved divergence - based on eBV estimated from G0, G1 and G2 infected by H.contortus
###===============================================================================================
colnames(gen1)[2] = 'sexe'
colnames(gen2)[12] = 'Gene'
colnames(gen2)[6] = 'fropg1'
colnames(gen2)[7] = 'fropg2'

## Add ancestral G0 individuals
nucleus = read.csv(file='./data/G0G1G2.csv',header=T,sep=';')
colnames(nucleus)[1] = 'IPG'
nucleusG0 = nucleus[nucleus$Gene == '_G0_',c('IPG','ipg_pe','ipg_me','fropg1','fropg2','line','pheno','sexe','mod_al','Gene')]

#####----- Bring 3 generations together
flock = rbind(nucleusG0[,c('IPG','ipg_pe','ipg_me','fropg1','fropg2','line','pheno','sexe','Gene','mod_al')],
              gen1[,c('IPG','ipg_pe','ipg_me','fropg1','fropg2','line','pheno','sexe','Gene','mod_al')],
              gen2[,c('IPG','ipg_pe','ipg_me','fropg1','fropg2','line','pheno','sexe','Gene','mod_al')])

flock$Gene = factor(flock$Gene)
flock$line = factor(flock$line)
flock$sexe = factor(flock$sexe)
dim(flock)
#[1] 570  10
table(flock$Gene)
# _G0_ _G1_ G2_1 
#  243  213  114 

#####----- Add individual eBVs and standardized relative to unselected G0
## Read eBVs
#ebvf = read.csv(file = './data/eBV_G0G1G2.csv',header=T,sep=';')
ebvf = read.csv(file = './data/solopg12_G0G1G2_1.csv', header = T, sep=';')
dim(ebvf)
#[1] 977   31
head(ebvf)
# fropg1 fropg2 fropg    sol_opg1  sd_opg1    sol_opg2  sd_opg2  sol_opgm12 sd_opgm12  sol_opg1G0
# 1      .      .     .  0.01838718 1.480186  0.00305889 1.327369  0.01227801  1.331044 -0.06192151
# 2      .      .     . -0.01462008 1.377549  0.00073714 1.236413  0.00849796  1.248544 -0.08235176
# 3      .      .     .  0.32570329 1.475434  0.26905392 1.323166  0.20644586  1.328067  0.48489949
# 4      .      .     .  1.17076065 1.250863  1.03729404 1.124168  1.08992846  1.135215  1.29476833
# 5      .      .     . -1.23677582 1.365019 -1.13514347 1.225201 -1.15776952  1.228204 -1.29118393
# 6      .      .     .  0.73721446 1.414123  0.66563365 1.268827  0.67666279  1.277719  0.65498462
# sd_opg1G0  sol_opg2G0  sd_opg2G0 sol_opgm12G0 sd_opgm12G0 sol_opg1_G01 sd_opg1_G01 sol_opg2_G01
# 1  1.6210377 -0.04901481 1.28306902  -0.04551125  1.22718839   0.03854319  1.68529799   0.00718625
# 2 1.49542043 -0.06517939   1.183642  -0.03157545  1.16100046  -0.05348916  1.55625188  -0.01993448
# 3 1.61256549  0.38379782 1.27636321   0.19592702  1.22416117   0.48055725  1.67817537   0.31842565
# 4 1.35675084  1.02482003   1.073884   1.05140144  1.06283971   1.32279142  1.41069719   0.96857961
# 5 1.49427306 -1.02199029 1.18273385  -0.98754045  1.13793443  -1.31380324  1.55326035  -1.00842808
# 6 1.54278707  0.51843024 1.22113305   0.49479961  1.18699206   0.80428904  1.60231715   0.60237076
# sd_opg2_G01 sol_opgm12_G01 sd_opgm12_G01 sol_opg1g_G01 sd_opg1g_G01 sol_opg2g_G01 sd_opg2g_G01
# 1  1.24510717     0.00919532    1.37383885    0.00830168   1.74406793   -0.01045111   1.23586249
# 2  1.15325119     0.01017007    1.28573632   -0.02769174   1.61295237    -0.0263695   1.14572981
# 3  1.24007105     0.22765369     1.3702726    0.42945489   1.73679539    0.28289039    1.2310484
# 4  1.04945987     1.13121949    1.16790885    1.39782146   1.46979919    0.98234761    1.0464842
# 5  1.15062096    -1.13631008    1.26772025   -1.33767375   1.61322096   -0.97140187   1.14534735
# 6  1.18602483     0.70825854    1.31662327     0.5869421   1.65857355    0.43106289   1.17719058
# sol_opgm12g_G01 sd_opgm12g_G01         ipg      ipg_pe       ipg_me Gene
# 1     -0.05141174     1.24597424 20000100038 20000180782 2.000016e+10     
# 2      0.00166487      1.1846006 20000100073 20000180816 1.817490e+13     
# 3      0.18237674     1.24494036 20000100077 20000181417 1.817490e+13     
# 4      1.17579788     1.08461379 20000100110 20000180097 2.000016e+10     
# 5     -1.01927874     1.15865016 20000100145 20000180299 2.000018e+10     
# 6       0.4159122     1.20516556 20000100165 20000181345 2.000016e+10  

colnames(ebvf)[ncol(ebvf)-3] = 'IPG'
colnames(ebvf)[4] = 'ebv_fec1'
colnames(ebvf)[6] = 'ebv_fec2'
colnames(ebvf)[8] = 'ebv_fec'
ebvf$fropg1[ebvf$fropg1=='.']=NA
ebvf$fropg2[ebvf$fropg2=='.']=NA
ebvf$fropg1=as.numeric(as.character(ebvf$fropg1))
ebvf$fropg2=as.numeric(as.character(ebvf$fropg2))

## add pedigree information for G0, G1 and G2
flock = merge(flock[,-c(4,5)],ebvf[,c('IPG','ebv_fec','ebv_fec1','ebv_fec2','fropg1','fropg2')],by='IPG')
dim(flock)
#[1] 369  13

## Estimate correlation between G0 geBV and eBV
idx0 = read.table(file = './data/G0_selection.txt',header = T)
idx0$IPG = idx0$animal
idx0 = idx0[idx0$sol_50_opg1!='.' & idx0$sol_50_opg2!='.' & idx0$sol_10_opg1!='.' & idx0$sol_10_opg2!='.',]
idx0$sol_10_opg1 = as.numeric(as.character(idx0$sol_10_opg1))
idx0$sol_10_opg2 = as.numeric(as.character(idx0$sol_10_opg2))
idx0$sol_50_opg1 = as.numeric(as.character(idx0$sol_10_opg1))
idx0$sol_50_opg2 = as.numeric(as.character(idx0$sol_50_opg2))
idx0$gebv = (idx0$sol_50_opg1+idx0$sol_50_opg2)/2
idx0 = idx0[,c('IPG','gebv','sol_50_opg1','sol_50_opg2')]
colnames(idx0) = c('IPG','gebv_fec','gebv_fec1','gebv_fec2')
  
dfcor = merge(flock[flock$Gene=='_G0_',c('IPG','ebv_fec','ebv_fec1','ebv_fec2')],
              idx0)
head(dfcor)
#           IPG    ebv_fec   ebv_fec1   ebv_fec2   gebv_fec  gebv_fec1   gebv_fec2
# 1 20000132302 -0.3475557 -0.3451051 -0.3065160 -0.5370424 -0.2996770 -0.77440776
# 2 20000132303  0.5094056  0.6884498  0.5856168  0.1001331  0.2512541 -0.05098792
# 3 20000132304  0.4751515  0.5036359  0.4581971 -0.1597739 -0.2144469 -0.10510098
# 4 20000132305  0.5674178  0.6241729  0.5136685  0.5922808  0.6519814  0.53258018
# 5 20000132306  1.9122748  2.2436112  1.9801954  1.4708703  1.5328501  1.40889043
# 6 20000132307 -1.8405980 -2.2609993 -1.9757226 -1.2077925 -1.6518506 -0.76373447

rcorr(dfcor$gebv_fec1,dfcor$ebv_fec1)
#      x    y
# x 1.00  0.9
# y  0.9 1.00
# 
# n = 241  

rcorr(dfcor$gebv_fec2,dfcor$ebv_fec2)
#      x    y
# x 1.00 0.76
# y 0.76 1.00
# 
# n = 241 

rcorr(dfcor$gebv_fec,dfcor$ebv_fec)
#      x    y
# x 1.00 0.95
# y 0.95 1.00
# 
# n = 241 

#####----- Standardized selection differential relative to unselected G0
flock$ebvs = NA
flock$ebvs1 = NA
flock$ebvs2 = NA
flock$ebvs[flock$Gene=='_G0_'] = (flock$ebv_fec[flock$Gene=='_G0_'] - mean(flock$ebv_fec[flock$Gene=='_G0_']))/sd(flock$ebv_fec[flock$Gene=='_G0_'])
flock$ebvs[flock$Gene=='_G1_'] = (flock$ebv_fec[flock$Gene=='_G1_'] - mean(flock$ebv_fec[flock$Gene=='_G0_']))/sd(flock$ebv_fec[flock$Gene=='_G0_'])
flock$ebvs[flock$Gene=='G2_1'] = (flock$ebv_fec[flock$Gene=='G2_1'] - mean(flock$ebv_fec[flock$Gene=='_G0_']))/sd(flock$ebv_fec[flock$Gene=='_G0_'])

flock$ebvs1[flock$Gene=='_G0_'] = (flock$ebv_fec1[flock$Gene=='_G0_'] - mean(flock$ebv_fec1[flock$Gene=='_G0_']))/sd(flock$ebv_fec1[flock$Gene=='_G0_'])
flock$ebvs1[flock$Gene=='_G1_'] = (flock$ebv_fec1[flock$Gene=='_G1_'] - mean(flock$ebv_fec1[flock$Gene=='_G0_']))/sd(flock$ebv_fec1[flock$Gene=='_G0_'])
flock$ebvs1[flock$Gene=='G2_1'] = (flock$ebv_fec1[flock$Gene=='G2_1'] - mean(flock$ebv_fec1[flock$Gene=='_G0_']))/sd(flock$ebv_fec1[flock$Gene=='_G0_'])

flock$ebvs2[flock$Gene=='_G0_'] = (flock$ebv_fec2[flock$Gene=='_G0_'] - mean(flock$ebv_fec2[flock$Gene=='_G0_']))/sd(flock$ebv_fec2[flock$Gene=='_G0_'])
flock$ebvs2[flock$Gene=='_G1_'] = (flock$ebv_fec2[flock$Gene=='_G1_'] - mean(flock$ebv_fec2[flock$Gene=='_G0_']))/sd(flock$ebv_fec2[flock$Gene=='_G0_'])
flock$ebvs2[flock$Gene=='G2_1'] = (flock$ebv_fec2[flock$Gene=='G2_1'] - mean(flock$ebv_fec2[flock$Gene=='_G0_']))/sd(flock$ebv_fec2[flock$Gene=='_G0_'])

aggregate(ebvs ~ Gene, FUN=summary, data=flock)
#   Gene     ebvs.Min.  ebvs.1st Qu.   ebvs.Median     ebvs.Mean  ebvs.3rd Qu.     ebvs.Max.
# 1 _G0_ -2.922200e+00 -5.543863e-01  1.030462e-02  8.279242e-18  7.264267e-01  2.559025e+00
# 2 _G1_ -3.499375e+00 -1.777529e+00 -8.815682e-01 -2.790938e-01  1.370243e+00  2.926240e+00
# 3 G2_1 -3.147376e+00 -1.884497e+00  1.185234e+00  6.344309e-03  1.959353e+00  2.998757e+00

aggregate(ebvs1 ~ Gene, FUN=summary, data=flock)
#   Gene    ebvs1.Min. ebvs1.1st Qu.  ebvs1.Median    ebvs1.Mean ebvs1.3rd Qu.    ebvs1.Max.
# 1 _G0_ -2.922734e+00 -5.476042e-01 -1.600691e-03 -5.392010e-19  7.114556e-01  2.255558e+00
# 2 _G1_ -3.543297e+00 -1.775928e+00 -7.447825e-01 -3.222340e-01  1.288236e+00  2.648643e+00
# 3 G2_1 -3.036941e+00 -1.715921e+00  8.690947e-01 -5.366113e-02  1.754750e+00  2.850034e+00

aggregate(ebvs2 ~ Gene, FUN=summary, data=flock)
# Gene    ebvs2.Min. ebvs2.1st Qu.  ebvs2.Median    ebvs2.Mean ebvs2.3rd Qu.    ebvs2.Max.
# 1 _G0_ -2.837522e+00 -5.442363e-01  5.634774e-03 -1.000681e-17  7.290716e-01  2.270988e+00
# 2 _G1_ -3.525071e+00 -1.754825e+00 -7.239059e-01 -3.110439e-01  1.276259e+00  2.657170e+00
# 3 G2_1 -3.037247e+00 -1.730099e+00  9.288666e-01 -4.399064e-02  1.728604e+00  2.838061e+00

###---- Phenotype across generations

## For G2 retain the sole H.contortus infected individuals and associated FEC
vecids = g2$IPG[g2$inoc=='tc']
phenodiv = na.omit(flock[!(flock$IPG %in% vecids),])
phenodiv$fropg1[match(g2$IPG[g2$inoc=='hc'],phenodiv$IPG)] = sqrt(sqrt((g2$opg241[g2$inoc=='hc'] + g2$opg301[g2$inoc=='hc'])/2))
phenodiv$fropg2[match(g2$IPG[g2$inoc=='hc'],phenodiv$IPG)] = sqrt(sqrt(g2$opg302[g2$inoc=='hc']))
phenodiv$line = factor(phenodiv$line)
dim(phenodiv)
#[1] 367   16
phenodiv = na.omit(phenodiv)
dim(phenodiv)
#[1] 367   16

## Average FEC across two infections
phenodiv$mfec = sqrt(sqrt((phenodiv$fropg1^4 + phenodiv$fropg2^4)/2))
phenodiv$Gene = factor(substr(gsub('_','',phenodiv$Gene),1,2))
table(phenodiv$Gene)
#  G0  G1  G2 
# 241  87  39

## Check that 4th root performs better than log
shapiro.test(phenodiv$fropg1)
# data:  phenodiv$fropg1
# W = 0.96911, p-value = 5.024e-07
shapiro.test(phenodiv$fropg2)
# data:  phenodiv$fropg2
# W = 0.90632, p-value = 2.746e-14
shapiro.test(phenodiv$mfec)
# data:  phenodiv$mfec
# W = 0.94868, p-value = 5.518e-10
shapiro.test(log(50+phenodiv$fropg1^4))
# data:  log(50 + phenodiv$fropg1^4)
# W = 0.90696, p-value = 3.107e-14
shapiro.test(log(50+phenodiv$fropg2^4))
# data:  log(50 + phenodiv$fropg2^4)
# W = 0.84082, p-value < 2.2e-16
shapiro.test(log(50+phenodiv$mfec^4))
# data:  log(50 + phenodiv$mfec^4)
# W = 0.87106, p-value < 2.2e-16

######------ Correct for environmental factors in selected animals & scaled for G0 dist.
## This aims to get rid of differences in value scales between generations
## Has to be done across G0-G2 including unselected G0 individuals (base population)

## Mean FEC in selected animals
mpdiv = lme(mfec ~ mod_al + sexe + Gene,
            random =~ 1|IPG,
            data = phenodiv)
summary(mpdiv)
## Include individual random effects, standardized by G0 mean and std
phenodiv$fm = mpdiv$residuals[,2]
phenodiv$fms = (phenodiv$fm - mean(phenodiv$fm[phenodiv$Gene=='G0']))/sd(phenodiv$fm[phenodiv$Gene=='G0'])

## Mean FEC at 1st infection
mpdiv1 = lme(fropg1 ~ mod_al + sexe +Gene,
             random =~ 1|IPG,
             data = phenodiv)

## Add estimates to phenodiv and scaled by G0 values
phenodiv$f1 = mpdiv1$residuals[,2]
phenodiv$fs1 = (phenodiv$f1 - mean(phenodiv$f1[phenodiv$Gene=='G0']))/sd(phenodiv$f1[phenodiv$Gene=='G0'])

## Mean FEC at 2nd infection and scaled by G0 values
mpdiv2 = lme(fropg2 ~ mod_al + sexe + Gene,
             random =~ 1|IPG,
             data = phenodiv)

## Include individual random effects, standardized by G0 mean and std
phenodiv$f2 = mpdiv2$residuals[,2]
phenodiv$fs2 = (phenodiv$f2 - mean(phenodiv$f2[phenodiv$Gene=='G0']))/sd(phenodiv$f2[phenodiv$Gene=='G0'])

colnames(phenodiv)[4]='line'

#####------ PHENOTYPIC divergence standardized to G0 units
pdiv = aggregate(fms ~ line + Gene,data = phenodiv,FUN = mean)
pdiv$std = aggregate(fms ~ line + Gene,data = phenodiv,FUN = sd)[,3]
pdiv
#   line Gene         fms       std
# 1    0   G0  0.03143868 0.8085665
# 2    R   G0 -0.83950077 0.9866936
# 3    S   G0  0.70840300 0.7071977
# 4    R   G1 -0.89144640 0.9759958
# 5    S   G1  1.00015937 0.9477335
# 6    R   G2 -0.93844737 1.2447304
# 7    S   G2  0.89152500 0.5141770

G1div = pdiv$fms[pdiv$Gene=='G1' & pdiv$line=='S']-pdiv$fms[pdiv$Gene=='G1' & pdiv$line=='R']
G1div
#[1] 1.891606

G2div = pdiv$fms[pdiv$Gene=='G2' & pdiv$line=='S']-pdiv$fms[pdiv$Gene=='G2' & pdiv$line=='R']
G2div
#[1] 1.829972

###---- Test phenotypic divergence between lines
summary(lm(fms ~ line,data = phenodiv[phenodiv$line !='0' & phenodiv$Gene=='G1',]))
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -0.8914     0.1420  -6.280 1.38e-08 ***
# lineS         1.8916     0.2068   9.148 2.71e-14 ***

summary(lm(fms ~ line,data = phenodiv[phenodiv$line !='0' & phenodiv$Gene=='G2',]))
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -0.9384     0.2164  -4.337 0.000107 ***
# lineS         1.8300     0.3021   6.057 5.26e-07 ***
 
###---- Test within line divergence significance relative to unselected G0
summary(lm(fms ~ Gene,data = phenodiv[phenodiv$line %in% c('0','R'),]))
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.26717    0.07489  -3.567 0.000436 ***
# GeneG1      -0.62428    0.16416  -3.803 0.000182 ***
# GeneG2      -0.67128    0.23931  -2.805 0.005449 ** 

summary(lm(fms ~ Gene,data = phenodiv[phenodiv$line %in% c('0','S'),]))
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.27829    0.06217   4.476 1.18e-05 ***
# GeneG1       0.72187    0.14467   4.990 1.16e-06 ***
# GeneG2       0.61324    0.19709   3.111  0.00209 ** 

t.test(fms ~ Gene,data = phenodiv[phenodiv$line %in% c('0','R') & phenodiv$Gene %in% c('G0','G1'),])
#df = 69.88
t.test(fms ~ Gene,data = phenodiv[phenodiv$line %in% c('0','R') & phenodiv$Gene %in% c('G0','G2'),])
#df = 20.415
t.test(fms ~ Gene,data = phenodiv[phenodiv$line %in% c('0','S') & phenodiv$Gene %in% c('G0','G1'),])
#df = 55.021
t.test(fms ~ Gene,data = phenodiv[phenodiv$line %in% c('0','S') & phenodiv$Gene %in% c('G0','G2'),])
#df 31.487
  
###---- From now on keep the only selected individuals
phenodiv = phenodiv[phenodiv$line!='0',]
phenodiv$line = factor(phenodiv$line)

######------ Observed GENETIC divergence
tab = aggregate(ebvs ~ line + Gene, data = phenodiv, FUN = mean)
tab$std = aggregate(ebvs ~ line + Gene, data = phenodiv, FUN = sd)[,3]

# Gap between genetic values of both lines 
#// these values are scaled (X-meanG0)/sdG0: so given in genetic std.
tab
#   line Gene       ebvs       std
# 1    R   G0 -0.9227285 0.8479890
# 2    S   G0  0.8769739 0.6924633
# 3    R   G1 -1.7606129 0.6581164
# 4    S   G1  1.3830982 0.8345231
# 5    R   G2 -1.9353769 0.9090914
# 6    S   G2  1.8509795 0.4922208

## G1 genetic divergence
(tab$ebvs[tab$Gene=='G1' & tab$line=='S'] - tab$ebvs[tab$Gene=='G1' & tab$line=='R'])
#[1] 3.143711

## G2 genetic divergence // inflated because of the contribution of G1 individuals ?
(tab$ebvs[tab$Gene=='G2' & tab$line=='S'] - tab$ebvs[tab$Gene=='G2' & tab$line=='R'])
#[1] 3.786356

###--- Expected genetic gain
### Estimated heritability in the whole G0, G1 and G2 population (369 individuals)
## h2 = 0.44 for FEC1
## h2 = 0.33 for FEC2
## h2 = 0.49 for avg. FEC across two infections

## in G1: R = h2*S = h2*i*sp = h*i*sg with h = sqrt(h2)
## This is the maximal value it can takes; with ebv, this is constrained by accuracy not h
# i = 1/2*(if + im)
# G1_R: 2% selection on males (i = 2.421) // 47% selection on females (i = 0.846)
# G1_S: 2% selection on males (i = 2.421) // 53% selection females ( i = (1-p)/p*i ; p = 0.47/0.53*0.846)

###-- Average of respective gains across infections
ER_R1 = sqrt(.44)*0.5*(2.421+0.846)
ER_S1 = sqrt(.44)*0.5*(2.421+0.47/0.53*0.846) 
ER_R1 + ER_S1
#[1] 2.135318

### Expected response in G2 = sum from the G1 x G1 mating and G1 x G0 matings
# G2_R: 5% selection on males (i = 2.063) //35% on G0 females (i = 1.058) & 96% selection on G1 females G1 (i = (1-0.96)/0.96*2.154)
# G2_S: 5% selection on males (i = 2.063) //30% on G0 females (i = 1.259) & 100% selection G1 females G1 (i = 0)
ER_R1.1 = sqrt(0.44)*0.5*(2.063 + (1-0.96)/0.96*2.154) #+ ## G1 matings FEC1
#             sqrt(0.33)*0.5*(2.063 + (1-0.96)/0.96*2.154))/2 ## G1 matings FEC2
ER_R1.0 = sqrt(0.44)*0.5*(2.063) + sqrt(0.44)*0.5*(1.058) #+ ## G0 matings FEC1
#             sqrt(0.33)*0.5*(2.063) + sqrt(0.33)*0.5*(1.058))/2 ## G0 matings FEC2
ER_S1.1 = sqrt(0.44)*0.5*(2.063+0) #+ sqrt(0.33)*0.5*(2.063+0))/2
ER_S1.0 = sqrt(0.44)*0.5*(2.063) + sqrt(0.44)*0.5*(1.259) #+
          #   sqrt(0.33)*0.5*(2.063) + sqrt(0.33)*0.5*(1.259))/2

mean(c(ER_R1.1, ER_R1.0)) + mean(c(ER_S1.1 + ER_S1.0))
#[1] 2.660555

#####----- Observed phenotypic divergence
tab2 = aggregate(fms ~ line + Gene, data = phenodiv, FUN = mean)
tab2$std = aggregate(fms ~ line + Gene, data = phenodiv, FUN = sd)[,3]

# Gap between genetic values of both lines 
#// these values are scaled (X-meanG0)/sdG0: so given in genetic std.
tab2
#   line Gene        fms       std
# 1    R   G0 -0.8395008 0.9866936
# 2    S   G0  0.7084030 0.7071977
# 3    R   G1 -0.8914464 0.9759958
# 4    S   G1  1.0001594 0.9477335
# 5    R   G2 -0.9384474 1.2447304
# 6    S   G2  0.8915250 0.5141770

## G1 phenotypic divergence
(tab2$fms[tab2$Gene=='G1' & tab2$line=='S'] - tab2$fms[tab2$Gene=='G1' & tab2$line=='R'])
#[1] 1.891606

summary(lm(fms ~ line,data = phenodiv[phenodiv$Gene=='G1',]))
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -0.8914     0.1420  -6.280 1.38e-08 ***
# lineS         1.8916     0.2068   9.148 2.71e-14 ***
 
## G2 genetic divergence 
(tab2$fms[tab2$Gene=='G2' & tab2$line=='S'] - tab2$fms[tab2$Gene=='G2' & tab2$line=='R'])
#[1] 1.829972
summary(lm(fms ~ line,data = phenodiv[phenodiv$Gene=='G2',]))
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -0.9384     0.2164  -4.337 0.000107 ***
# lineS         1.8300     0.3021   6.057 5.26e-07 ***
  
colnames(phenodiv)[4] = 'Line'
#####--- Plot phenotypic divergence through time (corrected phenotypes)
require(dplyr)
df = data.frame(phenodiv %>% 
                  group_by(Gene, Line) %>% 
                  summarize(avg_fecs = mean(fms), n=n(), sd =sd(fms), se=sd/sqrt(n)))
df
#   Gene Line   avg_fecs  n        sd         se
# 1   G0    R -0.8416823 60 0.9860541 0.12729904
# 2   G0    S  0.7061850 66 0.7074358 0.08707933
# 3   G1    R -0.8918226 46 0.9760485 0.14391051
# 4   G1    S  1.0005814 41 0.9471925 0.14792661
# 5   G2    R -0.9635068 19 1.2647033 0.29014285
# 6   G2    S  0.9153315 20 0.5146123 0.11507081

pp = ggplot(df,aes(x = Gene, y = avg_fecs, col = Line)) +
  geom_hline(yintercept = 0, col='grey',lty = 2) +
  geom_point(alpha = .6,size = 4) +
  geom_errorbar(aes(ymin = df$avg_fecs-df$se, ymax = df$avg_fecs + df$se,width = 0.02)) +
  scale_y_continuous(limits = c(-1.5,1.5), breaks = seq(-1.5,1.5,.5)) +
  scale_color_manual(values = c(colsRS[1],colsRS[2])) +
  theme(text = element_text(size = 14)) +
  xlab('Generation') +
  ylab(bquote('Standardized average FEC ('*sigma[p]*')'))

##--- Plot genetic divergence through time
df2 = data.frame(phenodiv[phenodiv$Line!='0',] %>% 
                  group_by(Gene, Line) %>% 
                  summarize(avg_ebvs = mean(ebvs), n=n(), sd =sd(ebvs), se=sd/sqrt(n)))
df2
#   Gene Line   avg_ebvs  n        sd         se
# 1   G0    R -0.9227285 60 0.8479890 0.10947491
# 2   G0    S  0.8769739 66 0.6924633 0.08523635
# 3   G1    R -1.7606129 46 0.6581164 0.09703396
# 4   G1    S  1.3830982 41 0.8345231 0.13033061
# 5   G2    R -1.9353769 19 0.9090914 0.20855987
# 6   G2    S  1.8509795 20 0.4922208 0.11006392

pg = ggplot(df2,aes(x = Gene, y = avg_ebvs, col = Line)) +
  geom_hline(yintercept = 0, col='grey',lty = 2) +
  geom_point(alpha = .6,size = 4) +
  geom_errorbar(aes(ymin = df2$avg_ebvs-df2$se,ymax = df2$avg_ebvs + df2$se,width = 0.02)) +
  scale_y_continuous(limits = c(-2.2,2), breaks = seq(-2.2,2,.5)) +
  scale_color_manual(values = c(colsRS[1],colsRS[2])) +
  theme(text = element_text(size = 14)) +
  xlab('Generation') +
  ylab(bquote('Standardized Breeding Value ('*sigma[g]*')'))

###-- Figure 1 - Achieved divergence
pdf(file = paste0(draft,'/Figure1.pdf'))
multiplot(pp,pg,cols=1)
dev.off()
colnames(phenodiv)[4]='line'

## Summary statistics for raw faecal egg counts
aggregate((fropg1^4+fropg2^4)/2 ~ line + Gene,FUN = mean,data = phenodiv)
#   line Gene (fropg1^4 + fropg2^4)/2
# 1    R   G0               4578.3958
# 2    S   G0              13381.4394
# 3    R   G1                424.0122
# 4    S   G1               3177.2834
# 5    R   G2               1000.6579
# 6    S   G2               3999.3750

## FEC variance across generations
aggregate(fms ~ line + Gene,FUN = var,data = phenodiv)
# line Gene       fms
# 1    R   G0 0.9735642
# 2    S   G0 0.5001286
# 3    R   G1 0.9525679
# 4    S   G1 0.8981988
# 5    R   G2 1.5493538
# 6    S   G2 0.2643779

aggregate(fropg1 ~ line + Gene,FUN = sd,data = phenodiv)
# line Gene   fropg1
# 1    0   G0 1.754939
# 2    R   G0 2.587705
# 3    S   G0 1.927481
# 4    R   G1 2.277235
# 5    S   G1 1.967647
# 6    R   G2 2.341706
# 7    S   G2 1.527282

aggregate(fropg2 ~ line + Gene,FUN = sd,data = phenodiv)
# line Gene   fropg2
# 1    0   G0 1.934833
# 2    R   G0 2.497168
# 3    S   G0 1.573784
# 4    R   G1 2.017517
# 5    S   G1 2.730001
# 6    R   G2 3.497754
# 7    S   G2 1.074861

aggregate((fropg1+fropg2)/2 ~ Gene,FUN = sd,data = phenodiv)
#   Gene (fropg1 + fropg2)/2
# 1   G0            2.258878
# 2   G1            2.458581
# 3   G2            2.640388

table(phenodiv$line,phenodiv$Gene,phenodiv$sexe)
# , ,  = 1
# 
# 
#    G0  G1  G2
# 0 115   0   0
# R   3  21  10
# S   3  21  10
# 
# , ,  = 2
# 
# 
#    G0  G1  G2
# 0   0   0   0
# R  57  25   9
# S  63  20  10

#####---- Plot midparents vs. offspring G1
#Initialize midparent value
phenodiv$midparent_fopg1 = -99
phenodiv$midparent_fopg2 = -99
phenodiv$midparent_fopg = -99

## G1 Mid parent
phenodiv$midparent_fopg1[phenodiv$Gene=='G1'] = (phenodiv$fs1[match(phenodiv$ipg_pe[phenodiv$Gene=='G1'],phenodiv$IPG)]+phenodiv$fs1[match(phenodiv$ipg_me[phenodiv$Gene=='G1'],phenodiv$IPG)])/2
phenodiv$midparent_fopg2[phenodiv$Gene=='G1'] = (phenodiv$fs2[match(phenodiv$ipg_pe[phenodiv$Gene=='G1'],phenodiv$IPG)]+phenodiv$fs2[match(phenodiv$ipg_me[phenodiv$Gene=='G1'],phenodiv$IPG)])/2
phenodiv$midparent_fopg[phenodiv$Gene=='G1'] = (phenodiv$fms[match(phenodiv$ipg_pe[phenodiv$Gene=='G1'],phenodiv$IPG)]+phenodiv$fms[match(phenodiv$ipg_me[phenodiv$Gene=='G1'],phenodiv$IPG)])/2

## G2 Mid parent
phenodiv$midparent_fopg1[phenodiv$Gene=='G2'] = (phenodiv$fs1[match(phenodiv$ipg_pe[phenodiv$Gene=='G2'],phenodiv$IPG)]+phenodiv$fs1[match(phenodiv$ipg_me[phenodiv$Gene=='G2'],phenodiv$IPG)])/2
phenodiv$midparent_fopg2[phenodiv$Gene=='G2'] = (phenodiv$fs2[match(phenodiv$ipg_pe[phenodiv$Gene=='G2'],phenodiv$IPG)]+phenodiv$fs2[match(phenodiv$ipg_me[phenodiv$Gene=='G2'],phenodiv$IPG)])/2
phenodiv$midparent_fopg[phenodiv$Gene=='G2'] = (phenodiv$fms[match(phenodiv$ipg_pe[phenodiv$Gene=='G2'],phenodiv$IPG)]+phenodiv$fms[match(phenodiv$ipg_me[phenodiv$Gene=='G2'],phenodiv$IPG)])/2

G1 = phenodiv[phenodiv$Gene=='G1',]
G2 = phenodiv[phenodiv$Gene=='G2',]

psel1 = ggplot(G1,aes(x=midparent_fopg,y=fms)) + 
  scale_color_manual(values = colsRS) +
  geom_point(size =3,alpha = .4,col = colsRS[match(G1$line,levels(G1$line))]) + 
  geom_smooth(method = 'lm',col='black',lwd = .5,alpha =.2) +
  xlab(bquote('Midparent FEC across infection ('*sigma[p]*')')) +
  ylab(bquote('G1 offspring FEC across infection ('*sigma[p]*')')) +
  theme(text = element_text(size = 14)) +
  geom_abline(slope=1, intercept=0,alpha = 0.4, lty = 3) +
  scale_x_continuous(limits = c(-3.5,3),breaks=seq(-3.5,3,1)) +
  scale_y_continuous(limits = c(-3.5,3),breaks=seq(-3.5,3,1)) 

psel2 = ggplot(G2,aes(x=midparent_fopg,y=fms)) + 
  scale_color_manual(values = colsRS) +
  geom_point(size =3,alpha = .4,col = colsRS[2+match(G2$line,levels(G2$line))]) + 
  geom_smooth(method = 'lm',col='black',lwd = .5,alpha =.2) +
  xlab(bquote('Midparent FEC across infection ('*sigma[p]*')')) +
  ylab(bquote('G2 offspring FEC across infection ('*sigma[p]*')')) +
  theme(text = element_text(size = 14)) +
  geom_abline(slope=1, intercept=0,alpha = 0.4, lty = 3) +
  scale_x_continuous(limits = c(-3.5,3),breaks=seq(-3.5,3,1)) +
  scale_y_continuous(limits = c(-3.5,3),breaks=seq(-3.5,3,1)) 

## Same plots by infection rank
psel10 = ggplot(G1,aes(x=midparent_fopg1,y=fs1)) +
  scale_color_manual(values = colsRS) +
  geom_point(size =3,alpha = .4,col = colsRS[match(G1$line,levels(G1$line))]) +
  geom_smooth(method = 'lm',col='black',lwd = .5,alpha =.2) +
  xlab(bquote('Midparent FEC at 1st infection ('*sigma[p]*')')) +
  ylab(bquote('G1 offspring FEC at 1st infection ('*sigma[p]*')')) +
  theme(text = element_text(size = 14)) +
  geom_abline(slope=1, intercept=0,alpha = 0.4, lty = 3) +
  scale_x_continuous(limits = c(-3.5,3),breaks=seq(-3.5,3,1)) +
  scale_y_continuous(limits = c(-3.5,3),breaks=seq(-3.5,3,1))

psel20 = ggplot(G1,aes(x=midparent_fopg2,y=fs2)) +
  scale_color_manual(values = colsRS) +
  geom_point(size =3,alpha = .4,col = colsRS[match(G1$line,levels(G1$line))]) +
  geom_smooth(method = 'lm',col='black',lwd = .5,alpha =.2) +
  xlab(bquote('Midparent FEC at 2nd infection ('*sigma[p]*')')) +
  ylab(bquote('G1 offspring FEC at 2nd infection ('*sigma[p]*')')) +
  theme(text = element_text(size = 14)) +
  geom_abline(slope=1, intercept=0,alpha = 0.4, lty = 3) +
  scale_x_continuous(limits = c(-3.5,3),breaks=seq(-3.5,3,1)) +
  scale_y_continuous(limits = c(-3.5,3),breaks=seq(-3.5,3,1))

###---- Plot midparents vs. offspring G2
psel21 = ggplot(G2,aes(x=midparent_fopg1,y=fs1)) +
  scale_color_manual(values = colsRS[3:4]) +
  geom_point(size =3,alpha = .4,col = colsRS[2+match(G2$line,levels(G2$line))]) +
  geom_smooth(method = 'lm',col='black',lwd = .5,alpha =.2) +
  xlab(bquote('Midparent FEC at 1st infection ('*sigma[p]*')')) +
  ylab(bquote('G2 offspring FEC at 1st infection ('*sigma[p]*')')) +
  theme(text = element_text(size = 14)) +
  geom_abline(slope=1, intercept=0,alpha = 0.4, lty = 3) +
  scale_x_continuous(limits = c(-3.5,3),breaks=seq(-3.5,3,1)) +
  scale_y_continuous(limits = c(-3.5,3),breaks=seq(-3.5,3,1))

psel22 = ggplot(G2,aes(x=midparent_fopg2,y=fs2)) +
  geom_smooth(method = 'lm',col='black',lwd = .5,alpha =.2) +
  scale_color_manual(values = colsRS[3:4]) +
  geom_point(size = 3,alpha = .4,col = colsRS[2+match(G2$line,levels(G2$line))]) +
  xlab(bquote('Midparent FEC at 2nd infection ('*sigma[p]*')')) +
  ylab(bquote('G2 offspring FEC at 2nd infection ('*sigma[p]*')')) +
  theme(text = element_text(size = 14)) +
  geom_abline(slope=1, intercept=0,alpha = 0.4, lty = 3) +
  scale_x_continuous(limits = c(-3.5,3),breaks=seq(-3.5,3,1)) +
  scale_y_continuous(limits = c(-3.5,3),breaks=seq(-3.5,3,1))

#### Plot midparents vs. offsprings
pdf(file='supplementary_Figure1.pdf',width = 12,height = 10)
multiplot(psel1,psel2,
          psel10,psel21,
          psel20,psel22,cols=3)
dev.off()

###----- Regression on mid-parents by infection
### Expected divergence for G1 and G2 on average FEC across infections
#G1
summary(lm(fms ~ midparent_fopg,data = G1))
# Call:
#   lm(formula = fms ~ midparent_fopg, data = G1)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.68750 -0.90116  0.09951  0.77876  2.55009 
# 
# Coefficients:
#                Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      0.1237     0.1170   1.057    0.293    
# midparent_fopg   0.8578     0.1222   7.022 5.04e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.079 on 85 degrees of freedom
# Multiple R-squared:  0.3671,	Adjusted R-squared:  0.3596 
# F-statistic:  49.3 on 1 and 85 DF,  p-value: 5.037e-10

#G2
summary(lm(fms ~ midparent_fopg,data = G2))   
# Call:
#   lm(formula = fms ~ midparent_fopg, data = G2)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.63877 -0.37213  0.02968  0.41775  1.90964 
# 
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     -0.1121     0.1600  -0.700    0.488    
# midparent_fopg   0.5865     0.1076   5.451 3.47e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.9911 on 37 degrees of freedom
# Multiple R-squared:  0.4454,	Adjusted R-squared:  0.4304 
# F-statistic: 29.72 on 1 and 37 DF,  p-value: 3.473e-06

## By infection
summary(lm(fs1 ~ midparent_fopg1,data = G1))
# Call:
#   lm(formula = fs1 ~ midparent_fopg1, data = G1)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.98401 -0.80678  0.06622  0.77281  2.19214 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)      0.37043    0.11467   3.230  0.00176 ** 
#   midparent_fopg1  0.57524    0.07281   7.901 8.96e-12 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.9761 on 85 degrees of freedom
# Multiple R-squared:  0.4234,	Adjusted R-squared:  0.4166 
# F-statistic: 62.42 on 1 and 85 DF,  p-value: 8.959e-12

summary(lm(fs2 ~ midparent_fopg2,data = G1))
# Call:
#   lm(formula = fs2 ~ midparent_fopg2, data = G1)
# 
#      Min       1Q   Median       3Q      Max 
# -1.94157 -0.76898 -0.08049  0.92239  2.20674 
# 
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)   
# (Intercept)     -0.03138    0.11572  -0.271    0.787   
# midparent_fopg2  0.57055    0.16744   3.407    0.001 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.076 on 85 degrees of freedom
# Multiple R-squared:  0.1202,	Adjusted R-squared:  0.1098 
# F-statistic: 11.61 on 1 and 85 DF,  p-value: 0.001004

summary(lm(fs1 ~ midparent_fopg1,data = G2))
# Call:
#   lm(formula = fs1 ~ midparent_fopg1, data = G2)
# 
# Residuals:
# Min       1Q   Median       3Q      Max 
# -1.69281 -0.37932 -0.04355  0.57145  1.79644 
# 
# Coefficients:
#                Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     -0.01558    0.13718  -0.114     0.91    
# midparent_fopg1  0.40638    0.08735   4.652 4.11e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.8565 on 37 degrees of freedom
# Multiple R-squared:  0.3691,	Adjusted R-squared:  0.352 
# F-statistic: 21.65 on 1 and 37 DF,  p-value: 4.105e-05

summary(lm(fs2 ~ midparent_fopg2,data = G2))   
# Call:
#   lm(formula = fs2 ~ midparent_fopg2, data = G2)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.8581 -0.4873  0.1960  0.7594  2.0350 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)     -0.09201    0.20308  -0.453  0.65316   
# midparent_fopg2  0.66152    0.18886   3.503  0.00122 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.258 on 37 degrees of freedom
# Multiple R-squared:  0.249,	Adjusted R-squared:  0.2287 
# F-statistic: 12.27 on 1 and 37 DF,  p-value: 0.001222

#####------ Asymmetry of response
t.test(phenodiv$fms[phenodiv$Gene=='G1' & phenodiv$line=='S'],
       phenodiv$fms[phenodiv$Gene=='G2' & phenodiv$line=='S'])

# Welch Two Sample t-test
# data:  phenodiv$fms[phenodiv$Gene == "G1" & phenodiv$line == "S"] and phenodiv$fms[phenodiv$Gene == "G2" & phenodiv$line == "S"]
# t = 0.57963, df = 58.214, p-value = 0.5644
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.2664978  0.4837665
# sample estimates:
#   mean of x mean of y 
#    1.000159  0.891525 

summary(lm(fms ~ midparent_fopg:line,data = G1))
# Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            0.1232     0.2342   0.526 0.600179    
# midparent_fopg:lineR   0.8573     0.2465   3.478 0.000803 ***
# midparent_fopg:lineS   0.8585     0.3144   2.731 0.007697 ** 

summary(lm(fms ~ midparent_fopg:line,data = G2))
# Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)   
# (Intercept)           -0.3930     0.3598  -1.092  0.28187   
# midparent_fopg:lineR   0.3750     0.2653   1.414  0.16604   
# midparent_fopg:lineS   0.7966     0.2638   3.019  0.00464 **
summary.aov(lm(fms ~ midparent_fopg:line,data = G2))
#                     Df Sum Sq Mean Sq F value   Pr(>F)    
# midparent_fopg:line  2  29.95  14.974   15.14 1.69e-05 ***
# Residuals           36  35.59   0.989

summary(lm(fs1 ~ midparent_fopg1:line,data = G1))
# Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             0.3389     0.2449   1.384   0.1701    
# midparent_fopg1:lineR   0.5589     0.1335   4.187 6.94e-05 ***
# midparent_fopg1:lineS   0.6200     0.3151   1.968   0.0524 .  

summary(lm(fs2 ~ midparent_fopg2:line,data = G1))
# Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            -0.2537     0.1565  -1.621 0.108804    
# midparent_fopg2:lineR   0.1852     0.2487   0.745 0.458495    
# midparent_fopg2:lineS   1.1245     0.3147   3.573 0.000587 ***

summary(lm(fs1 ~ midparent_fopg1:line,data = G2))
# Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)   
# (Intercept)           -0.55599    0.47953  -1.159   0.2539  
# midparent_fopg1:lineR  0.04698    0.31783   0.148   0.8833  
# midparent_fopg1:lineS  0.76486    0.31707   2.412   0.0211 *
  
summary(lm(fs2 ~ midparent_fopg2:line,data = G2))   
# Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)   
# (Intercept)            -0.4691     0.2867  -1.636  0.11051   
# midparent_fopg2:lineR   0.2882     0.2758   1.045  0.30313   
# midparent_fopg2:lineS   1.2544     0.3752   3.344  0.00194 **
  

### Check inbreeding within line and generation
inbreeding = inbreeding[inbreeding$ipg %in% c(G1$IPG[G1$pheno=='OUI'],G2$IPG),]

# Inbreeding within produced lambs
aggregate(coef ~ line + Gene,data=inbreeding,FUN=summary)
#   line Gene   coef.Min. coef.1st Qu. coef.Median   coef.Mean coef.3rd Qu.   coef.Max.
# 1    R _G1_ 0.000000000  0.000000000 0.000000000 0.029891304  0.000000000 0.250000000
# 2    S _G1_ 0.000000000  0.000000000 0.000000000 0.006097561  0.000000000 0.125000000
# 3    R G2_1 0.000000000  0.000000000 0.000000000 0.003472222  0.000000000 0.062500000
# 4    S G2_1 0.000000000  0.000000000 0.000000000 0.014045000  0.031200000 0.062500000

# Inbreeding within selected lambs
table(inbreeding$Gene,inbreeding$line)
#       R  S
# _G1_ 46 41
# _G3_  0  0
# G2_1 18 20
# G2_2  0  0

a = aggregate(coef ~ line + Gene,data=inbreeding,FUN=summary)
a$std = aggregate(coef ~ line + Gene,data=inbreeding,FUN=sd)[,3]
a
#   line Gene   coef.Min. coef.1st Qu. coef.Median   coef.Mean coef.3rd Qu.   coef.Max.        std
# 1    R _G1_ 0.000000000  0.000000000 0.000000000 0.029891304  0.000000000 0.250000000 0.07539149
# 2    S _G1_ 0.000000000  0.000000000 0.000000000 0.006097561  0.000000000 0.125000000 0.02726060
# 3    R G2_1 0.000000000  0.000000000 0.000000000 0.003289474  0.000000000 0.062500000 0.01433848
# 4    S G2_1 0.000000000  0.000000000 0.000000000 0.014045000  0.031200000 0.062500000 0.01888342

t.test(coef~ line,data = inbreeding[inbreeding$Gene=='G2_1',])

# Welch Two Sample t-test
# 
# data:  coef by line
# t = -2.0094, df = 35.326, p-value = 0.05218
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.0216182015  0.0001071489
# sample estimates:
#   mean in group R mean in group S 
# 0.003289474     0.014045000 

## Reduced inbreeding in the susceptible line

###===============================================================================================
## Divergent lines performances: 
## individual trajectories of selected individuals after H.contortus infection
###===============================================================================================

####------=========== Modeling resistance vs. H.contortus G1 & G2 data ======----------#####

####------ FEC modeling G1: without the stressed individuals
fec$time = factor(paste(fec$inf,fec$day,sep="-"))
fec$mod_al = factor(fec$mod_al)

## Plot raw data per line & time-point
pfec_g1 = ggplot(fec,aes(x = time,y = value + 1, fill = gen)) +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
               geom = "crossbar", width = 0.4, color = rep(c(colsRS[2],colsRS[1]),4),
               aes(color = gen), 
               position = position_dodge(width = .4)) +
  #geom_dotplot(binaxis = 'y', stackdir = 'center', position = position_dodge(.4), binwidth = 0.08) + 
  geom_boxplot() +
  scale_fill_manual(values=colsRS[1:2]) +
  annotation_logticks(scaled = TRUE,sides="l") +
  scale_y_log10(breaks = c(0.000000001, 10^(-2:5)),
                labels = c(0, math_format()(-2:5))) +
  xlab('Time point (Infection - Day post-infection)') + ylab('Faecal Egg Count (eggs/g)') +
  theme(legend.position = 'none',text= element_text(size=14),axis.text.x = element_text(angle=45, vjust=0.5))

shapiro.test(sqrt(sqrt(fec$value)))
# data:  sqrt(sqrt(fec$value))
# W = 0.93339, p-value = 3.819e-07

shapiro.test(log(50+fec$value))
# data:  log(50 + fec$value)
# W = 0.90738, p-value = 6.102e-09

m1 = lme(sqrt(sqrt(value)) ~ sexe + mod_al + gen*time, 
         random =~ 1|IPG,
        data = fec)

summary(m1)
# Linear mixed-effects model fit by REML
# Data: fec 
#      AIC      BIC    logLik
# 817.7307 854.7819 -396.8654
# 
# Random effects:
#   Formula: ~1 | IPG
# (Intercept) Residual
# StdDev:    1.357368 2.285629
# 
# Fixed effects: sqrt(sqrt(value)) ~ sexe + mod_al + gen * time 
#                      Value Std.Error  DF   t-value p-value
# (Intercept)       3.524762 0.6658832 123  5.293364  0.0000
# sexe2            -0.009259 0.5432116  39 -0.017045  0.9865
# mod_al2          -1.245452 0.5983132  39 -2.081605  0.0440
# genS-G1           3.782718 0.8247270  39  4.586630  0.0000
# time1-30          1.614800 0.6739953 123  2.395863  0.0181
# time2-24         -1.196111 0.6739953 123 -1.774658  0.0784
# time2-30          0.574258 0.6739953 123  0.852021  0.3959
# genS-G1:time1-30 -0.428413 0.9882711 123 -0.433497  0.6654
# genS-G1:time2-24 -2.208654 0.9882711 123 -2.234866  0.0272
# genS-G1:time2-30 -0.185993 0.9882711 123 -0.188201  0.8510
#  Correlation: 
#                  (Intr) sexe2  mod_l2 gnS-G1 tm1-30 tm2-24 tm2-30 gS-G1:1 gS-G1:2-2
# sexe2            -0.394                                                            
# mod_al2          -0.357 -0.080                                                     
# genS-G1          -0.626  0.001  0.169                                              
# time1-30         -0.506  0.000  0.000  0.409                                       
# time2-24         -0.506  0.000  0.000  0.409  0.500                                
# time2-30         -0.506  0.000  0.000  0.409  0.500  0.500                         
# genS-G1:time1-30  0.345  0.000  0.000 -0.599 -0.682 -0.341 -0.341                  
# genS-G1:time2-24  0.345  0.000  0.000 -0.599 -0.341 -0.682 -0.341  0.500           
# genS-G1:time2-30  0.345  0.000  0.000 -0.599 -0.341 -0.341 -0.682  0.500   0.500   
# 
# Standardized Within-Group Residuals:
#        Min         Q1        Med         Q3        Max 
# -2.3232027 -0.6342790  0.1155742  0.6224571  2.2581325 
# 
# Number of Observations: 172
# Number of Groups: 43

##-- Diff in phen. std
3.782718/sd(sqrt(sqrt(fec$value)))
#[1] 1.101783

####------ FEC modeling G2
## Plot raw data per line & time-point
pfec_g2 = ggplot(fec2,aes(x = time,y = value+1, fill = lot)) +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
               geom = "crossbar", width = 0.4, color = rep(c(colsRS[4],colsRS[3]),3),
               aes(color = lot), 
               position = position_dodge(width = .4)) +
  #geom_dotplot(binaxis = 'y', stackdir = 'center', position = position_dodge(.4), binwidth = 0.08) + 
  geom_boxplot() +
  scale_fill_manual(values = colsRS[3:4]) +
  annotation_logticks(scaled = TRUE,sides="l") +
  scale_y_log10(breaks = c(0.000000001, 10^(-2:5)),
                labels = c(0, math_format()(-2:5))) +
  xlab('Time point (Infection - Day post-infection)') + ylab('Faecal Egg Count (eggs/g)') +
  theme(legend.position = 'none',text= element_text(size=14),axis.text.x = element_text(angle=45, vjust=0.5))

fec2$time = factor(paste(fec2$inf,fec2$day,sep="-"))
fec2$mod_al = factor(fec2$mod_al) ## not considered further as all 0

shapiro.test(log(50+fec2$value))
# data:  log(50 + fec2$value)
# W = 0.8644, p-value = 8.387e-09

shapiro.test(sqrt(sqrt(fec2$value)))
# data:  log(50 + fec2$value)
# W = 0.88026, p-value = 4.061e-08

plot(density(sqrt(sqrt(fec2$value))))
plot(density(log(50+fec2$value)))

m2 = lme(sqrt(sqrt(value)) ~ sexe + gen*time,
         random =~ 1|IPG,
          data = fec2)
summary(m2)
# Linear mixed-effects model fit by REML
# Data: fec2 
# AIC      BIC    logLik
# 512.2509 536.3063 -247.1254
# 
# Random effects:
#   Formula: ~1 | IPG
# (Intercept) Residual
# StdDev:    1.075184 1.993445
# 
# Fixed effects: sqrt(sqrt(value)) ~ sexe + gen + time + gen:time 
#                      Value Std.Error DF   t-value p-value
# (Intercept)       1.126750 0.5918370 72  1.903819  0.0609
# sexe2            -1.124497 0.5109978 35 -2.200591  0.0345
# genS-G2           4.076465 0.7358557 35  5.539762  0.0000
# time1-30          3.143682 0.6644817 72  4.731029  0.0000
# time2-30          3.319297 0.6644817 72  4.995317  0.0000
# genS-G2:time1-30 -0.769951 0.9159248 72 -0.840626  0.4033
# genS-G2:time2-30  0.624081 0.9159248 72  0.681367  0.4978
# Correlation: 
#                   (Intr) sexe2  gnS-G2 tm1-30 tm2-30 gS-G2:1
# sexe2            -0.432                                    
# genS-G2          -0.654  0.000                             
# time1-30         -0.561  0.000  0.452                      
# time2-30         -0.561  0.000  0.452  0.500               
# genS-G2:time1-30  0.407  0.000 -0.622 -0.725 -0.363        
# genS-G2:time2-30  0.407  0.000 -0.622 -0.363 -0.725  0.500 
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -1.91943478 -0.48364489  0.06600518  0.50302581  2.10111636 
# 
# Number of Observations: 114
# Number of Groups: 38

##-- Diff in phen. std
4.076465/sd(sqrt(sqrt(fec2$value)))
#[1] 1.189584

##-- Dispersion of FEC values in G2 at 2-30 dpi
summary(fec2$value[fec2$time=='2-30' & fec2$gen=='S-G2'])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1650    3962    5625    5918    6825   14250 

####------ WEIGHT modeling G1 // units is hectogram
wgt$w0 = wgt$value[match(wgt$IPG,wgt$IPG[wgt$variable=='pds011'])]
wgt$time = factor(paste(wgt$inf,wgt$day,sep="-"))
wgt$mod_al = factor(wgt$mod_al) 

shapiro.test(wgt$value[wgt$variable!='pds011'])
# data:  wgt$value[wgt$variable != "pds011"]
# W = 0.98453, p-value = 0.002506

m3 = lme(value/10 ~ mod_al + sexe + gen + time + time:gen + w0,
         random =~ 1|IPG,
         data =  wgt[wgt$variable!='pds011',])
anova.lme(m3)
#             numDF denDF   F-value p-value
# (Intercept)     1   246 22508.750  <.0001
# mod_al          1    38    77.875  <.0001
# sexe            1    38   372.627  <.0001
# gen             1    38    26.390  <.0001
# time            6   246   205.603  <.0001
# w0              1    38   272.660  <.0001
# gen:time        6   246     1.536  0.1669

summary(m3,type='marginal')
# Linear mixed-effects model fit by REML
# Data: wgt[wgt$variable != "pds011", ] 
# AIC      BIC    logLik
# 1445.404 1514.735 -703.7021
# 
# Random effects:
#   Formula: ~1 | IPG
# (Intercept) Residual
# StdDev:    1.883985  2.28681
# 
# Fixed effects: value/10 ~ mod_al + sexe + gen + time + time:gen + w0 
#                      Value Std.Error  DF   t-value p-value
# (Intercept)      10.113915 2.3827505 246  4.244639  0.0000
# mod_al2          -1.340479 0.7484441  38 -1.791020  0.0813
# sexe2            -5.958042 0.7402545  38 -8.048640  0.0000
# genS-G1          -0.568358 0.9242645  38 -0.614930  0.5423
# time1-24          6.613043 0.6743434 246  9.806640  0.0000
# time1-30          7.456522 0.6743434 246 11.057454  0.0000
# time2-0          11.317391 0.6743434 246 16.782830  0.0000
# time2-14         14.095652 0.6743434 246 20.902779  0.0000
# time2-24         13.913043 0.6743434 246 20.631985  0.0000
# time2-30         14.443478 0.6743434 246 21.418579  0.0000
# w0                0.091888 0.0055648  38 16.512422  0.0000
# genS-G1:time1-24 -1.103043 0.9887816 246 -1.115558  0.2657
# genS-G1:time1-30 -0.901522 0.9887816 246 -0.911750  0.3628
# genS-G1:time2-0  -2.312391 0.9887816 246 -2.338627  0.0202
# genS-G1:time2-14 -2.515652 0.9887816 246 -2.544194  0.0116
# genS-G1:time2-24 -1.768043 0.9887816 246 -1.788103  0.0750
# genS-G1:time2-30 -1.373478 0.9887816 246 -1.389061  0.1661
# Correlation: 
#                   (Intr) mod_l2 sexe2  gnS-G1 tm1-24 tm1-30 tim2-0 tm2-14 tm2-24 tm2-30 w0     gS-G1:1-2 gS-G1:1-3 gS-G1:2-0 gS-G1:2-1
# mod_al2          -0.448                                                                                                              
# sexe2            -0.599  0.120                                                                                                       
# genS-G1          -0.281  0.197  0.047                                                                                                
# time1-24         -0.142  0.000  0.000  0.365                                                                                         
# time1-30         -0.142  0.000  0.000  0.365  0.500                                                                                  
# time2-0          -0.142  0.000  0.000  0.365  0.500  0.500                                                                           
# time2-14         -0.142  0.000  0.000  0.365  0.500  0.500  0.500                                                                    
# time2-24         -0.142  0.000  0.000  0.365  0.500  0.500  0.500  0.500                                                             
# time2-30         -0.142  0.000  0.000  0.365  0.500  0.500  0.500  0.500  0.500                                                      
# w0               -0.949  0.358  0.515  0.090  0.000  0.000  0.000  0.000  0.000  0.000                                               
# genS-G1:time1-24  0.097  0.000  0.000 -0.535 -0.682 -0.341 -0.341 -0.341 -0.341 -0.341  0.000                                        
# genS-G1:time1-30  0.097  0.000  0.000 -0.535 -0.341 -0.682 -0.341 -0.341 -0.341 -0.341  0.000  0.500                                 
# genS-G1:time2-0   0.097  0.000  0.000 -0.535 -0.341 -0.341 -0.682 -0.341 -0.341 -0.341  0.000  0.500     0.500                       
# genS-G1:time2-14  0.097  0.000  0.000 -0.535 -0.341 -0.341 -0.341 -0.682 -0.341 -0.341  0.000  0.500     0.500     0.500             
# genS-G1:time2-24  0.097  0.000  0.000 -0.535 -0.341 -0.341 -0.341 -0.341 -0.682 -0.341  0.000  0.500     0.500     0.500     0.500   
# genS-G1:time2-30  0.097  0.000  0.000 -0.535 -0.341 -0.341 -0.341 -0.341 -0.341 -0.682  0.000  0.500     0.500     0.500     0.500   
# gS-G1:2-2
# mod_al2                   
# sexe2                     
# genS-G1                   
# time1-24                  
# time1-30                  
# time2-0                   
# time2-14                  
# time2-24                  
# time2-30                  
# w0                        
# genS-G1:time1-24          
# genS-G1:time1-30          
# genS-G1:time2-0           
# genS-G1:time2-14          
# genS-G1:time2-24          
# genS-G1:time2-30  0.500   
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -2.64396761 -0.60677078  0.02451746  0.59750142  3.03765585 
# 
# Number of Observations: 301
# Number of Groups: 43  

pwgt_g1 = ggplot(wgt[wgt$variable != 'pds011',],
                 aes(x = time,y = (value - w0)/10, fill = lot)) +
  # stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
  #              geom = "crossbar", width = 0.4, color = rep(c(colsRS[2],colsRS[1]),7),
  #              aes(color = lot),
  #              position = position_dodge(width = .6)) +
  #geom_dotplot(binaxis = 'y', stackdir = 'center', position = position_dodge(.6), binwidth = 0.4) +
  geom_boxplot(alpha = .5) +
  scale_fill_manual(values = colsRS[1:2]) +
  scale_y_continuous(breaks = seq(0,30,5),limits = c(0,30)) +
  xlab('Time point (Infection - Day post-infection)') + 
  ylab('Weight gain (Kg)') +
  theme(legend.position = 'none',text= element_text(size=16),
        axis.text.x = element_text(angle=45, vjust=0.5))

####------ WEIGHT modeling G2

## Significant difference in weight before infection betw. both lines
summary.aov(lm(value/10 ~ sexe + gen,data = wgt2[wgt2$variable=='pds011',]))
#             Df Sum Sq Mean Sq F value   Pr(>F)    
# sexe         1  507.1   507.1  26.292 1.09e-05 ***
# gen          1  173.3   173.3   8.987  0.00498 ** 
# Residuals   35  675.0    19.3 

summary(lm(value/10 ~ sexe + gen,data = wgt2[wgt2$variable=='pds011',]))
# Call:
#   lm(formula = value/10 ~ sexe + gen, data = wgt2[wgt2$variable == 
#                                                     "pds011", ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -9.1446 -1.8296  0.4276  2.1724 10.4081 
# 
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)   36.955      1.257  29.409  < 2e-16 ***
#   sexe2         -7.306      1.425  -5.128 1.09e-05 ***
#   genS-G2        4.277      1.427   2.998  0.00498 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 4.392 on 35 degrees of freedom
# Multiple R-squared:  0.502,	Adjusted R-squared:  0.4735 
# F-statistic: 17.64 on 2 and 35 DF,  p-value: 5.032e-06

## Correction using weight at d0 as a covariate
wgt2$w0 = wgt2$value[match(wgt2$IPG,wgt2$IPG[wgt2$variable=='pds011'])]
wgt2$time = factor(paste(wgt2$inf,wgt2$day,sep="-"))
wgt2$mod_al = factor(wgt2$mod_al) 

shapiro.test(wgt2$value[wgt2$variable!='pds011'])
# data:  wgt2$value[wgt2$variable != "pds011"]
# W = 0.99099, p-value = 0.2819

m4 = lme(value/10 ~ sexe + gen + time + time:gen + w0,
         random =~ 1|IPG,
         data =  wgt2[wgt2$variable!='pds011',])
summary(m4)
# Linear mixed-effects model fit by REML
# Data: wgt2[wgt2$variable != "pds011", ] 
# AIC     BIC    logLik
# 885.383 929.928 -428.6915
# 
# Random effects:
#   Formula: ~1 | IPG
# (Intercept) Residual
# StdDev:    1.938526 2.002001
# 
# Fixed effects: value/10 ~ sexe + gen + time + time:gen + w0 
#                      Value Std.Error  DF   t-value p-value
# (Intercept)       1.427604 3.1266482 144  0.456593  0.6487
# sexe2            -1.455738 0.9167745  34 -1.587891  0.1216
# genS-G2          -1.276065 0.9712435  34 -1.313846  0.1977
# time1-30          5.478889 0.6673335 144  8.210121  0.0000
# time2-0           5.846667 0.6673335 144  8.761236  0.0000
# time2-14          7.177778 0.6673335 144 10.755907  0.0000
# time2-30         11.683889 0.6673335 144 17.508319  0.0000
# w0                0.114204 0.0082187  34 13.895652  0.0000
# genS-G2:time1-30 -0.332389 0.9198558 144 -0.361349  0.7184
# genS-G2:time2-0   0.467333 0.9198558 144  0.508051  0.6122
# genS-G2:time2-14  0.520222 0.9198558 144  0.565548  0.5726
# genS-G2:time2-30 -1.469389 0.9198558 144 -1.597412  0.1124
# Correlation: 
#   (Intr) sexe2  gnS-G2 tm1-30 tim2-0 tm2-14 tm2-30 w0     gS-G2:1 gS-G2:2-0 gS-G2:2-1
# sexe2            -0.720                                                                             
# genS-G2           0.210 -0.237                                                                      
# time1-30         -0.107  0.000  0.344                                                               
# time2-0          -0.107  0.000  0.344  0.500                                                        
# time2-14         -0.107  0.000  0.344  0.500  0.500                                                 
# time2-30         -0.107  0.000  0.344  0.500  0.500  0.500                                          
# w0               -0.971  0.655 -0.362  0.000  0.000  0.000  0.000                                   
# genS-G2:time1-30  0.077  0.000 -0.474 -0.725 -0.363 -0.363 -0.363  0.000                            
# genS-G2:time2-0   0.077  0.000 -0.474 -0.363 -0.725 -0.363 -0.363  0.000  0.500                     
# genS-G2:time2-14  0.077  0.000 -0.474 -0.363 -0.363 -0.725 -0.363  0.000  0.500   0.500             
# genS-G2:time2-30  0.077  0.000 -0.474 -0.363 -0.363 -0.363 -0.725  0.000  0.500   0.500     0.500   
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -4.09082162 -0.50279373 -0.03925336  0.55357964  2.21729083 
# 
# Number of Observations: 190
# Number of Groups: 38 

pwgt_g2 = ggplot(wgt2[wgt2$variable!='pds011',],aes(x = time,y = (value - w0)/10, fill = lot)) +
  # stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
  #              geom = "crossbar", width = 0.5, color = rep(c(colsRS[4],colsRS[3]),5),
  #              aes(color = lot), 
  #              position = position_dodge(width = .6)) +
  #geom_dotplot(binaxis = 'y', stackdir = 'center', position = position_dodge(.6), binwidth = 0.4) + 
  geom_boxplot(alpha = .5) +
  scale_fill_manual(values=colsRS[3:4]) +
  scale_y_continuous(breaks = seq(0,30,5),limits = c(0,30)) +
  xlab('Time point (Infection - Day post-infection)') + 
  ylab('Weight gain (Kg)') +
  theme(legend.position = 'none',text= element_text(size=16),
        axis.text.x = element_text(angle=45, vjust=0.5))

## Supplementary Figure 2
pdf(file = paste0(draft,'supplementary_Figure2.pdf'),width = 14,height = 10)
multiplot(pwgt_g1,pwgt_g2,
          cols=2)
dev.off()

########---------------------------===== RESISTANCE potential across environments  =====--------------------------------------########

###---- G1 generation: stress x genotype

########---------------------------===== BEHAVIOUR data following stress treatment =====--------------------------------------########
behav = read.csv(file = './data/BEHAVIOUR/behaviour_data.csv',header=T,sep=';')
behav$IPG = 20000150000+behav$boucle
behav$boucle = NULL

## Put standing-lying recordings into a df (hobo system)
hobo = behav[,c(1:3,ncol(behav),grep('Standing',colnames(behav)))]
hobo = melt(hobo,c('sexe','traitement','genotype','IPG'))
hobo$day = factor(sapply(str_split(hobo$variable,'_'),function(x) x[1]))
hobo$var = factor(sapply(str_split(hobo$variable,'_'),function(x) x[2]))
hobo$d = 'd1'
hobo$d[hobo$day=='ATBI'] = 'd2'
hobo$d[hobo$day=='D14I1'] = 'd3'
hobo$d[hobo$day=='D14I2'] = 'd4'
hobo$d = factor(hobo$d)
hobo$variable = NULL
##Rename columns
colnames(hobo) = c('sex','treatment','line','ind','value','day','var','d')
hobo = na.omit(hobo) ## measured on females only

dim(hobo)
#[1] 320   8

## Summary stats - Note: value has to be divided by 5 to yield time in min/day
stats = aggregate(value/5 ~ var + line + treatment, FUN=mean, data=hobo)
stats$std = aggregate(value/5 ~ var + line + treatment, FUN=sd,data = hobo)[,4]
stats
#               var line treatment value/5       std
# 1 sumFreqStanding    R     nostr  18.385  4.958575
# 2     sumStanding    R     nostr 630.370 63.228954
# 3 sumFreqStanding    S     nostr  20.665  5.200274
# 4     sumStanding    S     nostr 638.145 62.216363
# 5 sumFreqStanding    R       str  20.715  5.296131
# 6     sumStanding    R       str 587.380 62.873578
# 7 sumFreqStanding    S       str  21.210  5.380678
# 8     sumStanding    S       str 627.345 93.476611

stats.d0 = aggregate(value/5 ~ var + line + treatment, FUN = mean, data = hobo[hobo$d=='d1',])
stats.d0$std = aggregate(value/5 ~ var + line + treatment, FUN = sd,data = hobo[hobo$d=='d1',])[,4]
stats.d0
#               var line treatment value/5       std
# 1 sumFreqStanding    R     nostr   23.36  4.831195
# 2     sumStanding    R     nostr  627.88 67.254290
# 3 sumFreqStanding    S     nostr   26.06  4.745103
# 4     sumStanding    S     nostr  616.36 57.515973
# 5 sumFreqStanding    R       str   24.72  5.818896
# 6     sumStanding    R       str  578.12 49.239027
# 7 sumFreqStanding    S       str   23.16  4.553680
# 8     sumStanding    S       str  573.96 53.547802

## Plot density distribution
ggplot(hobo[hobo$var=='sumFreqStanding',],aes(x = value)) + 
  geom_density(alpha =.3) +
  theme(legend.position = 'bottom')

ggplot(hobo[hobo$var=='sumStanding',],aes(x = value)) + 
  geom_density(alpha =.3) +
  theme(legend.position = 'bottom')

shapiro.test(hobo$value[hobo$var=='sumFreqStanding'])
# data:  hobo$value[hobo$var == "sumFreqStanding"]
# W = 0.97151, p-value = 0.002171
shapiro.test(hobo$value[hobo$var=='sumStanding'])
# data:  hobo$value[hobo$var == "sumStanding"]
# W = 0.89977, p-value = 5.537e-09

summary(hobo$value[hobo$var=='sumStanding']/5)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 360.0   576.9   619.5   620.8   660.3  1087.2  

####-----===== Test standing-lying before behaviour treatment
bh$sex = factor(bh$sex)

## Evaluate differences before treatment
m0.h1 = lm(value/5 ~ line*treatment,
            data = hobo[hobo$var=='sumStanding' & hobo$d=='d1',])
anova(m0.h1)
# Response: value/5
#                Df Sum Sq Mean Sq F value  Pr(>F)  
# line            1    615   614.7  0.1874 0.66771  
# treatment       1  21234 21233.7  6.4722 0.01539 *
# line:treatment  1    135   135.4  0.0413 0.84015  
# Residuals      36 118108  3280.8  

summary(m0.h1, type = 'marginal')
# lm(formula = value/5 ~ line * treatment, data = hobo[hobo$var == 
#                                                        "sumStanding" & hobo$d == "d1", ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -108.56  -50.49    6.04   42.19  131.92 
# 
# Coefficients:
#                    Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          627.88      18.11  34.665   <2e-16 ***
# lineS                -11.52      25.61  -0.450   0.6556    
# treatmentstr         -49.76      25.61  -1.943   0.0599 .  
# lineS:treatmentstr     7.36      36.23   0.203   0.8401    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 57.28 on 36 degrees of freedom
# Multiple R-squared:  0.1569,	Adjusted R-squared:  0.08667 
# F-statistic: 2.234 on 3 and 36 DF,  p-value: 0.101

m0.h2 = lm(value/5 ~ line*treatment,
            data = hobo[hobo$var=='sumFreqStanding' & hobo$d=='d1',])
anova(m0.h2)
# Response: value/5
#                 Df Sum Sq Mean Sq F value Pr(>F)
# line            1   3.25   3.249  0.1294 0.7212
# treatment       1   5.93   5.929  0.2361 0.6300
# line:treatment  1  45.37  45.369  1.8066 0.1873
# Residuals      36 904.07  25.113

summary(m0.h2)
# Coefficients:
#                    Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          23.360      1.585  14.741   <2e-16 ***
# lineS                 2.700      2.241   1.205    0.236    
# treatmentstr          1.360      2.241   0.607    0.548    
# lineS:treatmentstr   -4.260      3.169  -1.344    0.187 

####-----===== Test standing-lying throughout experiment

### SumStanding not considered further

### SumFreqStanding
mfull = lm(value ~ d*line*treatment,
           data = hobo[hobo$var=='sumFreqStanding',])
s = stepAIC(mfull,direction = 'both')
s$anova
# Final Model:
#   value ~ d + line + treatment
rm(mfull,s)

m.hfs = lme(value ~  d + line + treatment,
           random =~ 1|ind,
           data = hobo[hobo$var=='sumFreqStanding',])
anova.lme(m.hfs)
#             numDF denDF   F-value p-value
# (Intercept)     1   117 1598.8512  <.0001
# d               3   117   45.6497  <.0001
# line            1    37    1.8777  0.1789
# treatment       1    37    2.0155  0.1641

summary(m.hfs,type='marginal')
# Linear mixed-effects model fit by REML
# Data: hobo[hobo$var == "sumFreqStanding", ] 
# AIC      BIC    logLik
# 1387.051 1411.347 -685.5256
# 
# Random effects:
#   Formula: ~1 | ind
# (Intercept) Residual
# StdDev:    13.73343 16.45694
# 
# Fixed effects: value ~ d + line + treatment 
# Value Std.Error  DF    t-value p-value
# (Intercept)  114.5625  4.929675 117  23.239361  0.0000
# dd2          -16.6250  3.679884 117  -4.517806  0.0000
# dd3          -22.4000  3.679884 117  -6.087149  0.0000
# dd4          -42.6000  3.679884 117 -11.576453  0.0000
# lineS          6.9375  5.062755  37   1.370301  0.1789
# treatmentstr   7.1875  5.062755  37   1.419681  0.1641

####-----===== Other behavioural traits
behav = behav[,-c(grep('Standing',colnames(behav)))]
pc = ade4::dudi.mix(behav[,-ncol(behav)],nf=2,scannf=F)
ade4::s.corcircle(pc$co)

## Plot and test
bh = melt(behav,c('sexe','traitement','genotype','IPG'))
bh$day = factor(sapply(str_split(bh$variable,'_'),function(x) x[1]))
bh$var = factor(sapply(str_split(bh$variable,'_'),function(x) x[2]))
bh$d = 'd1'
bh$d[bh$day=='ATBI'] = 'd2'
bh$d[bh$day=='D14I1'] = 'd3'
bh$d[bh$day=='D14I2'] = 'd4'
bh$d = factor(bh$d)
bh$variable = NULL

## Rename columns in english
colnames(bh) = c('sex','treatment','line','ind','value','day','var','d')

## Summary stats
stats = aggregate(value ~ var + line + treatment, FUN=mean, data=bh)
stats$std = aggregate(value ~ var + line + treatment, FUN=sd,data = bh)[,4]
stats
#            var line treatment      value       std
# 1       bleat1    R     nostr  6.8977273 5.0923422
# 2       bleat2    R     nostr  5.4659091 4.2883131
# 3     contact2    R     nostr  0.5113636 0.9221953
# 4  locomotion1    R     nostr 18.0795455 8.2213929
# 5  locomotion2    R     nostr 11.2500000 6.3828974
# 6   vigilance1    R     nostr  7.5454545 2.4349424
# 7   vigilance2    R     nostr  7.3295455 2.3030294
# 8       bleat1    S     nostr  7.4000000 5.4970648
# 9       bleat2    S     nostr  5.6625000 4.5477599
# 10    contact2    S     nostr  0.2875000 0.6970780
# 11 locomotion1    S     nostr 17.7625000 7.9480015
# 12 locomotion2    S     nostr 10.2875000 5.8139809
# 13  vigilance1    S     nostr  7.1000000 2.0598636
# 14  vigilance2    S     nostr  7.4875000 2.3329980
# 15      bleat1    R       str  4.5454545 4.0624694
# 16      bleat2    R       str  3.1136364 3.5442419
# 17    contact2    R       str  0.4431818 1.0705886
# 18 locomotion1    R       str 14.9545455 7.3997656
# 19 locomotion2    R       str 10.0454545 5.6505703
# 20  vigilance1    R       str  7.1931818 2.6386079
# 21  vigilance2    R       str  6.3409091 2.6647394
# 22      bleat1    S       str  7.7375000 5.4950927
# 23      bleat2    S       str  5.0250000 4.4179582
# 24    contact2    S       str  0.2125000 0.5668869
# 25 locomotion1    S       str 19.9125000 9.2383265
# 26 locomotion2    S       str 10.9000000 7.0128273
# 27  vigilance1    S       str  6.7125000 2.3876992
# 28  vigilance2    S       str  7.1250000 2.6972795

## Plot density distribution
ggplot(bh,aes(x = value)) + 
  geom_density(alpha =.3) +
  facet_wrap(~var, scales = "free_y") +
  theme(legend.position = 'bottom')

##--- Estimate Shapiro Wilk by variable
shapiro.test(bh$value[bh$var=='bleat1'])
#W = 0.93964, p-value = 1.86e-10
plot(density(bh$value[bh$var=='bleat1']))

shapiro.test(bh$value[bh$var=='bleat2']) 
#W = 0.90801, p-value = 1.936e-13
plot(density(bh$value[bh$var=='bleat2']))
shapiro.test(sqrt(bh$value[bh$var=='bleat2']))
#W = 0.92496, p-value = 6.044e-12

shapiro.test(bh$value[bh$var=='contact2']) ## too many 0s // NB or logistic reg.
#W = 0.49989, p-value < 2.2e-16
table(bh$value[bh$var=='contact2'])
#  0   1   2   3   4   5   6 
#262  44  18   7   3   1   1 

require(fitdistrplus)
qqcomp(fitdist(data = bh$value[bh$var=='contact2'], distr = "norm", method = "mle"),main="Contact2 \n QQ-Plot - Normal")
qqcomp(fitdist(data = bh$value[bh$var=='contact2'], distr = "pois", method = "mle"),main="Contact2 \n QQ-Plot - Poisson")
qqcomp(fitdist(data = bh$value[bh$var=='contact2'], distr = "nbinom", method = "mle"),main="Contact2 \n QQ-Plot - Negative binomial")

shapiro.test(bh$value[bh$var=='locomotion1'])
#W = 0.96166, p-value = 1.026e-07
shapiro.test(bh$value[bh$var=='locomotion2'])
#W = 0.91323, p-value = 5.318e-13
shapiro.test(bh$value[bh$var=='vigilance1'])
#W = 0.97806, p-value = 5.297e-05
shapiro.test(bh$value[bh$var=='vigilance2'])
#W = 0.98202, p-value = 0.0003274

###------====== Behaviour data modelling
bh$sex = factor(bh$sex)

## Evaluate differences before treatment
m0.b1 = lme(value ~ sex*treatment + line*treatment,
            random =~ 1|ind,
            data = bh[bh$var=='bleat1' & bh$d=='d1',])
Anova(m0.b1,'III')
# Response: value
#                 Chisq Df Pr(>Chisq)   
# (Intercept)    5.6866  1    0.01710 * 
# sex            2.6185  1    0.10563   
# treatment      1.5148  1    0.21841   
# line           0.0548  1    0.81496   
# sex:treatment  0.5475  1    0.45933   
# treatment:line 7.2124  1    0.00724 **

summary(m0.b1)
# Linear mixed-effects model fit by REML
# Data: bh[bh$var == "bleat1" & bh$d == "d1", ] 
# AIC     BIC    logLik
# 478.2743 497.128 -231.1372
# 
# Random effects:
#   Formula: ~1 | ind
# (Intercept) Residual
# StdDev:      3.9371 1.476412
# 
# Fixed effects: value ~ sex * treatment + line * treatment 
#                        Value Std.Error DF    t-value p-value
# (Intercept)         5.247826  2.200666 78  2.3846533  0.0195
# sex                 2.104348  1.300455 78  1.6181630  0.1097
# treatmentstr       -3.830435  3.112212 78 -1.2307756  0.2221
# lineS              -0.304348  1.300455 78 -0.2340318  0.8156
# sex:treatmentstr    1.360870  1.839121 78  0.7399566  0.4615
# treatmentstr:lineS  4.939130  1.839121 78  2.6855932  0.0088

## Bleat1 is more important in susceptible individuals under chronic stress before treatment
## This should be accounted for in subsequent analyses
ggplot(bh[bh$var=='bleat1' & bh$d=='d1',],
       aes(x = line,y = value,group = paste0(line,treatment),col = treatment))+
  geom_boxplot() + facet_wrap(~ sex)

m0.b2 = lme(value ~ sex*treatment + line*treatment,
            random =~ 1|ind,
            data = bh[bh$var=='bleat2' & bh$d=='d1',])
Anova(m0.b2,'III')
# Response: value
#                 Chisq Df Pr(>Chisq)  
# (Intercept)    1.1783  1    0.27771  
# sex            3.7166  1    0.05387 .
# treatment      0.0231  1    0.87916  
# line           0.0949  1    0.75800  
# sex:treatment  0.0880  1    0.76669  
# treatment:line 2.2194  1    0.13629

summary(m0.b2)
# Linear mixed-effects model fit by REML
# Data: bh[bh$var == "bleat2" & bh$d == "d1", ] 
# AIC      BIC    logLik
# 498.6301 517.4838 -241.3151
# 
# Random effects:
#   Formula: ~1 | ind
# (Intercept) Residual
# StdDev:     4.48586 1.682198
# 
# Fixed effects: value ~ sex * treatment + line * treatment 
# Value Std.Error DF    t-value p-value
# (Intercept)         2.7217391  2.507400 78  1.0854828  0.2811
# sex                 2.8565217  1.481715 78  1.9278486  0.0575
# treatmentstr       -0.5391304  3.545998 78 -0.1520391  0.8795
# lineS              -0.4565217  1.481715 78 -0.3081037  0.7588
# sex:treatmentstr   -0.6217391  2.095461 78 -0.2967075  0.7675
# treatmentstr:lineS  3.1217391  2.095461 78  1.4897624  0.1403

m0.l1 = lme(value ~ sex*treatment + line*treatment,
            random =~ 1|ind,
            data = bh[bh$var=='locomotion1' & bh$d=='d1',])
Anova(m0.l1,'III')
# Response: value
#                 Chisq Df Pr(>Chisq)   
# (Intercept)    9.3708  1   0.002205 **
# sex            1.8488  1   0.173923   
# treatment      0.8646  1   0.352467   
# line           0.2756  1   0.599629   
# sex:treatment  0.5978  1   0.439434   
# treatment:line 1.0466  1   0.306298   

summary(m0.l1)
# Linear mixed-effects model fit by REML
# Data: bh[bh$var == "locomotion1" & bh$d == "d1", ] 
# AIC     BIC    logLik
# 519.8794 538.733 -251.9397
# 
# Random effects:
#   Formula: ~1 | ind
# (Intercept) Residual
# StdDev:    5.140464 1.927674
# 
# Fixed effects: value ~ sex * treatment + line * treatment 
# Value Std.Error DF    t-value p-value
# (Intercept)         8.795652  2.873295 78  3.0611732  0.0030
# sex                 2.308696  1.697936 78  1.3597074  0.1778
# treatmentstr        3.778261  4.063452 78  0.9298155  0.3553
# lineS               0.891304  1.697936 78  0.5249341  0.6011
# sex:treatmentstr   -1.856522  2.401244 78 -0.7731501  0.4418
# treatmentstr:lineS  2.456522  2.401244 78  1.0230206  0.3095

m0.l2 = lme(value ~ sex*treatment + line*treatment,
            random =~ 1|ind,
            data = bh[bh$var=='locomotion2' & bh$d=='d1',])
Anova(m0.l2,'III')
# Response: value
#                  Chisq Df Pr(>Chisq)    
# (Intercept)    13.4991  1  0.0002387 ***
# sex             0.0626  1  0.8025036    
# treatment       0.0718  1  0.7887208    
# line            2.2616  1  0.1326184    
# sex:treatment   0.1486  1  0.6999045    
# treatment:line  0.2015  1  0.6534784

summary(m0.l2)
# Linear mixed-effects model fit by REML
# Data: bh[bh$var == "locomotion2" & bh$d == "d1", ] 
# AIC      BIC    logLik
# 484.7915 503.6452 -234.3957
# 
# Random effects:
#   Formula: ~1 | ind
# (Intercept) Residual
# StdDev:    4.105059 1.539409
# 
# Fixed effects: value ~ sex * treatment + line * treatment 
# Value Std.Error DF   t-value p-value
# (Intercept)         8.430435  2.294551 78  3.674111  0.0004
# sex                 0.339130  1.355935 78  0.250108  0.8032
# treatmentstr       -0.869565  3.244985 78 -0.267972  0.7894
# lineS              -2.039130  1.355935 78 -1.503856  0.1367
# sex:treatmentstr    0.739130  1.917581 78  0.385449  0.7010
# treatmentstr:lineS  0.860870  1.917581 78  0.448935  0.6547

m0.v1 = lme(value ~ sex*treatment + line*treatment,
            random =~ 1|ind,
            data = bh[bh$var=='vigilance1' & bh$d=='d1',])
Anova(m0.v1,'III')
# Response: value
#                  Chisq Df Pr(>Chisq)    
# (Intercept)    28.9942  1   7.26e-08 ***
# sex             0.4159  1    0.51900    
# treatment       2.9812  1    0.08424 .  
# line            0.0321  1    0.85783    
# sex:treatment   3.4671  1    0.06260 .  
# treatment:line  0.0753  1    0.78374

summary(m0.v1)
# Linear mixed-effects model fit by REML
# Data: bh[bh$var == "vigilance1" & bh$d == "d1", ] 
# AIC      BIC    logLik
# 387.7971 406.6508 -185.8986
# 
# Random effects:
#   Formula: ~1 | ind
# (Intercept)  Residual
# StdDev:    2.204416 0.8266559
# 
# Fixed effects: value ~ sex * treatment + line * treatment 
# Value Std.Error DF   t-value p-value
# (Intercept)         6.634783 1.2321720 78  5.384624  0.0000
# sex                -0.469565 0.7281358 78 -0.644887  0.5209
# treatmentstr       -3.008696 1.7425543 78 -1.726601  0.0882
# lineS              -0.130435 0.7281358 78 -0.179135  0.8583
# sex:treatmentstr    1.917391 1.0297395 78  1.862016  0.0664
# treatmentstr:lineS  0.282609 1.0297395 78  0.274447  0.7845

m0.v2 = lme(value ~ sex*treatment + line*treatment,
            random =~ 1|ind,
            data = bh[bh$var=='vigilance2' & bh$d=='d1',])
Anova(m0.v2,'III')
# Response: value
# Chisq Df Pr(>Chisq)    
# (Intercept)    13.9818  1  0.0001846 ***
# sex             4.0158  1  0.0450766 *  
# treatment       1.7924  1  0.1806322    
# line            2.0442  1  0.1527904    
# sex:treatment   1.8390  1  0.1750686    
# treatment:line  0.3660  1  0.5452022 

summary(m0.v2)
# Linear mixed-effects model fit by REML
# Data: bh[bh$var == "vigilance2" & bh$d == "d1", ] 
# AIC      BIC    logLik
# 393.9039 412.7575 -188.9519
# 
# Random effects:
#   Formula: ~1 | ind
# (Intercept)  Residual
# StdDev:     2.29242 0.8596591
# 
# Fixed effects: value ~ sex * treatment + line * treatment 
#                        Value Std.Error DF   t-value p-value
# (Intercept)         4.791304 1.2813627 78  3.739226  0.0004
# sex                 1.517391 0.7572044 78  2.003939  0.0485
# treatmentstr        2.426087 1.8121205 78  1.338811  0.1845
# lineS               1.082609 0.7572044 78  1.429744  0.1568
# sex:treatmentstr   -1.452174 1.0708488 78 -1.356096  0.1790
# treatmentstr:lineS -0.647826 1.0708488 78 -0.604965  0.5470

## Contact phase 2 - Before treatment
m0.c2 = lme(value ~ sex*treatment + line*treatment,
            random =~ 1|ind,
            data = bh[bh$var=='contact2' & bh$d=='d1',])
Anova(m0.c2,'III')
# Response: value
#                 Chisq Df Pr(>Chisq)
# (Intercept)    0.7055  1     0.4009
# sex            0.2245  1     0.6356
# treatment      0.2941  1     0.5876
# line           0.8979  1     0.3433
# sex:treatment  0.2453  1     0.6204
# treatment:line 0.2599  1     0.6102

summary(m0.c2)
# Linear mixed-effects model fit by REML
# Data: bh[bh$var == "contact2" & bh$d == "d1", ] 
# AIC      BIC    logLik
# 194.6166 213.4702 -89.30828
# 
# Random effects:
#   Formula: ~1 | ind
# (Intercept)  Residual
# StdDev:   0.6389849 0.2396193
# 
# Fixed effects: value ~ sex * treatment + line * treatment 
#                         Value Std.Error DF    t-value p-value
# (Intercept)         0.3000000 0.3571646 78  0.8399489  0.4035
# sex                 0.1000000 0.2110617 78  0.4737951  0.6370
# treatmentstr       -0.2739130 0.5051070 78 -0.5422872  0.5892
# lineS              -0.2000000 0.2110617 78 -0.9475902  0.3463
# sex:treatmentstr    0.1478261 0.2984863 78  0.4952525  0.6218
# treatmentstr:lineS  0.1521739 0.2984863 78  0.5098187  0.6116

##---==== Evaluate differences after behaviour treatment

## BLEAT phase 1 Covariate at d0 considered for bleat1
df_bleat1 = bh[bh$var=='bleat1',]
d1 = df_bleat1[df_bleat1$d=='d1',c('ind','value')]
df_bleat1 = df_bleat1[df_bleat1$d!='d1',]
df_bleat1$bl0 = d1$value[match(df_bleat1$ind,d1$ind)]

## Variable selection procedure
mb1_full = lm(value ~ treatment*(sex + d + line) + bl0,
              data = df_bleat1)
s = stepAIC(mb1_full,method='both')
s$anova
# Final Model:
#   value ~ treatment + sex + d + line + bl0 + treatment:sex

mb1 = lme(value ~ treatment + sex + d + line + bl0 + treatment:sex,
          random =~ 1|ind,
          data = df_bleat1)
anova.lme(mb1,type='marginal')
#               numDF denDF   F-value p-value
# (Intercept)       1   166  0.684358  0.4093
# treatment         1    78  0.263037  0.6095
# sex               1    78  4.331718  0.0407
# d                 2   166 15.705748  <.0001
# line              1    78  2.259032  0.1369
# bl0               1    78  6.450457  0.0131
# treatment:sex     1    78  1.487034  0.2264

summary(mb1)
# Linear mixed-effects model fit by REML
# Data: df_bleat1 
# AIC      BIC   logLik
# 1444.962 1479.934 -712.481
# 
# Random effects:
#   Formula: ~1 | ind
# (Intercept) Residual
# StdDev:     3.42353 3.304745
# 
# Fixed effects: value ~ treatment + sex + d + line + bl0 + treatment:sex 
#                       Value Std.Error  DF   t-value p-value
# (Intercept)       1.7015865 2.0568960 166  0.827259  0.4093
# treatmentstr      1.4096043 2.7484579  78  0.512871  0.6095 ns
# sex               2.5641904 1.2320269  78  2.081278  0.0407
# dd3              -2.2023810 0.5099332 166 -4.318960  0.0000
# dd4              -2.6785714 0.5099332 166 -5.252789  0.0000
# lineS             1.3295088 0.8845655  78  1.503008  0.1369
# bl0               0.2564518 0.1009742  78  2.539775  0.0131
# treatmentstr:sex -2.0929376 1.7163102  78 -1.219440  0.2264

## BLEAT phase 2
mb2_full = lm(value ~ treatment*(sex + d + line),
              data = bh[bh$var=='bleat2',])
s = stepAIC(mb2_full,method='both')
s$anova
# Final Model:
#   value ~ treatment + sex + d + line + treatment:line
rm(s,mb2_full)

mb2 = lme(value ~ treatment + sex + d + line + treatment:line,
          random =~ 1|ind,
          data = bh[bh$var=='bleat2',])
anova.lme(mb2,type='marginal')
#                numDF denDF   F-value p-value
# (Intercept)        1   249 14.530424  0.0002
# treatment          1    79  6.098883  0.0157 *
# sex                1    79  6.012504  0.0164
# d                  3   249 25.620006  <.0001
# line               1    79  0.078495  0.7801
# treatment:line     1    79  1.543366  0.2178

summary(mb2)
# Linear mixed-effects model fit by REML
# Data: bh[bh$var == "bleat2", ] 
# AIC      BIC    logLik
# 1780.99 1818.921 -880.4952
# 
# Random effects:
#   Formula: ~1 | ind
# (Intercept) Residual
# StdDev:     2.84269 2.755955
# 
# Fixed effects: value ~ treatment + sex + d + line + treatment:line 
#                        Value Std.Error  DF   t-value p-value
# (Intercept)         4.913352 1.2889578 249  3.811879  0.0002
# treatmentstr       -2.352273 0.9524946  79 -2.469592  0.0157 *
# sex                 1.694022 0.6908620  79  2.452041  0.0164 *
# dd2                -1.869048 0.4252530 249 -4.395142  0.0000
# dd3                -3.404762 0.4252530 249 -8.006438  0.0000
# dd4                -2.988095 0.4252530 249 -7.026629  0.0000
# lineS               0.273592 0.9765216  79  0.280170  0.7801
# treatmentstr:lineS  1.714773 1.3802958  79  1.242323  0.2178

ggplot(bh[bh$var=='bleat2',],
       aes(x = d,y = value,fill = treatment)) +
  #geom_smooth(method = 'lm') +
  geom_boxplot( alpha =.3) + 
  ylab('Bleat count - Phase 2') + xlab('Day') +
  facet_wrap(~ line)

## LOCOMOTION phase 1
ml1_full = lm(value ~ treatment*(sex + d + line),
              data = bh[bh$var=='locomotion1',])
s = stepAIC(ml1_full,method='both')
s$anova
# Final Model:
#   value ~ treatment + d + line + treatment:line
rm(s,ml1_full)

ml1 = lme(value ~ treatment + d + line + treatment:line,
          random =~ 1|ind,
          data = bh[bh$var=='locomotion1',])
anova.lme(ml1,type='marginal')
#                numDF denDF   F-value p-value
# (Intercept)        1   249 123.29093  <.0001
# treatment          1    80   3.97014  0.0497 
# d                  3   249  14.60307  <.0001 *
# line               1    80   0.03892  0.8441
# treatment:line     1    80   5.38681  0.0228 *

summary(ml1)
# Linear mixed-effects model fit by REML
# Data: bh[bh$var == "locomotion1", ] 
# AIC      BIC    logLik
# 2313.929 2348.094 -1147.965
# 
# Random effects:
#   Formula: ~1 | ind
# (Intercept) Residual
# StdDev:    3.917194 6.844865
# 
# Fixed effects: value ~ treatment + d + line + treatment:line 
# Value Std.Error  DF   t-value p-value
# (Intercept)        14.255141  1.283825 249 11.103645  0.0000
# treatmentstr       -3.125000  1.568365  80 -1.992521  0.0497 
# dd2                 3.607143  1.056186 249  3.415255  0.0007
# dd3                 4.952381  1.056186 249  4.688930  0.0000
# dd4                 6.738095  1.056186 249  6.379650  0.0000
# lineS              -0.317045  1.607096  80 -0.197278  0.8441
# treatmentstr:lineS  5.275000  2.272777  80  2.320949  0.0228

ggplot(bh[bh$var=='locomotion1',],
       aes(x = d,y = value, fill = treatment))+
  geom_boxplot(alpha=.3) +
  facet_wrap(~ line) + 
  ylab('Locomotion - phase 1') + xlab('Day') +
  theme(legend.position = 'bottom')

## LOCOMOTION phase 2
ml2_full = lm(value ~ treatment*(sex + d + line),
              data = bh[bh$var=='locomotion2',])
s = stepAIC(ml2_full,method='both')
s$anova
# Final Model:
#   value ~ d
rm(s,ml2_full)

ml2 = lme(value ~ d,
          random =~ 1|ind,
          data = bh[bh$var=='locomotion2',])

anova.lme(ml2,type='marginal')
#             numDF denDF   F-value p-value
# (Intercept)     1   249 162.30502  <.0001
# d               3   249  14.88553  <.0001
summary(ml2)
# Linear mixed-effects model fit by REML
# Data: bh[bh$var == "locomotion2", ] 
# AIC      BIC    logLik
# 2120.142 2142.973 -1054.071
# 
# Random effects:
#   Formula: ~1 | ind
# (Intercept) Residual
# StdDev:    3.356168 4.946879
# 
# Fixed effects: value ~ d 
# Value Std.Error  DF   t-value p-value
# (Intercept) 8.309524 0.6522441 249 12.739899  0.0000
# dd2         2.714286 0.7633200 249  3.555895  0.0005
# dd3         1.583333 0.7633200 249  2.074272  0.0391
# dd4         4.952381 0.7633200 249  6.487949  0.0000

ggplot(bh[bh$var=='locomotion2',],
       aes(x = d,y = value,fill=treatment)) +
  geom_boxplot(alpha = .3) +
  facet_wrap(~ line) + 
  ylab('Locomotion - Phase 2') + xlab('Day')

## VIGILANCE phase 1
mv1_full = lm(value ~ treatment*(sex + d + line),
              data = bh[bh$var=='vigilance1',])
s = stepAIC(mv1_full,method='both')
s$anova
# Final Model:
#   value ~ treatment + sex + d + line
rm(s,mv1_full)

mv1 = lme(value ~ treatment + sex + d + line,
          random =~ 1|ind,
          data = bh[bh$var=='vigilance1',])
anova.lme(mv1,type='marginal')
#             numDF denDF  F-value p-value
# (Intercept)     1   249 73.10180  <.0001
# treatment       1    80  1.44975  0.2321
# sex             1    80  8.53301  0.0045
# d               3   249 20.74257  <.0001
# line            1    80  1.88995  0.1730

summary(mv1)
# Linear mixed-effects model fit by REML
# Data: bh[bh$var == "vigilance1", ] 
# AIC      BIC    logLik
# 1485.728 1519.893 -733.8642
# 
# Random effects:
#   Formula: ~1 | ind
# (Intercept) Residual
# StdDev:    1.001042 1.970536
# 
# Fixed effects: value ~ treatment + sex + d + line 
# Value Std.Error  DF   t-value p-value
# (Intercept)   4.899275 0.5730174 249  8.549959  0.0000
# treatmentstr -0.369048 0.3065039  80 -1.204055  0.2321
# sex           0.897283 0.3071695  80  2.921132  0.0045
# dd2           1.047619 0.3040603 249  3.445431  0.0007
# dd3           1.773810 0.3040603 249  5.833742  0.0000
# dd4           2.250000 0.3040603 249  7.399847  0.0000
# lineS        -0.422283 0.3071695  80 -1.374754  0.1730
ggplot(bh[bh$var=='vigilance1',],
       aes(x = d,y = value, fill=treatment))+
  geom_boxplot(alpha=.3) +
  facet_wrap(~ line) +
  ylab('Vigilance - Phase 1') + xlab('Day')

## VIGILANCE phase 2 // altered
mv2_full = lm(value ~ treatment*(sex + d + line),
              data = bh[bh$var=='vigilance2',])
s = stepAIC(mv2_full,method='both')
s$anova
# Final Model:
#   value ~ treatment + sex + d + line
rm(s,mv2_full)

mv2 = lme(value ~ treatment + sex + d + line,
          random =~ 1|ind,
          data = bh[bh$var=='vigilance2',])

anova.lme(mv2,type='marginal')
#             numDF denDF   F-value p-value
# (Intercept)     1   249 115.67281  <.0001
# treatment       1    80   5.05710  0.0273 *
# sex             1    80   8.79398  0.0040
# d               3   249   5.96092  0.0006
# line            1    80   2.77401  0.0997

summary(mv2)
# Linear mixed-effects model fit by REML
# Data: bh[bh$var == "vigilance2", ] 
# AIC     BIC    logLik
# 1558.946 1593.11 -770.4729
# 
# Random effects:
#   Formula: ~1 | ind
# (Intercept) Residual
# StdDev:   0.8237961 2.281344
# 
# Fixed effects: value ~ treatment + sex + d + line 
#                  Value Std.Error  DF   t-value p-value
# (Intercept)   6.282143 0.5841067 249 10.755129  0.0000
# treatmentstr -0.690476 0.3070422  80 -2.248799  0.0273 *
# sex           0.912500 0.3077090  80  2.965465  0.0040
# dd2          -1.369048 0.3520191 249 -3.889129  0.0001
# dd3          -0.500000 0.3520191 249 -1.420378  0.1567
# dd4          -0.178571 0.3520191 249 -0.507278  0.6124
# lineS         0.512500 0.3077090  80  1.665535  0.0997

ggplot(bh[bh$var=='vigilance2',],
       aes(x = d,y = value, fill = treatment)) +
  geom_boxplot(alpha = .3) +
  facet_wrap(~ line) + 
  xlab('Day') + ylab('Vigilance - Phase2')

## CONTACT phase2 // NS
mc2_full = glm.nb(value ~ treatment*(sex + d + line),
                  data = bh[bh$var=='contact2',])
s = stepAIC(mc2_full,method='both')
s$anova
# Final Model:
#   value ~ sex + d + line

rm(s,mc2_full)

require(lmerTest)
# mc2 = glmer.nb(value ~ sex + d + line + (1|ind),
#           data = bh[bh$var=='contact2',],
#           glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000)))
# #NB does not converge

ggplot(bh[bh$var=='contact2',],
       aes(x = d,y = value,fill=treatment)) +
  geom_boxplot() +
  facet_wrap(~ line)

#Consider contact or not
bh_cont = bh[bh$var=='contact2',]
bh_cont$cont = 0
bh_cont$cont[bh_cont$value>0]=1
bh_cont$cont = factor(bh_cont$cont)

table(bh_cont$cont,bh_cont$sex)
#     1   2
# 0 127 135
# 1  33  41

table(bh_cont$cont,bh_cont$d)
#   d1 d2 d3 d4
# 0 61 66 72 63
# 1 23 18 12 21

ml3 = glmer(cont ~ sex + d + line + (1|ind),
            family = binomial,
            data = bh_cont)

car::Anova(ml3,'III')
# Response: cont
#              Chisq Df Pr(>Chisq)
# (Intercept) 2.4631  1     0.1166
# sex         0.0858  1     0.7696
# d           6.0682  3     0.1083
# line        2.1307  1     0.1444 

summary(ml3)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: binomial  ( logit )
# Formula: cont ~ sex + d + line + (1 | ind)
# Data: bh_cont
# 
# AIC      BIC   logLik deviance df.resid 
# 337.3    364.1   -161.7    323.3      329 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -1.3934 -0.3942 -0.3058 -0.1772  2.2515 
# 
# Random effects:
#   Groups Name        Variance Std.Dev.
# ind    (Intercept) 2.176    1.475   
# Number of obs: 336, groups:  ind, 84
# 
# Fixed effects:
#             Estimate Std. Error z value Pr(>|z|)  
# (Intercept)  -1.2839     0.8181  -1.569   0.1166  
# sex           0.1362     0.4652   0.293   0.7696  
# dd2          -0.4431     0.4239  -1.045   0.2959  
# dd3          -1.0918     0.4650  -2.348   0.0189 *
# dd4          -0.1701     0.4128  -0.412   0.6803  
# lineS        -0.6873     0.4708  -1.460   0.1444 

##------====== G1 stress env
g1$sol0 = scale(g1$sol0) ## scale g1 ebv
xp15 = melt(g1,c('IPG','ipg_pe','lot','sexe','gen','cpt','sol0'))
xp15$ind = xp15$IPG-20000e6
xp15$sire = xp15$ipg_pe-20000e6
xp15$ipg_pe = NULL
#xp15$ebv = phenodiv$ebv_fec[match(xp15$IPG,phenodiv$IPG)]
xp15 = xp15[,c('ind','sire','sexe','gen','cpt','sol0','variable','value')]
xp15 = xp15[grep('opg',xp15$variable),]

##------====== G2 - T.colubriformis env
rm(wk2)
g2$sol0 = scale(g2$sol0) ## scale g2 ebv
wk2 = melt(g2,c('IPG','pere','lot','sexe','gen','inoc','sol0'))
xp17 = wk2[grep('opg',wk2$variable),]
xp17$ind = xp17$IPG-20000e6
xp17$lot = NULL
#xp17$ebv = phenodiv$ebv_fec[match(xp17$IPG,phenodiv$IPG)]
xp17 = xp17[,c('ind','pere','sexe','gen','inoc','sol0','variable','value')]
colnames(xp17) = colnames(xp15)
colnames(xp17)[5] = 'env'
colnames(xp15)[5] = 'env'
xp15$xp="Chronic Stress effect"
xp17$xp="Intestinal Species effect"

##-- Group together 2015 & 2017
gbe = rbind(xp15,xp17)
gbe$variable = factor(gbe$variable)
colnames(gbe)[6]='ebv'

## Raw Tcol FEC by line and day
aggregate(value ~ variable + gen + env, data=gbe, FUN=mean)
#    variable  gen env      value
# 1    opg241 R-G1 nos  454.78261
# 2    opg301 R-G1 nos 1497.82609
# 3    opg242 R-G1 nos  174.56522
# 4    opg302 R-G1 nos 1016.30435
# 5    opg241 S-G1 nos 3982.50000
# 6    opg301 S-G1 nos 6697.50000
# 7    opg242 S-G1 nos  671.50000
# 8    opg302 S-G1 nos 6069.00000
# 9    opg241 R-G1 str  628.91304
# 10   opg301 R-G1 str  717.39130
# 11   opg242 R-G1 str  269.56522
# 12   opg302 R-G1 str  202.82609
# 13   opg241 S-G1 str 3009.50000
# 14   opg301 S-G1 str 8517.50000
# 15   opg242 S-G1 str  420.00000
# 16   opg302 S-G1 str 3761.50000
# 17   opg241 R-G2  hc   36.11111
# 18   opg301 R-G2  hc  941.66667
# 19   opg302 R-G2  hc 1572.22222
# 20   opg241 S-G2  hc 1030.00000
# 21   opg301 S-G2  hc 3132.50000
# 22   opg302 S-G2  hc 5917.50000
# 23   opg241 R-G2  tc  302.50000
# 24   opg301 R-G2  tc  390.00000
# 25   opg302 R-G2  tc  400.00000
# 26   opg241 S-G2  tc  442.10526
# 27   opg301 S-G2  tc  352.63158
# 28   opg302 S-G2  tc  589.47368

##-----===== Check wether R and S lambs are different towards T. colubriformis =====-----####
tc = gbe[gbe$env=='tc',]

m = lme(value ~ sexe + gen*variable,
            random =~ 1|ind,
            data = tc)
summary(m)
# Linear mixed-effects model fit by REML
# Data: tc 
# AIC      BIC    logLik
# 1528.038 1552.343 -755.0192
# 
# Random effects:
#   Formula: ~1 | ind
# (Intercept) Residual
# StdDev:    136.5449  178.068
# 
# Fixed effects: value ~ sexe + gen * variable 
#                             Value Std.Error DF   t-value p-value
# (Intercept)             276.26126  57.16401 74  4.832783  0.0000
# sexe2                    52.47748  54.77544 36  0.958048  0.3444
# genS-G2                 138.22428  71.90173 36  1.922406  0.0625
# variableopg301           87.50000  56.31005 74  1.553897  0.1245
# variableopg302           97.50000  56.31005 74  1.731485  0.0875
# genS-G2:variableopg301 -176.97368  80.67545 74 -2.193650  0.0314
# genS-G2:variableopg302   49.86842  80.67545 74  0.618136  0.5384
# Correlation: 
#   (Intr) sexe2  gnS-G2 vrb301 vrb302 gS-G2:301
# sexe2                  -0.479                                      
# genS-G2                -0.603 -0.020                               
# variableopg301         -0.493  0.000  0.392                        
# variableopg302         -0.493  0.000  0.392  0.500                 
# genS-G2:variableopg301  0.344  0.000 -0.561 -0.698 -0.349          
# genS-G2:variableopg302  0.344  0.000 -0.561 -0.349 -0.698  0.500   
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -1.99976997 -0.56264236 -0.03683956  0.49159413  2.82935313 
# 
# Number of Observations: 117
# Number of Groups: 39 

## Normalized fec data
gbe$value = sqrt(sqrt(gbe$value))

## Data structure
table(gbe$variable,gbe$env)
#        nos str hc tc
# opg241  43  43 38 39
# opg301  43  43 38 39
# opg242  43  43  0  0
# opg302  43  43 38 39

##-----===== Normalize FEC data within experiment and time =====-----####
gbe = plyr::ddply(gbe, c("env", "variable"), transform, fecs = scale(value))

table(gbe$sire,gbe$env)
#        nos str hc tc
# 132336  60  56  0  0
# 132361  28  36  0  0
# 132453  16  24  0  0
# 132471  16  12  0  0
# 132497  40  36  0  0
# 132550  12   8  0  0
# 152226   0   0 21 21
# 152258   0   0 21 21
# 152337   0   0 30 30
# 152347   0   0 21 21
# 152365   0   0 12 18
# 152434   0   0  9  6

gbe$sire = factor(gbe$sire)
gbe$ind = factor(gbe$ind)
gbe$ram = match(gbe$sire,levels(gbe$sire))
gbe$id = match(gbe$ind,levels(gbe$ind))

###-------- Plot by sire
gbe$dinf = as.integer(sapply(str_split(as.character(gbe$variable),'opg'),
                             function(x) paste0(substr(x[2],3,3),substr(x[2],1,2))))

gbe$dinf = factor(gbe$dinf)
gbe$xp = factor(gbe$xp)
gbe$sexe = factor(gbe$sexe)
gbe$ind = factor(gbe$ind)
gbe$sire = factor(gbe$sire)
gbe$d = factor(match(gbe$dinf,levels(gbe$dinf)))

gbe$time = as.integer(substr(gbe$dinf,2,3))
gbe$inf = as.integer(substr(gbe$dinf,1,1))
gbe$time = paste0('D',gbe$time)
gbe$time2 = factor(paste0('Infection ',gbe$inf,' - ',gbe$time))

###===----- G x E plot
gbe$Environment = 'No stress'
gbe$Environment[gbe$env=='str'] = 'Chronic stress'
gbe$Environment[gbe$env=='hc'] = 'H.contortus'
gbe$Environment[gbe$env=='tc'] = 'T.colubriformis'
gbe$Environment = factor(gbe$Environment,levels=c('No stress','Chronic stress',
                                                  'H.contortus','T.colubriformis'))
gbe$Line = factor(substr(gbe$gen,1,1))

## Genotype expression across environment
ggplot(gbe,aes(x = ebv, y = fecs,col = Environment,fill = Environment)) +
  geom_smooth(method = 'lm',se = TRUE,alpha=.2,lwd=.9) +
  geom_point(alpha=.3) + facet_wrap(~ xp,ncol = 1) +
  scale_color_manual(values=colEnv) +
  scale_fill_manual(values=colEnv) +
  ylab('Scaled FEC') + xlab('Scaled genetic value') +
  theme(legend.position = 'bottom',text = element_text(size=14))

## Through time
# p1 = ggplot(gbe[gbe$xp=='Chronic Stress effect',],
#             aes(x = ebv,y = fecs,col = Environment,fill = Environment)) +
#   geom_smooth(method = 'lm',se = TRUE,alpha=.2,lwd=.9) +
#   geom_point(alpha=.3,
#              aes(x = ebv,y = fecs,col = Environment,fill = Environment,shape = Line)) + 
#   facet_wrap(~ time2,ncol = 4) +
#   xlab(bquote('Estimated Breeding Value ('*sigma[g]*')')) +
#   ylab(bquote('Faecal Egg Count ('*sigma[p]*')')) +
#   scale_color_manual(values = colEnv) +
#   scale_fill_manual(values = colEnv) +
#   theme(legend.position = 'bottom',text = element_text(size=14)) +
#   guides(shape = guide_legend(override.aes = list(size = 3))) +
#   scale_y_continuous(breaks = seq(-2,3, 1), limits = c(-2,3))

# #####
# p1 = ggplot(gbe[gbe$xp=='Chronic Stress effect',],
#             aes(x = Line,y = fecs,col = Environment,fill = Environment)) +
#   #geom_smooth(method = 'lm',se = TRUE,alpha=.2,lwd=.9) +
#   geom_boxplot(alpha = .4) +
#   facet_wrap(~ time2,ncol = 4) +
#   xlab(bquote('Estimated Breeding Value ('*sigma[g]*')')) +
#   ylab(bquote('Faecal Egg Count ('*sigma[p]*')')) +
#   scale_color_manual(values = colEnv) +
#   scale_fill_manual(values = colEnv) +
#   theme(legend.position = 'bottom',text = element_text(size=14)) +
#   guides(shape = guide_legend(override.aes = list(size = 3))) +
#      scale_y_continuous(breaks = seq(-3,3, 1), limits = c(-3,3))
# p2 = ggplot(gbe[gbe$xp=='Intestinal Species effect',],
#             aes(x = Line,y = fecs,col = Environment,fill = Environment)) +
#   #geom_smooth(method = 'lm',se = TRUE,alpha=.2,lwd=.9) +
#   geom_boxplot(alpha = .4) +
#   facet_wrap(~ time2,ncol = 4) +
#   xlab(bquote('Estimated Breeding Value ('*sigma[g]*')')) +
#   ylab(bquote('Faecal Egg Count ('*sigma[p]*')')) +
#   scale_color_manual(values = colEnv[3:4]) +
#   scale_fill_manual(values = colEnv[3:4]) +
#   theme(legend.position = 'bottom',text = element_text(size=14)) +
#   guides(shape = guide_legend(override.aes = list(size = 3))) +
#   scale_y_continuous(breaks = seq(-3,3, 1), limits = c(-3,3))
# !!!!!#####

p2 = ggplot(gbe[gbe$xp=='Intestinal Species effect',],
            aes(x = ebv,y = fecs, col = Environment,fill = Environment)) +
  geom_smooth(method = 'lm',se = TRUE,alpha=.2,lwd=.9) +
  geom_point(alpha=.3,
             aes(x = ebv,y = fecs,col = Environment,fill = Environment,shape = Line)) + 
  facet_wrap(~ time2,ncol = 4) +
  xlab(bquote('Estimated Breeding Value ('*sigma[g]*')')) +
  ylab(bquote('Faecal Egg Count ('*sigma[p]*')')) +
  scale_y_continuous(breaks = seq(-2,3, 1), limits = c(-2,3)) +
  theme(legend.position = 'bottom',text = element_text(size=14)) +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(values = colEnv[3:4],
                     labels = c(expression(italic('H. contortus')),
                                expression(italic('T. colubriformis')))) +
  scale_fill_manual(values = colEnv[3:4],
                    labels = c(expression(italic('H. contortus')),
                               expression(italic('T. colubriformis')))) 

pdf(file=paste0(draft,'./Figure2.pdf')) #,width=10,height=12)
multiplot(p1,p2)
dev.off()

###-------- What interactions ? When ?
gbe$d = as.integer(as.character(substr(gbe$time,2,3)))
gbe$d[gbe$inf == 2] = gbe$d[gbe$inf == 2] + 45
gbe$d = factor(gbe$d)
gbe$g = factor(substr(gbe$gen,1,1)) ## genetic group

## Genotype x stress 
mgst = lme(fecs ~ d*g + g*env, 
           random =~ 1|ind,
          data = gbe[gbe$xp=='Chronic Stress effect',])
anova(mgst)
#             numDF denDF  F-value p-value
# (Intercept)     1   252  0.00000  1.0000
# d               3   252  0.00000  1.0000
# g               1    82 75.91776  <.0001
# env             1    82  0.00000  1.0000
# d:g             3   252  4.62594  0.0036
# g:env           1    82  0.02803  0.8675

summary(mgst)
# Linear mixed-effects model fit by REML
# Data: gbe[gbe$xp == "Chronic Stress effect", ] 
# AIC      BIC    logLik
# 875.2522 920.9859 -425.6261
# 
# Random effects:
#   Formula: ~1 | ind
# (Intercept)  Residual
# StdDev:   0.4001002 0.7444159
# 
# Fixed effects: fecs ~ d * g + g * env 
#                  Value Std.Error  DF   t-value p-value
# (Intercept) -0.5610017 0.1483866 252 -3.780678  0.0002
# d30         -0.0522059 0.1552214 252 -0.336332  0.7369
# d69          0.3083369 0.1552214 252  1.986433  0.0481
# d75          0.0359739 0.1552214 252  0.231759  0.8169
# gS           1.2061536 0.2175774  82  5.543562  0.0000
# envstr       0.0183993 0.1611424  82  0.114180  0.9094
# d30:gS       0.1122427 0.2275993 252  0.493159  0.6223
# d69:gS      -0.6629244 0.2275993 252 -2.912682  0.0039
# d75:gS      -0.0773439 0.2275993 252 -0.339825  0.7343
# gS:envstr   -0.0395584 0.2362811  82 -0.167421  0.8675
# Correlation: 
#   (Intr) d30    d69    d75    gS     envstr d30:gS d69:gS d75:gS
# d30       -0.523                                                        
# d69       -0.523  0.500                                                 
# d75       -0.523  0.500  0.500                                          
# gS        -0.682  0.357  0.357  0.357                                   
# envstr    -0.543  0.000  0.000  0.000  0.370                            
# d30:gS     0.357 -0.682 -0.341 -0.341 -0.523  0.000                     
# d69:gS     0.357 -0.341 -0.682 -0.341 -0.523  0.000  0.500              
# d75:gS     0.357 -0.341 -0.341 -0.682 -0.523  0.000  0.500  0.500       
# gS:envstr  0.370  0.000  0.000  0.000 -0.543 -0.682  0.000  0.000  0.000
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -2.5440700 -0.6844306  0.0516469  0.6249162  2.7221794 
# 
# Number of Observations: 344
# Number of Groups: 86 

##-- Genotype x Species : reduced divergence
mgsp = lme(fecs ~ g*d + g*env,
           random =~ 1|ind,
         data = gbe[gbe$xp=='Intestinal Species effect',])

anova.lme(mgsp)
#             numDF denDF  F-value p-value
# (Intercept)     1   150  0.00000  1.0000
# g               1    73 34.93955  <.0001
# d               2   150  0.00000  1.0000
# env             1    73  0.05360  0.8176
# g:d             2   150  3.89051  0.0225
# g:env           1    73  9.76833  0.0025

summary(mgsp)
# Linear mixed-effects model fit by REML
# Data: gbe[gbe$xp == "Intestinal Species effect", ] 
# AIC      BIC   logLik
# 599.2139 633.2857 -289.607
# 
# Random effects:
#   Formula: ~1 | ind
# (Intercept)  Residual
# StdDev:   0.4742686 0.7267584
# 
# Fixed effects: fecs ~ g * d + g * env 
#                  Value Std.Error  DF   t-value p-value
# (Intercept) -0.7968186 0.1776049 150 -4.486467  0.0000
# gS           1.5218450 0.2462147  73  6.180967  0.0000
# d30          0.2967981 0.1667298 150  1.780114  0.0771
# d75          0.0215120 0.1667298 150  0.129023  0.8975
# envtc        0.4913003 0.2057346  73  2.388030  0.0195
# gS:d30      -0.5859860 0.2342752 150 -2.501272  0.0134
# gS:d75      -0.0424725 0.2342752 150 -0.181293  0.8564
# gS:envtc    -0.9030335 0.2889307  73 -3.125433  0.0025
# Correlation: 
#   (Intr) gS     d30    d75    envtc  gS:d30 gS:d75
# gS       -0.721                                          
# d30      -0.469  0.339                                   
# d75      -0.469  0.339  0.500                            
# envtc    -0.610  0.440  0.000  0.000                     
# gS:d30    0.334 -0.476 -0.712 -0.356  0.000              
# gS:d75    0.334 -0.476 -0.356 -0.712  0.000  0.500       
# gS:envtc  0.434 -0.595  0.000  0.000 -0.712  0.000  0.000
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -3.5983230 -0.4837153  0.1316387  0.5951407  1.9675975 
# 
# Number of Observations: 231
# Number of Groups: 77 

###-------- What genotype is more affected ? 
## By sire: some are more affected than others
msiresp = lme(fecs ~ d*sire + env*sire,
           random =~ 1|ind,
           data = gbe[gbe$xp=='Intestinal Species effect',])
summary(msiresp) ## no significant interaction between sire and environment
# Linear mixed-effects model fit by REML
# Data: gbe[gbe$xp == "Intestinal Species effect", ] 
# AIC      BIC    logLik
# 614.3631 701.0137 -281.1815
# 
# Random effects:
#   Formula: ~1 | ind
# (Intercept)  Residual
# StdDev:   0.4146541 0.7377842
# 
# Fixed effects: fecs ~ d * sire + env * sire 
#                       Value Std.Error  DF   t-value p-value
# (Intercept)      -0.9207072 0.2764111 142 -3.330934  0.0011
# d30               0.4483058 0.2788562 142  1.607659  0.1101
# d75               0.2766185 0.2788562 142  0.991975  0.3229
# sire152258        1.3724637 0.3909044  65  3.510996  0.0008
# sire152337        1.9533598 0.3603960  65  5.420037  0.0000
# sire152347        0.3160139 0.3909044  65  0.808418  0.4218
# sire152365       -0.0330896 0.4483709  65 -0.073800  0.9414
# sire152434        1.2309588 0.5165009  65  2.383266  0.0201
# envtc             0.0149356 0.3177510  65  0.047004  0.9627
# d30:sire152258   -0.7352285 0.3943622 142 -1.864348  0.0643
# d75:sire152258   -0.0895489 0.3943622 142 -0.227073  0.8207
# d30:sire152337   -0.8266205 0.3635840 142 -2.273534  0.0245
# d75:sire152337   -0.4770881 0.3635840 142 -1.312181  0.1916
# d30:sire152347   -0.1395171 0.3943622 142 -0.353779  0.7240
# d75:sire152347   -0.3517334 0.3943622 142 -0.891904  0.3740
# d30:sire152365   -0.3804052 0.4320022 142 -0.880563  0.3800
# d75:sire152365   -0.4769778 0.4320022 142 -1.104110  0.2714
# d30:sire152434   -0.3873286 0.5435908 142 -0.712537  0.4773
# d75:sire152434   -0.1620261 0.5435908 142 -0.298066  0.7661
# sire152258:envtc -0.5386071 0.4493677  65 -1.198589  0.2350
# sire152337:envtc -0.5316172 0.4142966  65 -1.283180  0.2040
# sire152347:envtc  0.8262358 0.4493677  65  1.838663  0.0705
# sire152365:envtc  0.7263519 0.4982041  65  1.457940  0.1497
# sire152434:envtc  0.3634222 0.6288473  65  0.577918  0.5653

msirestress = lme(fecs ~ d*sire + env*sire,
              random =~ 1|ind,
              data = gbe[gbe$xp=='Chronic Stress effect',])
summary(msirestress) ## Significant interaction for sire 132361, day 69 and d75
# Linear mixed-effects model fit by REML
# Data: gbe[gbe$xp == "Chronic Stress effect", ] 
# AIC      BIC    logLik
# 889.505 1009.486 -412.7525
# 
# Random effects:
#   Formula: ~1 | ind
# (Intercept)  Residual
# StdDev:   0.3434933 0.7438707
# 
# Fixed effects: fecs ~ d * sire + env * sire 
#                        Value Std.Error  DF   t-value p-value
# (Intercept)       -0.5050819 0.1771972 240 -2.850394  0.0047
# d30               -0.1065541 0.1953500 240 -0.545452  0.5859
# d69                0.1765587 0.1953500 240  0.903807  0.3670
# d75                0.0308272 0.1953500 240  0.157805  0.8747
# sire132361         0.8065036 0.3065200  74  2.631162  0.0103
# sire132453         0.1125256 0.3701088  74  0.304034  0.7620
# sire132471        -0.4781614 0.3934038  74 -1.215447  0.2281
# sire132497         1.3785400 0.2808447  74  4.908550  0.0000
# sire132550         1.0315288 0.4470364  74  2.307482  0.0238
# envstr            -0.0652795 0.1881411  74 -0.346971  0.7296
# d30:sire132361     0.0080861 0.3276119 240  0.024682  0.9803
# d69:sire132361    -0.7949202 0.3276119 240 -2.426408  0.0160
# d75:sire132361    -0.5586632 0.3276119 240 -1.705259  0.0894
# d30:sire132453    -0.1735745 0.3857854 240 -0.449925  0.6532
# d69:sire132453     0.1161487 0.3857854 240  0.301071  0.7636
# d75:sire132453     0.0124633 0.3857854 240  0.032306  0.9743
# d30:sire132471     0.6051085 0.4430122 240  1.365896  0.1733
# d69:sire132471     0.7000443 0.4430122 240  1.580192  0.1154
# d75:sire132471     0.0160163 0.4430122 240  0.036153  0.9712
# d30:sire132497     0.1908450 0.3104969 240  0.614644  0.5394
# d69:sire132497    -0.4555986 0.3104969 240 -1.467321  0.1436
# d75:sire132497     0.1466293 0.3104969 240  0.472241  0.6372
# d30:sire132550     0.5816405 0.5094105 240  1.141791  0.2547
# d69:sire132550     0.0258497 0.5094105 240  0.050744  0.9596
# d75:sire132550     0.6529531 0.5094105 240  1.281782  0.2012
# sire132361:envstr  0.7254361 0.3170096  74  2.288372  0.0250 *
# sire132453:envstr -0.0161830 0.3770922  74 -0.042915  0.9659
# sire132471:envstr  0.5570900 0.4300223  74  1.295491  0.1992
# sire132497:envstr -0.1941468 0.2991819  74 -0.648926  0.5184
# sire132550:envstr -1.0719037 0.4989990  74 -2.148108  0.0350 *

p3 = ggplot(gbe[gbe$xp=='Chronic Stress effect',],
       aes(x = paste0(g,'-',sire), y = fecs ,col = Environment, fill = Environment, group = paste0(sire,env))) +
  geom_boxplot(alpha = .3) + #facet_wrap(~ time,ncol = 4) +
  scale_color_manual(values = colEnv) +
  scale_fill_manual(values = colEnv) +
  ylab(bquote('Faecal Egg Count ('*sigma[p]*')')) +
  xlab('Sire') +
  theme(legend.position = 'bottom',text = element_text(size=14),axis.text.x = element_text(angle=45, vjust=0.5))

p4 = ggplot(gbe[gbe$xp=='Intestinal Species effect',],
            aes(x = paste0(g,'-',sire), y = fecs ,col = Environment, fill = Environment, group = paste0(sire,env))) +
  geom_boxplot(alpha = .3) + #facet_wrap(~ time,ncol = 4) +
  scale_color_manual(values = colEnv[3:4],labels = c(expression(italic('H. contortus')),
                                                     expression(italic('T. colubriformis')))) +
  scale_fill_manual(values = colEnv[3:4],labels = c(expression(italic('H. contortus')),
                                                   expression(italic('T. colubriformis')))) +
  ylab(bquote('Faecal Egg Count ('*sigma[p]*')')) +
  xlab('Sire') +
  theme(legend.position = 'bottom',text = element_text(size=14),axis.text.x = element_text(angle=45, vjust=0.5))

pdf(file=paste0(draft,'supplementary_Figure3.pdf'))
multiplot(p3,p4)
dev.off()

###-------- Any effect of the environment experienced by the dams ?
g2$dam = phenodiv$ipg_me[match(g2$IPG,phenodiv$IPG)]
transgen = na.omit(g2)

transgen$gp_dam = g1$cpt[match(transgen$dam,g1$IPG)]
table(transgen$gp_dam)
# nos str 
# 8  10 

transgen = reshape2::melt(transgen[,c('IPG','sexe','lot','gp_dam','opg241', 'opg301', 'opg302')],1:4)
table(transgen$gp_dam,transgen$lot)
#     R-hc S-hc
# nos    9   15
# str   12   18

ggplot(transgen,aes(x=variable,y=value+1,col=gp_dam)) +
  geom_boxplot() +
  scale_y_log10()

plot(density(transgen$value))

m = lme(sqrt(sqrt(value)) ~ sexe + lot + gp_dam + variable,
    random =~ 1|IPG,
    data = na.omit(transgen))

########---------------------------===== SUMMARY STATS for every considered traits =====--------------------------------------########

## FEC data
a = aggregate(fropg1 ~ Gene + line, data = phenodiv, FUN = mean)
a[,4] = aggregate(fropg1 ~ Gene + line, data = phenodiv, FUN = sd)[,3]
a[,5] = aggregate(fropg1 ~ Gene + line, data = phenodiv, FUN = min)[,3]
a[,6] = aggregate(fropg1 ~ Gene + line, data = phenodiv, FUN = max)[,3]
colnames(a) = c('Generation', 'Line','mean','sd','min','max')
a
# Generation Line      mean       sd      min       max
# 1         G0    R  7.039057 2.587705 1.000000 11.869860
# 2         G1    R  3.651333 2.277235 0.000000  8.130171
# 3         G2    R  3.253667 2.408326 0.000000  7.579289
# 4         G0    S 10.140820 1.927481 4.731071 13.989511
# 5         G1    S  7.741319 1.967647 2.524634 10.783829
# 6         G2    S  6.316833 1.527282 2.236068  8.882492

a = aggregate(fropg2 ~ Gene + line, data = phenodiv, FUN = mean)
a[,4] = aggregate(fropg2 ~ Gene + line, data = phenodiv, FUN = sd)[,3]
a[,5] = aggregate(fropg2 ~ Gene + line, data = phenodiv, FUN = min)[,3]
a[,6] = aggregate(fropg2 ~ Gene + line, data = phenodiv, FUN = max)[,3]
colnames(a) = c('Generation', 'Line','mean','sd','min','max')
a
#   Generation Line      mean       sd      min       max
# 1         G0    R  7.382098 2.497168 1.000000 11.464650
# 2         G1    R  2.006370 2.017517 0.000000  7.161833
# 3         G2    R  3.883798 3.593453 0.000000  8.745587
# 4         G0    S 10.546093 1.573784 4.898979 13.444889
# 5         G1    S  4.775108 2.730001 0.000000  9.048738
# 6         G2    S  8.584344 1.074861 6.373397 10.925812

phenodiv$fropgt = sqrt(sqrt((phenodiv$fropg1^4+phenodiv$fropg2^4)/2))
a = aggregate(fropgt ~ Gene + line, data = phenodiv, FUN = mean)
a[,4] = aggregate(fropgt ~ Gene + line, data = phenodiv, FUN = sd)[,3]
a[,5] = aggregate(fropgt ~ Gene + line, data = phenodiv, FUN = min)[,3]
a[,6] = aggregate(fropgt ~ Gene + line, data = phenodiv, FUN = max)[,3]
colnames(a) = c('Generation', 'Line','mean','sd','min','max')
a
#   Generation Line      mean        sd      min       max
# 1         G0    R  7.737187 1.7785072 2.672345 10.336053
# 2         G1    R  3.599005 1.7555414 0.000000  6.898579
# 3         G2    R  4.484873 2.3819911 0.000000  7.902530
# 4         G0    S 10.528183 1.2902965 7.362604 13.467974
# 5         G1    S  7.027353 1.6523722 3.756271  9.536822
# 6         G2    S  7.785662 0.9756063 5.533410  9.382648

a = aggregate(fecs ~ Gene + line, data = phenodiv, FUN = mean)
a[,4] = aggregate(fecs ~ Gene + line, data = phenodiv, FUN = sd)[,3]
a[,5] = aggregate(fecs ~ Gene + line, data = phenodiv, FUN = min)[,3]
a[,6] = aggregate(fecs ~ Gene + line, data = phenodiv, FUN = max)[,3]
colnames(a) = c('Generation', 'Line','mean','sd','min','max')
a
#   Generation Line       mean        sd        min       max
# 1         G0    R -0.7380955 0.8109464 -3.0434185 0.6675773
# 2         G1    R -0.6573451 0.6532137 -1.8079483 0.9808895
# 3         G2    R -0.7650879 0.9033903 -2.1014971 0.8511910
# 4         G0    S  0.6709959 0.6050677 -0.6732451 2.0608602
# 5         G1    S  0.7375091 0.7845638 -0.8765047 1.9026315
# 6         G2    S  0.6885791 0.4043688 -0.1829788 1.3647498

##---- Midparent / offspring scaled values
## G1
a = aggregate(fs1 ~ Gene + line, data = gen1[gen1$pheno=='OUI',], FUN = mean)
a[,4] = aggregate(fs1 ~ Gene + line, data = gen1[gen1$pheno=='OUI',], FUN = sd)[,3]
a[,5] = aggregate(fs1 ~ Gene + line, data = gen1[gen1$pheno=='OUI',], FUN = min)[,3]
a[,6] = aggregate(fs1 ~ Gene + line, data = gen1[gen1$pheno=='OUI',], FUN = max)[,3]
colnames(a) = c('Generation', 'Line','mean','sd','min','max')
a
# Generation Line       mean        sd       min       max
# 1       _G1_    R -0.6523139 0.7706872 -1.888038 0.8634641
# 2       _G1_    S  0.7318644 0.6659129 -1.033624 1.7615444

a = aggregate(fs2 ~ Gene + line, data = gen1[gen1$pheno=='OUI',], FUN = mean)
a[,4] = aggregate(fs2 ~ Gene + line, data = gen1[gen1$pheno=='OUI',], FUN = sd)[,3]
a[,5] = aggregate(fs2 ~ Gene + line, data = gen1[gen1$pheno=='OUI',], FUN = min)[,3]
a[,6] = aggregate(fs2 ~ Gene + line, data = gen1[gen1$pheno=='OUI',], FUN = max)[,3]
colnames(a) = c('Generation', 'Line','mean','sd','min','max')
a
# Generation Line       mean        sd       min       max
# 1       _G1_    R -0.6523139 0.7706872 -1.888038 0.8634641
# 2       _G1_    S  0.7318644 0.6659129 -1.033624 1.7615444

a = aggregate(midparent_fopg1 ~ Gene + line, data = gen1[gen1$pheno=='OUI',], FUN = mean)
a[,4] = aggregate(midparent_fopg1 ~ Gene + line, data = gen1[gen1$pheno=='OUI',], FUN = sd)[,3]
a[,5] = aggregate(midparent_fopg1 ~ Gene + line, data = gen1[gen1$pheno=='OUI',], FUN = min)[,3]
a[,6] = aggregate(midparent_fopg1 ~ Gene + line, data = gen1[gen1$pheno=='OUI',], FUN = max)[,3]
colnames(a) = c('Generation', 'Line','mean','sd','min','max')
a
#   Generation Line       mean        sd        min        max
# 1       _G1_    R -1.7086250 0.6765939 -2.9441503 -0.1894161
# 2       _G1_    S  0.7512839 0.4337007 -0.0130751  1.5694649

a = aggregate(midparent_fopg2 ~ Gene + line, data = gen1[gen1$pheno=='OUI',], FUN = mean)
a[,4] = aggregate(midparent_fopg2 ~ Gene + line, data = gen1[gen1$pheno=='OUI',], FUN = sd)[,3]
a[,5] = aggregate(midparent_fopg2 ~ Gene + line, data = gen1[gen1$pheno=='OUI',], FUN = min)[,3]
a[,6] = aggregate(midparent_fopg2 ~ Gene + line, data = gen1[gen1$pheno=='OUI',], FUN = max)[,3]
colnames(a) = c('Generation', 'Line','mean','sd','min','max')
a
#   Generation Line       mean        sd       min       max
# 1       _G1_    R -0.4480384 0.5472210 -1.721565 0.3490481
# 2       _G1_    S  0.5036376 0.3483212 -0.454955 1.1594208

a = aggregate(midparent_fopg ~ Gene + line, data = gen1[gen1$pheno=='OUI',], FUN = mean)
a[,4] = aggregate(midparent_fopg ~ Gene + line, data = gen1[gen1$pheno=='OUI',], FUN = sd)[,3]
a[,5] = aggregate(midparent_fopg ~ Gene + line, data = gen1[gen1$pheno=='OUI',], FUN = min)[,3]
a[,6] = aggregate(midparent_fopg ~ Gene + line, data = gen1[gen1$pheno=='OUI',], FUN = max)[,3]
colnames(a) = c('Generation', 'Line','mean','sd','min','max')
a
#   Generation Line       mean        sd       min       max
# 1       _G1_    R -1.2851357 0.4884900 -2.7283425 -0.4317925
# 2       _G1_    S  0.7515785 0.3897972 -0.1047734  1.4491418

## G2
a = aggregate(fs1 ~ Gene + line, data = gen2mid[gen2mid$pheno=='OUI',], FUN = mean)
a[,4] = aggregate(fs1 ~ Gene + line, data = gen2mid[gen2mid$pheno=='OUI',], FUN = sd)[,3]
a[,5] = aggregate(fs1 ~ Gene + line, data = gen2mid[gen2mid$pheno=='OUI',], FUN = min)[,3]
a[,6] = aggregate(fs1 ~ Gene + line, data = gen2mid[gen2mid$pheno=='OUI',], FUN = max)[,3]
colnames(a) = c('Generation', 'Line','mean','sd','min','max')
a
# Generation Line       mean        sd       min      max
# 1       G2_1    R -0.7211062 0.6647347 -1.488620 1.096188
# 2       G2_1    S  0.6850509 0.7532566 -1.118041 1.889172

a = aggregate(fs2 ~ Gene + line, data = gen2mid[gen2mid$pheno=='OUI',], FUN = mean)
a[,4] = aggregate(fs2 ~ Gene + line, data = gen2mid[gen2mid$pheno=='OUI',], FUN = sd)[,3]
a[,5] = aggregate(fs2 ~ Gene + line, data = gen2mid[gen2mid$pheno=='OUI',], FUN = min)[,3]
a[,6] = aggregate(fs2 ~ Gene + line, data = gen2mid[gen2mid$pheno=='OUI',], FUN = max)[,3]
colnames(a) = c('Generation', 'Line','mean','sd','min','max')
a
# Generation Line       mean        sd         min       max
# 1       G2_1    R -0.7021285 0.9854043 -1.76521403 0.7187086
# 2       G2_1    S  0.6670221 0.3446230 -0.04179956 1.4177900

a = aggregate(midparent_fopg1 ~ Gene + line, data = gen2mid[gen2mid$pheno=='OUI',], FUN = mean)
a[,4] = aggregate(midparent_fopg1 ~ Gene + line, data = gen2mid[gen2mid$pheno=='OUI',], FUN = sd)[,3]
a[,5] = aggregate(midparent_fopg1 ~ Gene + line, data = gen2mid[gen2mid$pheno=='OUI',], FUN = min)[,3]
a[,6] = aggregate(midparent_fopg1 ~ Gene + line, data = gen2mid[gen2mid$pheno=='OUI',], FUN = max)[,3]
colnames(a) = c('Generation', 'Line','mean','sd','min','max')
a
# Generation Line       mean        sd        min      max
# 1       G2_1    R -0.7157259 0.7604569 -1.8356777 0.670206
# 2       G2_1    S  0.7072678 0.7964003 -0.9499933 2.025686

a = aggregate(midparent_fopg2 ~ Gene + line, data = gen2mid[gen2mid$pheno=='OUI',], FUN = mean)
a[,4] = aggregate(midparent_fopg2 ~ Gene + line, data = gen2mid[gen2mid$pheno=='OUI',], FUN = sd)[,3]
a[,5] = aggregate(midparent_fopg2 ~ Gene + line, data = gen2mid[gen2mid$pheno=='OUI',], FUN = min)[,3]
a[,6] = aggregate(midparent_fopg2 ~ Gene + line, data = gen2mid[gen2mid$pheno=='OUI',], FUN = max)[,3]
colnames(a) = c('Generation', 'Line','mean','sd','min','max')
a
#   Generation Line       mean        sd       min      max
# 1       G2_1    R -0.2068000 0.8865360 -1.618166 1.089737
# 2       G2_1    S  0.4397687 0.6604678 -1.372049 1.445655

a = aggregate(midparent_fopg ~ Gene + line, data = gen2mid[gen2mid$pheno=='OUI',], FUN = mean)
a[,4] = aggregate(midparent_fopg ~ Gene + line, data = gen2mid[gen2mid$pheno=='OUI',], FUN = sd)[,3]
a[,5] = aggregate(midparent_fopg ~ Gene + line, data = gen2mid[gen2mid$pheno=='OUI',], FUN = min)[,3]
a[,6] = aggregate(midparent_fopg ~ Gene + line, data = gen2mid[gen2mid$pheno=='OUI',], FUN = max)[,3]
colnames(a) = c('Generation', 'Line','mean','sd','min','max')
a
# Generation Line       mean        sd       min       max
# 1       G2_1    R -0.5537916 0.8182359 -1.959672 0.4387635
# 2       G2_1    S  0.6670180 0.8115003 -1.291064 2.0488685


## G x E interactions
a = aggregate(fecs ~ g + env +  time2, data = gbe[gbe$xp=='Chronic Stress effect',], FUN = mean)
a[,5] = aggregate(fecs ~ g + env + time2, data = gbe[gbe$xp=='Chronic Stress effect',], FUN = sd)[,4]
a[,6] = aggregate(fecs ~ g + env + time2, data = gbe[gbe$xp=='Chronic Stress effect',], FUN = min)[,4]
a[,7] = aggregate(fecs ~ g + env + time2, data = gbe[gbe$xp=='Chronic Stress effect',], FUN = max)[,4]
colnames(a) = c('Line', 'Environment', 'Time','mean','sd','min','max')
a
#    Line Environment              Time        mean        sd        min         max
# 1     R         nos Infection 1 - D24 -0.55134215 0.3533418 -0.7041585  0.97594440
# 2     S         nos Infection 1 - D24  0.63404347 1.1290791 -0.7041585  3.21048120
# 3     R         str Infection 1 - D24 -0.50771081 0.3998655 -0.7960883  0.46487806
# 4     S         str Infection 1 - D24  0.58386744 1.1638328 -0.7960883  3.88095070
# 5     R         nos Infection 1 - D30 -0.54625008 0.5014089 -0.8845604  0.88849996
# 6     S         nos Infection 1 - D30  0.62818760 1.0690285 -0.8732670  2.58250670
# 7     R         str Infection 1 - D30 -0.59382281 0.1712711 -0.7112452  0.05804918
# 8     S         str Infection 1 - D30  0.68289623 1.1219655 -0.6130374  3.96180876
# 9     R         nos Infection 2 - D24 -0.37177694 0.4998305 -0.6525654  1.03636137
# 10    S         nos Infection 2 - D24  0.42754349 1.2496072 -0.6525654  3.04698847
# 11    R         str Infection 2 - D24 -0.09765218 0.7980624 -0.4738671  1.96849588
# 12    S         str Infection 2 - D24  0.11230001 1.2034286 -0.4738671  4.34107707
# 13    R         nos Infection 2 - D30 -0.50382458 0.3444203 -0.7217051  0.46813238
# 14    S         nos Infection 2 - D30  0.57939826 1.1884952 -0.7184893  2.85852663
# 15    R         str Infection 2 - D30 -0.51037321 0.1557284 -0.5729138 -0.01789123
# 16    S         str Infection 2 - D30  0.58692919 1.2265848 -0.5729138  3.08098474

a = aggregate(fecs ~ g + env +  time2, data = gbe[gbe$xp=='Intestinal Species effect',], FUN = mean)
a[,5] = aggregate(fecs ~ g + env + time2, data = gbe[gbe$xp=='Intestinal Species effect',], FUN = sd)[,4]
a[,6] = aggregate(fecs ~ g + env + time2, data = gbe[gbe$xp=='Intestinal Species effect',], FUN = min)[,4]
a[,7] = aggregate(fecs ~ g + env + time2, data = gbe[gbe$xp=='Intestinal Species effect',], FUN = max)[,4]
colnames(a) = c('Line', 'Environment', 'Time','mean','sd','min','max')
a
#    Line Environment              Time        mean        sd          min       max
# 1     R          hc Infection 1 - D24 -0.75781170 0.4835239 -0.957198505 0.7532966
# 2     S          hc Infection 1 - D24  0.68203053 0.8393340 -0.957198505 1.7691932
# 3     R          tc Infection 1 - D24 -0.34062447 0.9406705 -2.728711315 1.3152432
# 4     S          tc Infection 1 - D24  0.35855208 0.9552344 -0.865107723 2.1654393
# 5     R          hc Infection 1 - D30 -0.60556109 1.0346006 -1.895895764 1.1729846
# 6     S          hc Infection 1 - D30  0.54500498 0.5740870 -0.970593520 1.5394623
# 7     R          tc Infection 1 - D30  0.08626636 0.9753445 -1.614371243 2.2622672
# 8     S          tc Infection 1 - D30 -0.09080670 1.0440267 -2.365487787 1.6043235
# 9     R          hc Infection 2 - D30 -0.70877280 1.0294952 -1.821449460 0.6840902
# 10    S          hc Infection 2 - D30  0.63789552 0.3079389  0.004477087 1.3087069
# 11    R          tc Infection 2 - D30 -0.34388665 1.2540869 -3.599198649 0.9900835
# 12    S          tc Infection 2 - D30  0.36198594 0.4257838 -0.743137004 0.9316094


