library(TwoSampleMR)        #载入两个必要的R包
library(MRInstruments)
data(gwas_catalog)          #该部分是部分exposure GWAS summary statistic的数据
bmi_gwas <- subset(gwas_catalog, grepl("Speliotes", Author) & Phenotype == "Body mass index")      #利用gwas_catalog中的exposure数据
trait1.exposure_data <- format_data(bmi_gwas)
trait1.exposure_data <- read_exposure_data(                 #读取exposure 数据，利用自己的数据
  exposureFile,
  sep = ',',
  snp_col = "rsid",
    beta_col = "beta",
    se_col = "se",
    eaf_col = "maf_A2",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    pval_col = "P",
    units_col = "units",
    samplesize_col = "N",
    gene_col = "gene",
    id_col = "id",
    phenotype_col = 'phenotype'
  )
trait1.exposure_data = trait1.exposure_data[trait1.exposure_data$pval.exposure < 5e-8,]           #对得到的exposure数据进行按照P值筛选
trait1.exposure_data_clumped <- clump_data(trait1.exposure_data)                                  #clump数据，一般默认值就可以，主要去除过多的连锁不平衡位点
trait1.exposure_data_clumped['r2'] = get_r_from_pn(trait1.exposure_data_clumped$pval.exposure,trait1.exposure_data_clumped$samplesize.exposure)  ##计算R2值
trait1.trait2.outcome_data <- extract_outcome_data(snps = trait1.exposure_data_clumped$SNP,outcomes = "ieu-a-297")             #读取outcome数据，利用gwas_catalog里面的数据
trait1.trait2.outcome_data <- read_outcome_data(                                                            #读取自己的outcome的值
  snps = trait1.exposure_data_clumped$SNP,
  filename = outcomeFile,
  sep = ",",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "maf_A2",
  pval_col = "P",
  units_col = "Units",
  gene_col = "Gene",
  samplesize_col = "N",
  phenotype_col = "phenotype"
)

trait1.trait2.dat <- harmonise_data(exposure_dat = trait1.exposure_data_clumped, outcome_dat = trait1.trait2.outcome_data)     #harmonise exposure和outcome的数据，基本就是对其整合操作
trait1.trait2.results <- mr(trait1.trait2.dat)     #进行MR计算
trait1.trait2.single_snp_analysis <- mr_singlesnp(trait1.trait2.dat) #单个SNP分析
trait1.trait2.results.withOR = generate_odds_ratios(trait1.trait2.results) #MR结果整理并得出相应的置信区间
trait1.trait2.dat.mr_heterogeneity = mr_heterogeneity(trait1.trait2.dat) #敏感度分析
trait1.trait2.dat.mr_pleiotropy_test = mr_pleiotropy_test(trait1.trait2.dat) #敏感度分析
trait1.trait2.dat.run_mr_presso = run_mr_presso(trait1.trait2.dat)
trait1.trait2.single_snp_analysis.withOR = generate_odds_ratios(trait1.trait2.single_snp_analysis)
trait1.trait2.dat.mr_leaveoneout <- mr_leaveoneout(trait1.trait2.dat)
trait1.trait2.mr_steiger = mr_steiger2(r_exp = get_r_from_pn(trait1.trait2.dat$pval.exposure,trait1.trait2.dat$samplesize.exposure), r_out = get_r_from_pn(trait1.trait2.dat$pval.outcome,trait1.trait2.dat$samplesize.outcome), n_exp = trait1.trait2.dat$samplesize.exposure, n_out = trait1.trait2.dat$samplesize.outcome)
exp.R2 = trait1.trait2.mr_steiger$r2_exp
F.stat = (trait1.exposure_data_clumped$samplesize.exposure[1] - dim(trait1.exposure_data_clumped)[1] - 1) / (dim(trait1.exposure_data_clumped)[1]) * exp.R2 / ( 1 - exp.R2 )
trait1.trait2.results.withOR['F_stat'] = F.stat
trait1.trait2.results.withOR['power'] = ''
trait1.trait2.results.withOR['R2'] = exp.R2
trait1.trait2.results.withOR['egger_intercept'] = ''
trait1.trait2.results.withOR['MR_Egger.Q'] = ''
trait1.trait2.results.withOR['MR_Egger.Q_df'] = ''
trait1.trait2.results.withOR['MR_Egger.Q_pval'] = ''
trait1.trait2.results.withOR['Inverse_variance_weighted.Q'] = ''
trait1.trait2.results.withOR['Inverse_variance_weighted.Q_df'] = ''
trait1.trait2.results.withOR['Inverse_variance_weighted.Q_pval'] = ''
trait1.trait2.results.withOR['egger_intercept.pval'] = ''


trait1.trait2.results.withOR$egger_intercept[trait1.trait2.results.withOR$method == 'MR Egger'] = trait1.trait2.dat.mr_pleiotropy_test$egger_intercept
trait1.trait2.results.withOR$MR_Egger.Q[trait1.trait2.results.withOR$method == 'MR Egger'] = trait1.trait2.dat.mr_heterogeneity$Q[trait1.trait2.dat.mr_heterogeneity$method == 'MR Egger']
trait1.trait2.results.withOR$MR_Egger.Q_df[trait1.trait2.results.withOR$method == 'MR Egger'] = trait1.trait2.dat.mr_heterogeneity$Q_df[trait1.trait2.dat.mr_heterogeneity$method == 'MR Egger']
trait1.trait2.results.withOR$MR_Egger.Q_pval[trait1.trait2.results.withOR$method == 'MR Egger'] = trait1.trait2.dat.mr_heterogeneity$Q_pval[trait1.trait2.dat.mr_heterogeneity$method == 'MR Egger']
trait1.trait2.results.withOR$Inverse_variance_weighted.Q[trait1.trait2.results.withOR$method == 'Inverse variance weighted'] = trait1.trait2.dat.mr_heterogeneity$Q[trait1.trait2.dat.mr_heterogeneity$method == 'Inverse variance weighted']
trait1.trait2.results.withOR$Inverse_variance_weighted.Q_df[trait1.trait2.results.withOR$method == 'Inverse variance weighted'] = trait1.trait2.dat.mr_heterogeneity$Q_df[trait1.trait2.dat.mr_heterogeneity$method == 'Inverse variance weighted']
trait1.trait2.results.withOR$Inverse_variance_weighted.Q_pval[trait1.trait2.results.withOR$method == 'Inverse variance weighted'] = trait1.trait2.dat.mr_heterogeneity$Q_pval[trait1.trait2.dat.mr_heterogeneity$method == 'Inverse variance weighted']
trait1.trait2.results.withOR$egger_intercept.pval[trait1.trait2.results.withOR$method == 'MR Egger'] = trait1.trait2.dat.mr_pleiotropy_test$pval
trait1.trait2.results.withOR$Main_MR_results.Sd_Raw[trait1.trait2.results.withOR$method == 'Inverse variance weighted'] = trait1.trait2.dat.run_mr_presso[[1]]$`Main MR results`$Sd[1]
trait1.trait2.results.withOR$Main_MR_results.Sd_Corre[trait1.trait2.results.withOR$method == 'Inverse variance weighted'] = trait1.trait2.dat.run_mr_presso[[1]]$`Main MR results`$Sd[2]
trait1.trait2.results.withOR$Main_MR_results.pval_Raw[trait1.trait2.results.withOR$method == 'Inverse variance weighted'] = trait1.trait2.dat.run_mr_presso[[1]]$`Main MR results`$`P-value`[1]
trait1.trait2.results.withOR$Main_MR_results.pval_Corre[trait1.trait2.results.withOR$method == 'Inverse variance weighted'] = trait1.trait2.dat.run_mr_presso[[1]]$`Main MR results`$`P-value`[2]
trait1.trait2.results.withOR$MR_PRESSO.Results.Global.pval[trait1.trait2.results.withOR$method == 'Inverse variance weighted'] = trait1.trait2.dat.run_mr_presso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue
trait1.trait2.results.withOR$MR_PRESSO.Results.Distortion.pval[trait1.trait2.results.withOR$method == 'Inverse variance weighted'] = trait1.trait2.dat.run_mr_presso[[1]]$`MR-PRESSO results`$`Distortion Test`$Pvalue
trait1.trait2.results.withOR$isq[trait1.trait2.results.withOR$method == 'MR Egger'] = Isq(trait1.trait2.results$b[trait1.trait2.results$method == 'MR Egger'],trait1.trait2.results$se[trait1.trait2.results$method == 'MR Egger'])




mr_report(trait1.trait2.dat,output_path = reportSavePath)
