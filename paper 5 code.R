setwd("/Users/arthurauyeung/Desktop/Arthur/code/paper5")

library(ChemmineR)
library(data.table)

OBABEL <- "/opt/homebrew/bin/obabel"

# ============================================
# STEP 1: PRE-COMPUTE FINGERPRINTS (RUN ONCE)
# ============================================

precompute_fingerprints <- function(smiles_vector, output_file = "tcm_fps.rds") {
  
  if (file.exists(output_file)) {
    cat("✓ 指纹缓存已存在，加载中...\n")
    return(readRDS(output_file))
  }
  
  cat("预计算", length(smiles_vector), "个指纹...\n")
  cat("（只需运行一次）\n")
  
  writeLines(smiles_vector, "all_tcm.smi")
  
  # ONE OpenBabel call (not 310!)
  cmd <- paste(OBABEL, "all_tcm.smi -O all_tcm.sdf -gen2D")
  system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  sdf_set <- read.SDFset("all_tcm.sdf")
  ap_set <- sdf2ap(sdf_set)
  
  saveRDS(ap_set, output_file)
  unlink(c("all_tcm.smi", "all_tcm.sdf"))
  
  cat("✓ 缓存保存至:", output_file, "\n")
  return(ap_set)
}

# ============================================
# STEP 2: SIMILARITY CALCULATION (FIXED)
# ============================================

smiles_to_ap_query <- function(smiles, prefix = "query") {
  # Only for the query drug (1 compound)
  smi_file <- paste0(prefix, ".smi")
  sdf_file <- paste0(prefix, ".sdf")
  
  writeLines(smiles, smi_file)
  cmd <- paste(OBABEL, smi_file, "-O", sdf_file)
  system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  sdf <- read.SDFset(sdf_file)
  fp <- sdf2ap(sdf)
  
  unlink(c(smi_file, sdf_file))
  return(fp)
}

# FIXED: Calculate similarities one by one (proper S4 handling)
calculate_similarities <- function(query_fp, target_fps) {
  n <- length(target_fps)
  sim <- numeric(n)
  
  for (i in 1:n) {
    sim[i] <- cmp.similarity(query_fp, target_fps[i])
  }
  
  return(sim)
}

estimate_properties <- function(smiles) {
  c_count <- nchar(gsub("[^C]", "", smiles))
  o_count <- nchar(gsub("[^O]", "", smiles))
  n_count <- nchar(gsub("[^N]", "", smiles))
  logP_est <- (c_count - o_count - n_count) / 10
  
  hbd <- nchar(gsub("[^OH]", "", smiles)) + nchar(gsub("[^NH]", "", smiles))
  hba <- nchar(gsub("[^NO]", "", smiles))
  
  return(list(logP = round(logP_est, 2), HBD = min(hbd, 10), HBA = min(hba, 10)))
}

property_similarity <- function(p1, p2) {
  distance <- sqrt((abs(p1$logP - p2$logP)/5)^2 + 
                     (abs(p1$HBD - p2$HBD)/10)^2 + 
                     (abs(p1$HBA - p2$HBA)/10)^2)
  return(1 / (1 + distance))
}

# ============================================
# STEP 3: ANALYZE DRUG (FAST)
# ============================================

analyze_drug <- function(query_smiles, query_name, 
                         target_data, target_fps,
                         top_n = 15) {
  
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("分析药物:", query_name, "\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  start_time <- Sys.time()
  
  # Query properties
  query_props <- estimate_properties(query_smiles)
  cat("  logP:", query_props$logP, " HBD:", query_props$HBD, " HBA:", query_props$HBA, "\n")
  
  # Query fingerprint (ONE OpenBabel call per drug)
  cat("  生成查询指纹...\n")
  query_fp <- smiles_to_ap_query(query_smiles, paste0("q_", gsub(" ", "_", query_name)))
  
  # Calculate structural similarity (fixed)
  cat("  计算结构相似度...\n")
  struct_sim <- calculate_similarities(query_fp, target_fps)
  
  # Calculate property similarity
  cat("  计算理化性质相似度...\n")
  prop_sim <- numeric(nrow(target_data))
  for (i in 1:nrow(target_data)) {
    target_props <- estimate_properties(target_data$Canonical_smiles[i])
    prop_sim[i] <- property_similarity(query_props, target_props)
  }
  
  combined <- 0.7 * struct_sim + 0.3 * prop_sim
  
  top_idx <- order(combined, decreasing = TRUE)[1:top_n]
  
  results <- data.frame(
    rank = 1:top_n,
    name = target_data$Ingredient_name[top_idx],
    struct_sim = round(struct_sim[top_idx], 4),
    prop_sim = round(prop_sim[top_idx], 4),
    combined = round(combined[top_idx], 4)
  )
  
  elapsed <- round(Sys.time() - start_time, 2)
  cat("  完成! 用时:", elapsed, "秒\n")
  
  print(results)
  
  return(results)
}

# ============================================
# STEP 4: MAIN EXECUTION
# ============================================

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("TCM化合物相似性管道\n")
cat(paste(rep("=", 60), collapse=""), "\n")

# 加载数据
cat("\n[1/3] 加载数据...\n")
ingredients <- fread("HERB_ingredient_info_v2.txt", sep = "\t", data.table = FALSE)
ingredients_valid <- ingredients[!is.na(ingredients$Canonical_smiles) & 
                                   ingredients$Canonical_smiles != "", ]
cat("  ✓ 有效TCM化合物:", nrow(ingredients_valid), "\n")

# 预计算指纹 (只运行一次)
cat("\n[2/3] 预计算指纹...\n")
target_fps <- precompute_fingerprints(ingredients_valid$Canonical_smiles)

# 准备目标数据
target_data <- ingredients_valid[, c("Canonical_smiles", "Ingredient_name")]

# 分析药物
cat("\n[3/3] 分析药物...\n")

# Metformin
metformin_smiles <- "CN=C(N)N=C(N)N"
metformin_results <- analyze_drug(metformin_smiles, "Metformin", target_data, target_fps, top_n = 15)

# Warfarin
warfarin_smiles <- "CC(=O)C(CC1=CC=CC=C1)C2=C(OC3=CC=CC=C3C2=O)O"
warfarin_results <- analyze_drug(warfarin_smiles, "Warfarin", target_data, target_fps, top_n = 15)

# Berberine
berberine_smiles <- "COC1=C(C2=C[N+]3=C(C=C2C=C1)C4=CC5=C(C=C4C3)OCO5)OC"
berberine_results <- analyze_drug(berberine_smiles, "Berberine", target_data, target_fps, top_n = 15)

# 保存结果
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("保存结果\n")
cat(paste(rep("=", 60), collapse=""), "\n")

write.csv(metformin_results, "metformin_results.csv", row.names = FALSE)
write.csv(warfarin_results, "warfarin_results.csv", row.names = FALSE)
write.csv(berberine_results, "berberine_results.csv", row.names = FALSE)

cat("✓ metformin_results.csv\n")
cat("✓ warfarin_results.csv\n")
cat("✓ berberine_results.csv\n")

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("✓ 完成！\n")
cat(paste(rep("=", 60), collapse=""), "\n")

acenocoumarol_smiles <- "CC(=O)CC(C1=CC=C(C=C1)[N+](=O)[O-])C2=C(C3=CC=CC=C3OC2=O)O"
acenocoumarol_results <- analyze_drug(acenocoumarol_smiles, "acenocoumarol", target_data, target_fps, top_n = 15)

Anisindione_smiles <- "COC1=CC=C(C=C1)C1C(=O)C2=CC=CC=C2C1=O"
Anisindione_results <- analyze_drug(Anisindione_smiles, "Anisindione", target_data, target_fps, top_n = 15)

Apixaban_smiles <- "COC1=CC=C(C=C1)N2N=C(C(N)=O)C3=C2C(=O)N(CC3)C4=CC=C(C=C4)N5CCCCC5=O"
Apixaban_results <- analyze_drug(Apixaban_smiles, "Apixaban", target_data, target_fps, top_n = 15)

Rivaroxaban_smiles <- "ClC1=CC=C(S1)C(=O)NC[C@H]2CN(C(=O)O2)C3=CC=C(C=C3)N4CCOCC4=O"
Rivaroxaban_results <- analyze_drug(Rivaroxaban_smiles, "Rivaroxaban", target_data, target_fps, top_n = 15)

Milvexian_smiles <- "O=C1Nc2cnn(c2c2ccnc(c2)[C@H](CCC[C@H]1C)n1cnc(cc1=O)c1cc(Cl)ccc1n1nnc(c1)Cl)C(F)F"
Milvexian_results <- analyze_drug(Milvexian_smiles, "Milvexian", target_data, target_fps, top_n = 15)

Dabigatran_smiles <- "OC(=O)CCN(C(=O)c1ccc2c(c1)nc(n2C)CNc1ccc(cc1)C(=N)N)c1ccccn1"
Dabigatran_results <- analyze_drug(Dabigatran_smiles, "Dabigatran", target_data, target_fps, top_n = 15)

Edoxaban_smiles <- "CN(C)C(=O)[C@H]1CC[C@H](NC(=O)C(=O)NC2=NC=C(Cl)C=C2)[C@@H](C1)NC(=O)C3=NC4=C(CN(C)CC4)S3"
Edoxaban_results <- analyze_drug(Edoxaban_smiles, "Edoxaban", target_data, target_fps, top_n = 15)

Clopidogrel_smiles <- "COC(=O)[C@H](c1ccccc1Cl)N1CCc2sccc2C1"
Clopidogrel_results <- analyze_drug(Clopidogrel_smiles, "Clopidogrel", target_data, target_fps, top_n = 15)

Aspirin_smiles <- "CC(=O)Oc1ccccc1C(=O)O"
Aspirin_results <- analyze_drug(Aspirin_smiles, "Aspirin", target_data, target_fps, top_n = 15)

Dipyridamole_smiles <- "OCCN(CCO)C1=NC2=C(N=C(N=C2N2CCCCC2)N(CCO)CCO)C(=N1)N1CCCCC1"
Dipyridamole_results <- analyze_drug(Dipyridamole_smiles, "Dipyridamole", target_data, target_fps, top_n = 15)

Prasugrel_smiles <- "CC(=O)OC1=CC2=C(S1)CCN(C2)C(C3=CC=CC=C3F)C(=O)C4CC4"
Prasugrel_results <- analyze_drug(Prasugrel_smiles, "Prasugrel", target_data, target_fps, top_n = 15)

Ticagrelor_smiles <- "CCCSc1nc(N[C@@H]2C[C@H]2c2ccc(c(c2)F)F)c2c(n1)n(nn2)[C@@H]1C[C@@H]([C@H]([C@H]1O)O)OCCO"
Ticagrelor_results <- analyze_drug(Ticagrelor_smiles, "Ticagrelor", target_data, target_fps, top_n = 15)

Cilostazol_smiles <- "O=C1CCc2c(N1)ccc(c2)OCCCCc1nnnn1C1CCCCC1"
Cilostazol_results <- analyze_drug(Cilostazol_smiles, "Cilostazol", target_data, target_fps, top_n = 15)

Argatroban_smiles <- "C[C@@H]1CCN([C@H](C1)C(O)=O)C(=O)[C@H](CCCNC(N)=N)NS(=O)(=O)C1=CC=CC2=C1NCC(C)C2"
Argatroban_results <- analyze_drug(Argatroban_smiles, "Argatroban", target_data, target_fps, top_n = 15)

Tirofiban_smiles <- "CCCCS(=O)(=O)N[C@H](C(=O)O)Cc1ccc(cc1)OCCCCC1CCNCC1"
Tirofiban_results <- analyze_drug(Tirofiban_smiles, "Tirofiban", target_data, target_fps, top_n = 15)

Tanogitran_smiles <- "CC(C1=CC2=C(C=C1)N(C(=N2)CNC3=CC=C(C=C3)C(=N)N)C)(C(=O)N4CCCC4)NCC(=O)O"
Tanogitran_results <- analyze_drug(Tanogitran_smiles, "Tanogitran", target_data, target_fps, top_n = 15)

Asundexian_smiles <- "CCC(C(=O)NC1=CC(=C(C=C1)C(=O)N)F)N2C=C(C(=CC2=O)C3=C(C=CC(=C3)Cl)N4C=C(N=N4)C(F)(F)F)OC"
Asundexian_results <- analyze_drug(Asundexian_smiles, "Asundexian", target_data, target_fps, top_n = 15)

BMS_smiles <- "CN1CCN(C(=O)C1)C2=CC=CC3=C2CCN([C@@H]3C(=O)NC4=CC=C(C=C4)C(=O)O)C(=O)/C=C/C5=C(C=CC(=C5F)Cl)N6C=NN=N6"
BMS_results <- analyze_drug(BMS_smiles, "BMS-962212", target_data, target_fps, top_n = 15)

Bevifibatide_smiles <- "C1=CC2=C(C=C1)NC=C2C[C@H]3C(=O)N4CCC[C@H]4C(=O)N[C@@H](CSSCCC(=O)N[C@@H](CCCN=C(N)N)C(=O)NCC(=O)N[C@@H](CC(=O)O)C(=O)N3)C(=O)N"
Bevifibatide_results <- analyze_drug(Bevifibatide_smiles, "Bevifibatide", target_data, target_fps, top_n = 15)

Amlodipine_smiles <- "CCOC(=O)C1=C(NC(=C(C1C2=CC=CC=C2Cl)C(=O)OC)C)COCCN"
Amlodipine_results <- analyze_drug(Amlodipine_smiles, "Amlodipine", target_data, target_fps, top_n = 15)

Nifedipine_smiles <- "COC(=O)C1=C(C)NC(C)=C(C(=O)OC)C1c1ccccc1[N+](=O)[O-]"
Nifedipine_results <- analyze_drug(Nifedipine_smiles, "Nifedipine", target_data, target_fps, top_n = 15)

FelodipineER_smiles <- "CCOC(=O)C1=C(C)NC(C)=C(C(=O)OC)C1c1cccc(Cl)c1Cl"
FelodipineER_results <- analyze_drug(FelodipineER_smiles, "FelodipineER", target_data, target_fps, top_n = 15)

Diltiazem_smiles <- "COc1ccc([C@@H]2Sc3ccccc3N(CCN(C)C)C(=O)[C@@H]2OC(C)=O)cc1"
Diltiazem_results <- analyze_drug(Diltiazem_smiles, "Diltiazem", target_data, target_fps, top_n = 15)

Verapamil_smiles <- "COc1ccc(CCN(C)CCCC(C#N)(c2ccc(OC)c(OC)c2)C(C)C)cc1OC"
Verapamil_results <- analyze_drug(Verapamil_smiles, "Verapamil", target_data, target_fps, top_n = 15)

Lisinopril_smiles <- "NCCCC[C@H](N[C@@H](CCc1ccccc1)C(=O)O)C(=O)N1CCC[C@H]1C(=O)O.O.O"
Lisinopril_results <- analyze_drug(Lisinopril_smiles, "Lisinopril", target_data, target_fps, top_n = 15)

Ramipril_smiles <- "CCOC(=O)[C@H](CCc1ccccc1)N[C@@H](C)C(=O)N1[C@H](C(=O)O)C[C@@H]2CCC[C@@H]21"
Ramipril_results <- analyze_drug(Ramipril_smiles, "Ramipril", target_data, target_fps, top_n = 15)

Enalapril_smiles <- "CCOC(=O)[C@H](CCc1ccccc1)N[C@@H](C)C(=O)N1CCC[C@H]1C(=O)O"
Enalapril_results <- analyze_drug(Enalapril_smiles, "Enalapril", target_data, target_fps, top_n = 15)

Perindopril_smiles <- "CCC[C@H](N[C@@H](C)C(=O)N1[C@H](C(=O)O)C[C@@H]2CCCC[C@@H]21)C(=O)OCC"
Perindopril_results <- analyze_drug(Perindopril_smiles, "Perindopril", target_data, target_fps, top_n = 15)

Captopril_smiles <- "C[C@H](CS)C(=O)N1CCC[C@H]1C(=O)O"
Captopril_results <- analyze_drug(Captopril_smiles, "Captopril", target_data, target_fps, top_n = 15)

Losartan_smiles <- "CCCCc1nc(Cl)c(CO)n1Cc1ccc(-c2ccccc2-c2nnn[nH]2)cc1"
Losartan_results <- analyze_drug(Losartan_smiles, "Losartan", target_data, target_fps, top_n = 15)

Valsartan_smiles <- "CCCCC(=O)N(Cc1ccc(-c2ccccc2-c2nnn[nH]2)cc1)[C@H](C(=O)O)C(C)C"
Valsartan_results <- analyze_drug(Valsartan_smiles, "Valsartan", target_data, target_fps, top_n = 15)

Irbesartan_smiles <- "CCCCC1=NC2(CCCC2)C(=O)N1Cc1ccc(-c2ccccc2-c2nnn[nH]2)cc1"
Irbesartan_results <- analyze_drug(Irbesartan_smiles, "Irbesartan", target_data, target_fps, top_n = 15)

Candesartan_smiles <- "CCOc1nc2cccc(C(=O)O)c2n1Cc1ccc(-c2ccccc2-c2nnn[nH]2)cc1"
Candesartan_results <- analyze_drug(Candesartan_smiles, "Candesartan", target_data, target_fps, top_n = 15)

Telmisartan_smiles <- "CCCc1nc2c(C)cc(-c3nc4ccccc4n3C)cc2n1Cc1ccc(-c2ccccc2C(=O)O)cc1"
Telmisartan_results <- analyze_drug(Telmisartan_smiles, "Telmisartan", target_data, target_fps, top_n = 15)

Metoprolol_smiles <- "COCCc1ccc(OCC(O)CNC(C)C)cc1"
Metoprolol_results <- analyze_drug(Metoprolol_smiles, "Metoprolol", target_data, target_fps, top_n = 15)

Bisoprolol_smiles <- "CC(C)NCC(O)COc1ccc(COCCOC(C)C)cc1"
Bisoprolol_results <- analyze_drug(Bisoprolol_smiles, "Bisoprolol", target_data, target_fps, top_n = 15)

Atenolol_smiles <- "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1"
Atenolol_results <- analyze_drug(Atenolol_smiles, "Atenolol", target_data, target_fps, top_n = 15)

Carvedilol_smiles <- "COc1ccccc1OCCNCC(O)COc1cccc2[nH]c3ccccc3c12"
Carvedilol_results <- analyze_drug(Carvedilol_smiles, "Carvedilol", target_data, target_fps, top_n = 15)

Propranolol_smiles <- "CC(C)NCC(O)COc1cccc2ccccc12"
Propranolol_results <- analyze_drug(Propranolol_smiles, "Propranolol", target_data, target_fps, top_n = 15)

Gliclazide_smiles <- "Cc1ccc(S(=O)(=O)NC(=O)NN2CC3CCCC3C2)cc1"
Gliclazide_results <- analyze_drug(Gliclazide_smiles, "Gliclazide", target_data, target_fps, top_n = 15)

Glimepiride_smiles <- "CCC1=C(C)CN(C(=O)NCCc2ccc(S(=O)(=O)NC(=O)N[C@H]3CC[C@H](C)CC3)cc2)C1=O"
Glimepiride_results <- analyze_drug(Glimepiride_smiles, "Glimepiride", target_data, target_fps, top_n = 15)

Glipizide_smiles <- "Cc1cnc(C(=O)NCCc2ccc(S(=O)(=O)NC(=O)NC3CCCCC3)cc2)cn1"
Glipizide_results <- analyze_drug(Glipizide_smiles, "Glipizide", target_data, target_fps, top_n = 15)

Glyburide_smiles <- "COc1ccc(Cl)cc1C(=O)NCCc1ccc(S(=O)(=O)NC(=O)NC2CCCCC2)cc1"
Glyburide_results <- analyze_drug(Glyburide_smiles, "Glyburide", target_data, target_fps, top_n = 15)

Repaglinide_smiles <- "CCOc1cc(CC(=O)N[C@@H](CC(C)C)c2ccccc2N2CCCCC2)ccc1C(=O)O"
Repaglinide_results <- analyze_drug(Repaglinide_smiles, "Repaglinide", target_data, target_fps, top_n = 15)

Pioglitazone_smiles <- "CCc1ccc(CCOc2ccc(CC3SC(=O)NC3=O)cc2)nc1"
Pioglitazone_results <- analyze_drug(Pioglitazone_smiles, "Pioglitazone", target_data, target_fps, top_n = 15)

Acarbose_smiles <- "C[C@H]1O[C@H](O[C@H]2[C@H](O)[C@@H](O)[C@@H](O[C@H]3[C@H](O)[C@@H](O)C(O)O[C@@H]3CO)O[C@@H]2CO)[C@H](O)[C@@H](O)[C@@H]1N[C@H]1C=C(CO)[C@@H](O)[C@H](O)[C@H]1O"
Acarbose_results <- analyze_drug(Acarbose_smiles, "Acarbose", target_data, target_fps, top_n = 15)

Sitagliptin_smiles <- "N[C@@H](CC(=O)N1CCn2c(nnc2C(F)(F)F)C1)Cc1cc(F)c(F)cc1F"
Sitagliptin_results <- analyze_drug(Sitagliptin_smiles, "Sitagliptin", target_data, target_fps, top_n = 15)

Linagliptin_smiles <- "CC#CCn1c(N2CCC[C@@H](N)C2)nc2c1c(=O)n(Cc1nc(C)c3ccccc3n1)c(=O)n2C"
Linagliptin_results <- analyze_drug(Linagliptin_smiles, "Linagliptin", target_data, target_fps, top_n = 15)

Empagliflozin_smiles <- "OC[C@H]1O[C@@H](c2ccc(Cl)c(Cc3ccc(O[C@H]4CCOC4)cc3)c2)[C@H](O)[C@@H](O)[C@@H]1O"
Empagliflozin_results <- analyze_drug(Empagliflozin_smiles, "Empagliflozin", target_data, target_fps, top_n = 15)

Dapagliflozin_smiles <- "CCOc1ccc(Cc2cc([C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)ccc2Cl)cc1"
Dapagliflozin_results <- analyze_drug(Dapagliflozin_smiles, "Dapagliflozin", target_data, target_fps, top_n = 15)

Canagliflozin_smiles <- "Cc1ccc([C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)cc1Cc1ccc(-c2ccc(F)cc2)s1"
Canagliflozin_results <- analyze_drug(Canagliflozin_smiles, "Canagliflozin", target_data, target_fps, top_n = 15)

Semaglutide_smiles <- "CC[C@H](C)[C@@H](C(=O)N[C@@H](C)C(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CCCNC(=N)N)C(=O)NCC(=O)N[C@@H](CCCNC(=N)N)C(=O)NCC(=O)O)NC(=O)[C@H](CC3=CC=CC=C3)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CCCCNC(=O)COCCOCCNC(=O)COCCOCCNC(=O)CC[C@H](C(=O)O)NC(=O)CCCCCCCCCCCCCCCCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](C)NC(=O)[C@H](CCC(=O)N)NC(=O)CNC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC4=CC=C(C=C4)O)NC(=O)[C@H](CO)NC(=O)[C@H](CO)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](CC5=CC=CC=C5)NC(=O)[C@H]([C@@H](C)O)NC(=O)CNC(=O)[C@H](CCC(=O)O)NC(=O)C(C)(C)NC(=O)[C@H](CC6=CN=CN6)N"
Semaglutide_results <- analyze_drug(Semaglutide_smiles, "Semaglutide", target_data, target_fps, top_n = 15)

Liraglutide_smiles <- "CCCCCCCCCCCCCCCC(=O)N[C@@H](CCC(=O)NCCCC[C@H](NC(=O)[C@H](C)NC(=O)[C@H](C)NC(=O)[C@H](CCC(N)=O)NC(=O)CNC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](Cc1ccc(O)cc1)NC(=O)[C@H](CO)NC(=O)[C@H](CO)NC(=O)[C@@H](NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@@H](NC(=O)CNC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@@H](N)Cc1cnc[nH]1)[C@@H](C)O)[C@@H](C)O)C(C)C)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@H](C(=O)N[C@@H](C)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](CC(C)C)C(=O)N[C@H](C(=O)N[C@@H](CCCNC(=N)N)C(=O)NCC(=O)N[C@@H](CCCNC(=N)N)C(=O)NCC(=O)O)C(C)C)[C@@H](C)CC)C(=O)O"
Liraglutide_results <- analyze_drug(Liraglutide_smiles, "Liraglutide", target_data, target_fps, top_n = 15)

Etoformin_smiles <- "CCCCN/C(=N\\CC)NC(=N)N"
Etoformin_results <- analyze_drug(Etoformin_smiles, "Etoformin", target_data, target_fps, top_n = 15)

Exenatide_smiles <- "CC[C@H](C)[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H](NC(=O)[C@H](C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CCSC)NC(=O)[C@H](CCC(N)=O)NC(=O)[C@H](CCCCN)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@@H](NC(=O)CNC(=O)[C@H](CCC(=O)O)NC(=O)CNC(=O)[C@@H](N)Cc1cnc[nH]1)[C@@H](C)O)[C@@H](C)O)C(C)C)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CC(N)=O)C(=O)NCC(=O)NCC(=O)N1CCC[C@H]1C(=O)N[C@@H](CO)C(=O)N[C@@H](CO)C(=O)NCC(=O)N[C@@H](C)C(=O)N1CCC[C@H]1C(=O)N1CCC[C@H]1C(=O)N1CCC[C@H]1C(=O)N[C@@H](CO)C(N)=O"
Exenatide_results <- analyze_drug(Exenatide_smiles, "Exenatide", target_data, target_fps, top_n = 15)

Acetohexamide_smiles <- "CC(=O)c1ccc(S(=O)(=O)NC(=O)NC2CCCCC2)cc1"
Acetohexamide_results <- analyze_drug(Acetohexamide_smiles, "Acetohexamide", target_data, target_fps, top_n = 15)

Vildagliptin_smiles <- "N#C[C@@H]1CCCN1C(=O)CNC12CC3CC(CC(O)(C3)C1)C2"
Vildagliptin_results <- analyze_drug(Vildagliptin_smiles, "Vildagliptin", target_data, target_fps, top_n = 15)

Anagliptin_smiles <- "Cc1cc2ncc(C(=O)NCC(C)(C)NCC(=O)N3CCC[C@H]3C#N)cn2n1"
Anagliptin_results <- analyze_drug(Anagliptin_smiles, "Anagliptin", target_data, target_fps, top_n = 15)

Atorvastatin_smiles <- "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O"
Atorvastatin_results <- analyze_drug(Atorvastatin_smiles, "Atorvastatin", target_data, target_fps, top_n = 15)

Rosuvastatin_smiles <- "CC(C)c1nc(N(C)S(C)(=O)=O)nc(-c2ccc(F)cc2)c1/C=C/[C@@H](O)C[C@@H](O)CC(=O)O"
Rosuvastatin_results <- analyze_drug(Rosuvastatin_smiles, "Rosuvastatin", target_data, target_fps, top_n = 15)

Simvastatin_smiles <- "CCC(C)(C)C(=O)O[C@H]1C[C@@H](C)C=C2C=C[C@H](C)[C@H](CC[C@@H]3C[C@@H](O)CC(=O)O3)[C@H]21"
Simvastatin_results <- analyze_drug(Simvastatin_smiles, "Simvastatin", target_data, target_fps, top_n = 15)

Pravastatin_smiles <- "CC[C@H](C)C(=O)O[C@H]1C[C@H](O)C=C2C=C[C@H](C)[C@H](CC[C@@H](O)C[C@@H](O)CC(=O)O)[C@H]21"
Pravastatin_results <- analyze_drug(Pravastatin_smiles, "Pravastatin", target_data, target_fps, top_n = 15)

Fluvastatin_smiles <- "CC(C)n1c(/C=C/C(O)CC(O)CC(=O)O)c(-c2ccc(F)cc2)c2ccccc21"
Fluvastatin_results <- analyze_drug(Fluvastatin_smiles, "Fluvastatin", target_data, target_fps, top_n = 15)

Lovastatin_smiles <- "CC[C@H](C)C(=O)O[C@H]1C[C@@H](C)C=C2C=C[C@H](C)[C@H](CC[C@@H]3C[C@@H](O)CC(=O)O3)[C@H]21"
Lovastatin_results <- analyze_drug(Lovastatin_smiles, "Lovastatin", target_data, target_fps, top_n = 15)

Pitavastatin_smiles <- "O=C(O)C[C@H](O)C[C@H](O)/C=C/c1c(C2CC2)nc2ccccc2c1-c1ccc(F)cc1"
Pitavastatin_results <- analyze_drug(Pitavastatin_smiles, "Pitavastatin", target_data, target_fps, top_n = 15)

Mevastatin_smiles <- "CC[C@H](C)C(=O)O[C@H]1CCC=C2C=C[C@H](C)[C@H](CC[C@@H]3C[C@@H](O)CC(=O)O3)[C@H]21"
Mevastatin_results <- analyze_drug(Mevastatin_smiles, "Mevastatin", target_data, target_fps, top_n = 15)

Cerivastatin_smiles <- "COCc1c(C(C)C)nc(C(C)C)c(/C=C/[C@@H](O)C[C@@H](O)CC(=O)O)c1-c1ccc(F)cc1"
Cerivastatin_results <- analyze_drug(Cerivastatin_smiles, "Cerivastatin", target_data, target_fps, top_n = 15)

Bervastatin_smiles <- "CCOC(=O)CC(O)CC(O)/C=C/C1=C(c2ccc(F)cc2)c2ccccc2OC12CCCC2"
Bervastatin_results <- analyze_drug(Bervastatin_smiles, "Bervastatin", target_data, target_fps, top_n = 15)

Fenofibrate_smiles <- "CC(C)OC(=O)C(C)(C)Oc1ccc(C(=O)c2ccc(Cl)cc2)cc1"
Fenofibrate_results <- analyze_drug(Fenofibrate_smiles, "Fenofibrate", target_data, target_fps, top_n = 15)

Gemfibrozil_smiles <- "Cc1ccc(C)c(OCCCC(C)(C)C(=O)O)c1"
Gemfibrozil_results <- analyze_drug(Gemfibrozil_smiles, "Gemfibrozil", target_data, target_fps, top_n = 15)

Ezetimibe_smiles <- "O=C1[C@H](CC[C@H](O)c2ccc(F)cc2)[C@@H](c2ccc(O)cc2)N1c1ccc(F)cc1"
Ezetimibe_results <- analyze_drug(Ezetimibe_smiles, "Ezetimibe", target_data, target_fps, top_n = 15)

Bezafibrate_smiles <- "CC(C)(Oc1ccc(CCNC(=O)c2ccc(Cl)cc2)cc1)C(=O)O"
Bezafibrate_results <- analyze_drug(Bezafibrate_smiles, "Bezafibrate", target_data, target_fps, top_n = 15)

Laroprovstat_smiles <- "O=c1ccccn1-c1ccc(N[C@H]2CC[C@H](Nc3ncc(OC(F)F)cn3)C2)nc1"
Laroprovstat_results <- analyze_drug(Laroprovstat_smiles, "Laroprovstat", target_data, target_fps, top_n = 15)

Fluoxetine_smiles <- "CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1"
Fluoxetine_results <- analyze_drug(Fluoxetine_smiles, "Fluoxetine", target_data, target_fps, top_n = 15)

Sertraline_smiles <- "CN[C@H]1CC[C@@H](c2ccc(Cl)c(Cl)c2)c2ccccc21"
Sertraline_results <- analyze_drug(Sertraline_smiles, "Sertraline", target_data, target_fps, top_n = 15)

Paroxetine_smiles <- "Fc1ccc([C@@H]2CCNC[C@H]2COc2ccc3c(c2)OCO3)cc1"
Paroxetine_results <- analyze_drug(Paroxetine_smiles, "Paroxetine", target_data, target_fps, top_n = 15)

Paroxetine_smiles <- "Fc1ccc([C@@H]2CCNC[C@H]2COc2ccc3c(c2)OCO3)cc1"
Paroxetine_results <- analyze_drug(Paroxetine_smiles, "Paroxetine", target_data, target_fps, top_n = 15)

Citalopram_smiles <- "CN(C)CCCC1(c2ccc(F)cc2)OCc2cc(C#N)ccc21"
Citalopram_results <- analyze_drug(Citalopram_smiles, "Citalopram", target_data, target_fps, top_n = 15)

Escitalopram_smiles <- "CN(C)CCC[C@@]1(c2ccc(F)cc2)OCc2cc(C#N)ccc21"
Escitalopram_results <- analyze_drug(Escitalopram_smiles, "Escitalopram", target_data, target_fps, top_n = 15)

Fluvoxamine_smiles <- "COCCCC/C(=N\\OCCN)c1ccc(C(F)(F)F)cc1"
Fluvoxamine_results <- analyze_drug(Fluvoxamine_smiles, "Fluvoxamine", target_data, target_fps, top_n = 15)

Vortioxetine_smiles <- "Cc1ccc(Sc2ccccc2N2CCNCC2)c(C)c1"
Vortioxetine_results <- analyze_drug(Vortioxetine_smiles, "Vortioxetine", target_data, target_fps, top_n = 15)

Duloxetine_smiles <- "CNCC[C@H](Oc1cccc2ccccc12)c1cccs1"
Duloxetine_results <- analyze_drug(Duloxetine_smiles, "Duloxetine", target_data, target_fps, top_n = 15)

Venlafaxine_smiles <- "COc1ccc(C(CN(C)C)C2(O)CCCCC2)cc1"
Venlafaxine_results <- analyze_drug(Venlafaxine_smiles, "Venlafaxine", target_data, target_fps, top_n = 15)

Desvenlafaxine_smiles <- "CN(C)CC(c1ccc(O)cc1)C1(O)CCCCC1"
Desvenlafaxine_results <- analyze_drug(Desvenlafaxine_smiles, "Desvenlafaxine", target_data, target_fps, top_n = 15)

Mirtazapine_smiles <- "CN1CCN2c3ncccc3Cc3ccccc3C2C1"
Mirtazapine_results <- analyze_drug(Mirtazapine_smiles, "Mirtazapine", target_data, target_fps, top_n = 15)

Bupropion_smiles <- "CC(NC(C)(C)C)C(=O)c1cccc(Cl)c1"
Bupropion_results <- analyze_drug(Bupropion_smiles, "Bupropion", target_data, target_fps, top_n = 15)

Trazodone_smiles <- "O=c1n(CCCN2CCN(c3cccc(Cl)c3)CC2)nc2ccccn12"
Trazodone_results <- analyze_drug(Trazodone_smiles, "Trazodone", target_data, target_fps, top_n = 15)

Agomelatine_smiles <- "COc1ccc2cccc(CCNC(C)=O)c2c1"
Agomelatine_results <- analyze_drug(Agomelatine_smiles, "Agomelatine", target_data, target_fps, top_n = 15)

Tianeptine_smiles <- "CN1c2ccccc2C(NCCCCCCC(=O)O)c2ccc(Cl)cc2S1(=O)=O"
Tianeptine_results <- analyze_drug(Tianeptine_smiles, "Tianeptine", target_data, target_fps, top_n = 15)

Ibuprofen_smiles <- "CC(C)Cc1ccc(C(C)C(=O)O)cc1"
Ibuprofen_results <- analyze_drug(Ibuprofen_smiles, "Ibuprofen", target_data, target_fps, top_n = 15)

Naproxen_smiles <- "COc1ccc2cc([C@H](C)C(=O)O)ccc2c1"
Naproxen_results <- analyze_drug(Naproxen_smiles, "Naproxen", target_data, target_fps, top_n = 15)

Diclofenac_smiles <- "O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl"
Diclofenac_results <- analyze_drug(Diclofenac_smiles, "Diclofenac", target_data, target_fps, top_n = 15)

Indomethacin_smiles <- "COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1"
Indomethacin_results <- analyze_drug(Indomethacin_smiles, "Indomethacin", target_data, target_fps, top_n = 15)

Ketorolac_smiles <- "O=C(c1ccccc1)c1ccc2n1CCC2C(=O)O"
Ketorolac_results <- analyze_drug(Ketorolac_smiles, "Ketorolac", target_data, target_fps, top_n = 15)

Piroxicam_smiles <- "CN1C(C(=O)Nc2ccccn2)=C(O)c2ccccc2S1(=O)=O"
Piroxicam_results <- analyze_drug(Piroxicam_smiles, "Piroxicam", target_data, target_fps, top_n = 15)

Meloxicam_smiles <- "Cc1cnc(NC(=O)C2=C(O)c3ccccc3S(=O)(=O)N2C)s1"
Meloxicam_results <- analyze_drug(Meloxicam_smiles, "Meloxicam", target_data, target_fps, top_n = 15)

Celecoxib_smiles <- "Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1"
Celecoxib_results <- analyze_drug(Celecoxib_smiles, "Celecoxib", target_data, target_fps, top_n = 15)

Etoricoxib_smiles <- "Cc1ccc(-c2ncc(Cl)cc2-c2ccc(S(C)(=O)=O)cc2)cn1"
Etoricoxib_results <- analyze_drug(Etoricoxib_smiles, "Etoricoxib", target_data, target_fps, top_n = 15)

Parecoxib_smiles <- "CCC(=O)NS(=O)(=O)c1ccc(-c2c(-c3ccccc3)noc2C)cc1"
Parecoxib_results <- analyze_drug(Parecoxib_smiles, "Parecoxib", target_data, target_fps, top_n = 15)

MefenamicAcid_smiles <- "Cc1cccc(Nc2ccccc2C(=O)O)c1C"
MefenamicAcid_results <- analyze_drug(MefenamicAcid_smiles, "Mefenamic Acid", target_data, target_fps, top_n = 15)

Sulindac_smiles <- "CC1=C(CC(=O)O)c2cc(F)ccc2/C1=C\\c1ccc([S+](C)[O-])cc1"
Sulindac_results <- analyze_drug(Sulindac_smiles, "Sulindac", target_data, target_fps, top_n = 15)

Ketoprofen_smiles <- "CC(C(=O)O)c1cccc(C(=O)c2ccccc2)c1"
Ketoprofen_results <- analyze_drug(Ketoprofen_smiles, "Ketoprofen", target_data, target_fps, top_n = 15)

Flurbiprofen_smiles <- "CC(C(=O)O)c1ccc(-c2ccccc2)c(F)c1"
Flurbiprofen_results <- analyze_drug(Flurbiprofen_smiles, "Flurbiprofen", target_data, target_fps, top_n = 15)

TiaprofenicAcid_smiles <- "CC(C(=O)O)c1ccc(C(=O)c2ccccc2)s1"
TiaprofenicAcid_results <- analyze_drug(TiaprofenicAcid_smiles, "Tiaprofenic Acid", target_data, target_fps, top_n = 15)

Tenoxicam_smiles <- "CN1C(C(=O)Nc2ccccn2)=C(O)c2sccc2S1(=O)=O"
Tenoxicam_results <- analyze_drug(Tenoxicam_smiles, "Tenoxicam", target_data, target_fps, top_n = 15)

Lornoxicam_smiles <- "CN1C(C(=O)Nc2ccccn2)=C(O)c2sc(Cl)cc2S1(=O)=O"
Lornoxicam_results <- analyze_drug(Lornoxicam_smiles, "Lornoxicam", target_data, target_fps, top_n = 15)

Nabumetone_smiles <- "COc1ccc2cc(CCC(C)=O)ccc2c1"
Nabumetone_results <- analyze_drug(Nabumetone_smiles, "Nabumetone", target_data, target_fps, top_n = 15)

Acemetacin_smiles <- "COc1ccc2c(c1)c(CC(=O)OCC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1"
Acemetacin_results <- analyze_drug(Acemetacin_smiles, "Acemetacin", target_data, target_fps, top_n = 15)

Diflunisal_smiles <- "O=C(O)c1cc(-c2ccc(F)cc2F)ccc1O"
Diflunisal_results <- analyze_drug(Diflunisal_smiles, "Diflunisal", target_data, target_fps, top_n = 15)

Bucindolol_smiles <- "CC(C)(Cc1c[nH]c2ccccc12)NCC(O)COc1ccccc1C#N"
Bucindolol_results <- analyze_drug(Bucindolol_smiles, "Bucindolol", target_data, target_fps, top_n = 15)

Bevantolol_smiles <- "COc1ccc(CCNCC(O)COc2cccc(C)c2)cc1OC"
Bevantolol_results <- analyze_drug(Bevantolol_smiles, "Bevantolol", target_data, target_fps, top_n = 15)

Penbutolol_smiles <- "CC(C)(C)NC[C@H](O)COc1ccccc1C1CCCC1"
Penbutolol_results <- analyze_drug(Penbutolol_smiles, "Penbutolol", target_data, target_fps, top_n = 15)

Nebivolol_smiles <- "OC(CNCC(O)C1CCc2cc(F)ccc2O1)C1CCc2cc(F)ccc2O1"
Nebivolol_results <- analyze_drug(Nebivolol_smiles, "Nebivolol", target_data, target_fps, top_n = 15)

Esmolol_smiles <- "COC(=O)CCc1ccc(OCC(O)CNC(C)C)cc1"
Esmolol_results <- analyze_drug(Esmolol_smiles, "Esmolol", target_data, target_fps, top_n = 15)

Oxprenolol_smiles <- "C=CCOc1ccccc1OCC(O)CNC(C)C"
Oxprenolol_results <- analyze_drug(Oxprenolol_smiles, "Oxprenolol", target_data, target_fps, top_n = 15)

Carteolol_smiles <- "CC(C)(C)NCC(O)COc1cccc2c1CCC(=O)N2"
Carteolol_results <- analyze_drug(Carteolol_smiles, "Carteolol", target_data, target_fps, top_n = 15)

Labetalol_smiles <- "CC(CCc1ccccc1)NCC(O)c1ccc(O)c(C(N)=O)c1"
Labetalol_results <- analyze_drug(Labetalol_smiles, "Labetalol", target_data, target_fps, top_n = 15)

Nadolol_smiles <- "CC(C)(C)NCC(O)COc1cccc2c1CC(O)C(O)C2"
Nadolol_results <- analyze_drug(Nadolol_smiles, "Nadolol", target_data, target_fps, top_n = 15)

Sotalol_smiles <- "CC(C)NCC(O)c1ccc(NS(C)(=O)=O)cc1"
Sotalol_results <- analyze_drug(Sotalol_smiles, "Sotalol", target_data, target_fps, top_n = 15)

Pindolol_smiles <- "CC(C)NCC(O)COc1cccc2[nH]ccc12"
Pindolol_results <- analyze_drug(Pindolol_smiles, "Pindolol", target_data, target_fps, top_n = 15)

Timolol_smiles <- "CC(C)(C)NC[C@H](O)COc1nsnc1N1CCOCC1"
Timolol_results <- analyze_drug(Timolol_smiles, "Timolol", target_data, target_fps, top_n = 15)

Metipranolol_smiles <- "CC(=O)Oc1c(C)cc(OCC(O)CNC(C)C)c(C)c1C"
Metipranolol_results <- analyze_drug(Metipranolol_smiles, "Metipranolol", target_data, target_fps, top_n = 15)

Acebutolol_smiles <- "CCCC(=O)Nc1ccc(OCC(O)CNC(C)C)c(C(C)=O)c1"
Acebutolol_results <- analyze_drug(Acebutolol_smiles, "Acebutolol", target_data, target_fps, top_n = 15)

Celiprolol_smiles <- "CCN(CC)C(=O)Nc1ccc(OCC(O)CNC(C)(C)C)c(C(C)=O)c1"
Celiprolol_results <- analyze_drug(Celiprolol_smiles, "Celiprolol", target_data, target_fps, top_n = 15)

Betaxolol_smiles <- "CC(C)NCC(O)COc1ccc(CCOCC2CC2)cc1"
Betaxolol_results <- analyze_drug(Betaxolol_smiles, "Betaxolol", target_data, target_fps, top_n = 15)

Levobunolol_smiles <- "CC(C)(C)NC[C@H](O)COc1cccc2c1CCCC2=O"
Levobunolol_results <- analyze_drug(Levobunolol_smiles, "Levobunolol", target_data, target_fps, top_n = 15)

Mepindolol_smiles <- "Cc1cc2c(OCC(O)CNC(C)C)cccc2[nH]1"
Mepindolol_results <- analyze_drug(Mepindolol_smiles, "Mepindolol", target_data, target_fps, top_n = 15)

Bopindolol_smiles <- "Cc1cc2c(OCC(CNC(C)(C)C)OC(=O)c3ccccc3)cccc2[nH]1"
Bopindolol_results <- analyze_drug(Bopindolol_smiles, "Bopindolol", target_data, target_fps, top_n = 15)

Tertatolol_smiles <- "CC(C)(C)NCC(O)COc1cccc2c1SCCC2"
Tertatolol_results <- analyze_drug(Tertatolol_smiles, "Tertatolol", target_data, target_fps, top_n = 15)

Amoxicillin_smiles <- "CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccc(O)cc3)C(=O)N2[C@H]1C(=O)O"
Amoxicillin_results <- analyze_drug(Amoxicillin_smiles, "Amoxicillin", target_data, target_fps, top_n = 15)

Penicillin_G_smiles <- "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O"
Penicillin_G_results <- analyze_drug(Penicillin_G_smiles, "Penicillin G", target_data, target_fps, top_n = 15)

Cephalexin_smiles <- "CC1=C(C(=O)O)N2C(=O)[C@@H](NC(=O)[C@H](N)c3ccccc3)[C@H]2SC1.O"
Cephalexin_results <- analyze_drug(Cephalexin_smiles, "Cephalexin", target_data, target_fps, top_n = 15)

Cefuroxime_smiles <- "CO/N=C(\\C(=O)N[C@@H]1C(=O)N2C(C(=O)O)=C(COC(N)=O)CS[C@H]12)c1ccco1"
Cefuroxime_results <- analyze_drug(Cefuroxime_smiles, "Cefuroxime", target_data, target_fps, top_n = 15)

Ceftriaxone_smiles <- "CO/N=C(\\C(=O)N[C@@H]1C(=O)N2C(C(=O)O)=C(CSc3nc(=O)c(O)nn3C)CS[C@H]12)c1csc(N)n1"
Ceftriaxone_results <- analyze_drug(Ceftriaxone_smiles, "Ceftriaxone", target_data, target_fps, top_n = 15)

Ceftazidime_smiles <- "CC(C)(O/N=C(\\C(=O)N[C@@H]1C(=O)N2C(C(=O)[O-])=C(C[n+]3ccccc3)CS[C@H]12)c1csc(N)n1)C(=O)O"
Ceftazidime_results <- analyze_drug(Ceftazidime_smiles, "Ceftazidime", target_data, target_fps, top_n = 15)

Cefepime_smiles <- "CO/N=C(\\C(=O)N[C@@H]1C(=O)N2C(C(=O)[O-])=C(C[N+]3(C)CCCC3)CS[C@H]12)c1csc(N)n1"
Cefepime_results <- analyze_drug(Cefepime_smiles, "Cefepime", target_data, target_fps, top_n = 15)

Azithromycin_smiles <- "CC[C@H]1OC(=O)[C@H](C)[C@@H](O[C@H]2C[C@@](C)(OC)[C@@H](O)[C@H](C)O2)[C@H](C)[C@@H](O[C@@H]2O[C@H](C)C[C@H](N(C)C)[C@H]2O)[C@](C)(O)C[C@@H](C)CN(C)[C@H](C)[C@@H](O)[C@]1(C)O"
Azithromycin_results <- analyze_drug(Azithromycin_smiles, "Azithromycin", target_data, target_fps, top_n = 15)

Clarithromycin_smiles <- "CC[C@H]1OC(=O)[C@H](C)[C@@H](O[C@H]2C[C@@](C)(OC)[C@@H](O)[C@H](C)O2)[C@H](C)[C@@H](O[C@@H]2O[C@H](C)C[C@H](N(C)C)[C@H]2O)[C@](C)(OC)C[C@@H](C)C(=O)[C@H](C)[C@@H](O)[C@]1(C)O"
Clarithromycin_results <- analyze_drug(Clarithromycin_smiles, "Clarithromycin", target_data, target_fps, top_n = 15)

Erythromycin_smiles <- "CC[C@H]1OC(=O)[C@H](C)[C@@H](O[C@H]2C[C@@](C)(OC)[C@@H](O)[C@H](C)O2)[C@H](C)[C@@H](O[C@@H]2O[C@H](C)C[C@H](N(C)C)[C@H]2O)[C@](C)(O)C[C@@H](C)C(=O)[C@H](C)[C@@H](O)[C@]1(C)O"
Erythromycin_results <- analyze_drug(Erythromycin_smiles, "Erythromycin", target_data, target_fps, top_n = 15)

Ciprofloxacin_smiles <- "O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O"
Ciprofloxacin_results <- analyze_drug(Ciprofloxacin_smiles, "Ciprofloxacin", target_data, target_fps, top_n = 15)

Levofloxacin_smiles <- "C[C@H]1COc2c(N3CCN(C)CC3)c(F)cc3c(=O)c(C(=O)O)cn1c23.C[C@H]1COc2c(N3CCN(C)CC3)c(F)cc3c(=O)c(C(=O)O)cn1c23.O"
Levofloxacin_results <- analyze_drug(Levofloxacin_smiles, "Levofloxacin", target_data, target_fps, top_n = 15)

Moxifloxacin_smiles <- "COc1c(N2C[C@@H]3CCCN[C@@H]3C2)c(F)cc2c(=O)c(C(=O)O)cn(C3CC3)c12"
Moxifloxacin_results <- analyze_drug(Moxifloxacin_smiles, "Moxifloxacin", target_data, target_fps, top_n = 15)

Doxycycline_smiles <- "C[C@H]1c2cccc(O)c2C(=O)C2=C(O)[C@]3(O)C(=O)C(C(N)=O)=C(O)[C@@H](N(C)C)[C@@H]3[C@@H](O)[C@@H]21.O"
Doxycycline_results <- analyze_drug(Doxycycline_smiles, "Doxycycline", target_data, target_fps, top_n = 15)

Gentamicin_smiles <- "C(C1CCC(C(O1)OC2C(CC(C(C2O)OC3C(C(C(CO3)(C)O)NC)O)N)N)N)NC"
Gentamicin_results <- analyze_drug(Gentamicin_smiles, "Gentamicin", target_data, target_fps, top_n = 15)

Metronidazole_smiles <- "Cc1ncc([N+](=O)[O-])n1CCO"
Metronidazole_results <- analyze_drug(Metronidazole_smiles, "Metronidazole", target_data, target_fps, top_n = 15)

Vancomycin_smiles <- "CN[C@H](CC(C)C)C(=O)N[C@H]1C(=O)N[C@@H](CC(N)=O)C(=O)N[C@H]2C(=O)N[C@H]3C(=O)N[C@H](C(=O)N[C@H](C(=O)O)c4cc(O)cc(O)c4-c4cc3ccc4O)[C@H](O)c3ccc(c(Cl)c3)Oc3cc2cc(c3O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O[C@H]2C[C@](C)(N)[C@H](O)[C@H](C)O2)Oc2ccc(cc2Cl)[C@H]1O"
Vancomycin_results <- analyze_drug(Vancomycin_smiles, "Vancomycin", target_data, target_fps, top_n = 15)

Clindamycin_smiles <- "CCC[C@@H]1C[C@@H](C(=O)N[C@@H]([C@H]2O[C@H](SC)[C@H](O)[C@@H](O)[C@H]2O)[C@H](C)Cl)N(C)C1"
Clindamycin_results <- analyze_drug(Clindamycin_smiles, "Clindamycin", target_data, target_fps, top_n = 15)

Sulfamethoxazole_smiles <- "Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1"
Sulfamethoxazole_results <- analyze_drug(Sulfamethoxazole_smiles, "Sulfamethoxazole", target_data, target_fps, top_n = 15)
