library(data.table)

Cases = fread("~/SFARI/cases.txt", header= T)

setnames(Cases, "id", "individual")

prs1 = fread("~/ALSPAC/PRSice2results/Sfarimergedfriendshipmtagprsice.all.score", header = TRUE)
prs1 = prs1[,c(2,10)]
setnames(prs1, 2, "friendship")

prs2 = fread("~/ALSPAC/PRSice2results/Sfarimergedempathyprsice.all.score", header = TRUE)
prs2 = prs2[,c(2, 10)]
setnames(prs2, 2, "empathy")


prs3 = fread("~/ALSPAC/PRSice2results/SfarimergedSQprsice.all.score", header = TRUE)
prs3 = prs3[,c(2, 10)]
setnames(prs3, 2, "systemizing")

phenotype = fread("~/SFARI/ssccore.txt")
phenotype2 = fread("~/SFARI/rrbsscore_proband_sfari.csv")

pca = fread("~/SFARI/liftOverPlink/files_imputed/SSC_cases_pca.eigenvec")
setnames(pca, 2, "IID")

srs_subscale = fread("~/SFARI/srs_subscale_score.txt")

prs_merged = merge(prs1, prs2, by = "IID")
prs_merged = merge(prs_merged, prs3, by = "IID")
merged = merge(Cases, prs_merged, by = "IID")
merged = merge(merged, pca, by = "IID")
merged = merge(merged, phenotype, by = "individual")
merged = merged[!duplicated(merged$FID),] 
merged = merge(merged, srs_subscale)


#ADOS-social-communication

a = lm(scale(ados_communication_social) ~ age_at_ados + ados_module + as.character(Sex) + as.character(dataset) + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + V15 + V16 + V17 + scale(ssc_diagnosis_full_scale_iq) , data = merged)

b = lm(scale(ados_communication_social) ~ age_at_ados + ados_module + as.character(Sex) + as.character(dataset) + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + V15 + V16 + V17  + scale(ssc_diagnosis_full_scale_iq) + scale(systemizing) , data = merged)

c = lm(scale(ados_communication_social) ~ age_at_ados + ados_module + as.character(Sex) + as.character(dataset) + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + V15 + V16 + V17  + scale(ssc_diagnosis_full_scale_iq) + scale(friendship) + scale(empathy) , data = merged)

d = lm(scale(ados_communication_social) ~ age_at_ados + ados_module + as.character(Sex) + as.character(dataset) + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + V15 + V16 + V17 + scale(ssc_diagnosis_full_scale_iq) + scale(friendship) + scale(empathy) + scale(systemizing) , data = merged)



anova(a,b)

anova(a,c)

anova(a,d)

anova(c,d)


#ADOS-restricted-repetitive

a = lm(scale(ados_restricted_repetitive) ~ age_at_ados + ados_module + as.character(Sex) + as.character(dataset) + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + V15 + V16 + V17  + scale(ssc_diagnosis_full_scale_iq) , data = merged)

b = lm(scale(ados_restricted_repetitive) ~ age_at_ados + ados_module + as.character(Sex) + as.character(dataset) + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + V15 + V16 + V17 + scale(ssc_diagnosis_full_scale_iq) + scale(systemizing) , data = merged)

c = lm(scale(ados_restricted_repetitive) ~ age_at_ados + ados_module + as.character(Sex) + as.character(dataset) + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12  + V13 + V14 + V15 + V16 + V17+ scale(ssc_diagnosis_full_scale_iq) + scale(friendship) + scale(empathy) , data = merged)

d = lm(scale(ados_restricted_repetitive) ~ age_at_ados + ados_module + as.character(Sex) + as.character(dataset) + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + V15 + V16 + V17 + scale(ssc_diagnosis_full_scale_iq) + scale(friendship) + scale(empathy) + scale(systemizing) , data = merged)



anova(a,b)

anova(a,c)

anova(a,d)




#RBS_r_overall_score

a = lm(scale(rbs_r_overall_score) ~ as.character(Sex) + as.character(dataset) + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12  + V13 + V14 + V15 + V16 + V17 + scale(ssc_diagnosis_full_scale_iq) , data = merged)

b = lm(scale(rbs_r_overall_score) ~ as.character(Sex) + as.character(dataset) + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12  +  V13 + V14 + V15 + V16 + V17 + scale(ssc_diagnosis_full_scale_iq) + scale(systemizing) , data = merged)

c = lm(scale(rbs_r_overall_score) ~ as.character(Sex) + as.character(dataset) + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12  + V13 + V14 + V15 + V16 + V17 + scale(ssc_diagnosis_full_scale_iq) + scale(friendship) + scale(empathy) , data = merged)

d = lm(scale(rbs_r_overall_score) ~ as.character(Sex) + as.character(dataset) + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12  + V13 + V14 + V15 + V16 + V17 +scale(ssc_diagnosis_full_scale_iq) + scale(friendship) + scale(empathy) + scale(systemizing) , data = merged)



anova(a,b)

anova(a,c)

anova(a,d)

anova(b,d)
