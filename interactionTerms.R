# script to demonstrate each case to construct DESeq2 design term

set.seed(1234)

library(DESeq2)

# Case 1: One factor (2 groups) --------------
dds <- makeExampleDESeqDataSet(n=10000,m=12)
names(colData(dds))[1] <- 'genotype'
dds$genotype <- as.factor(rep(c(rep("WT", 3), rep("Mutant",3)),2))
dds$genotype <- relevel(dds$genotype, "WT")
df <- as.data.frame(colData(dds))

model.matrix(~ genotype, data = df)

unique(dds$genotype)

design(dds) <- ~ genotype
dds <- DESeq(dds)
resultsNames(dds)


# Case 2: One factor (multiple groups) --------------
dds <- makeExampleDESeqDataSet(n=10000,m=12)
names(colData(dds))[1] <- 'genotype'
dds$genotype <- as.factor(rep(c(rep("I", 2), rep("II",2), rep("III", 2)),2))
dds$genotype <- relevel(dds$genotype, "I")
df <- as.data.frame(colData(dds))

model.matrix(~ genotype, data = df)

unique(dds$genotype)

design(dds) <- ~ genotype
dds <- DESeq(dds)
resultsNames(dds)



# Case 3: Two factors (2 groups each) --------------

dds <- makeExampleDESeqDataSet(n=10000,m=12)
names(colData(dds))[1] <- 'genotype'
dds$treatment <- as.factor(c(rep("Control",6), rep("Treatment",6)))
dds$genotype <- as.factor(rep(c(rep("WT", 3), rep("Mutant",3)),2))
dds$genotype <- relevel(dds$genotype, "WT")
df <- as.data.frame(colData(dds))

model.matrix(~ genotype + treatment, data = df)


design(dds) <- ~ genotype + treatment
dds <- DESeq(dds)
resultsNames(dds)



# Case 4: Two factors (multiple groups) --------------
dds <- makeExampleDESeqDataSet(n=10000,m=12)
names(colData(dds))[1] <- 'genotype'
dds$genotype <- as.factor(rep(c(rep("I", 2), rep("II",2), rep("III", 2)),2))
dds$treatment <- as.factor(c(rep("Control",6), rep("Treatment",6)))
dds$genotype <- relevel(dds$genotype, "I")

# create a merged group
dds$group <- factor(paste0(dds$genotype,'_', dds$treatment))
df <- as.data.frame(colData(dds))


design(dds) <- ~ group
dds <- DESeq(dds)
resultsNames(dds)


model.matrix(~ group, data = df)



## Example 2: two conditions, two genotypes, with an interaction term

dds <- makeExampleDESeqDataSet(n=100,m=12)
dds$genotype <- factor(rep(rep(c("I","II"),each=3),2))

design(dds) <- ~ genotype + condition + genotype:condition
dds <- DESeq(dds) 
resultsNames(dds)

colData(dds)

# the condition effect for genotype I (the main effect)
results(dds, contrast=c("condition","B","A"))

# the condition effect for genotype II
# this is, by definition, the main effect *plus* the interaction term
# (the extra condition effect in genotype II compared to genotype I).
results(dds, list( c("condition_B_vs_A","genotypeII.conditionB") ))

# the interaction term, answering: is the condition effect *different* across genotypes?
results(dds, name="genotypeII.conditionB")

model.matrix(~ genotype + condition + genotype:condition, data = colData(dds))


# Case 5: Two factors (multiple groups - interaction term) --------------
dds <- makeExampleDESeqDataSet(n=10000,m=12)
names(colData(dds))[1] <- 'genotype'
dds$genotype <- as.factor(rep(c(rep("I", 2), rep("II",2), rep("III", 2)),2))
dds$treatment <- as.factor(c(rep("Control",6), rep("Treatment",6)))
dds$genotype <- relevel(dds$genotype, "I")
df <- as.data.frame(colData(dds))

design(dds) <- ~ genotype + treatment + genotype:treatment
dds <- DESeq(dds)
resultsNames(dds)


model.matrix(~ genotype + treatment + genotype:treatment, data = df)


