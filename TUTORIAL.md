# Tutorial Prático: Expressão Diferencial e Análise Funcional de RNA-seq
## Workshop RNA-seq UFRJ 2026
**Instrutora:** Mariana Feitosa Cavalheiro (mariana.cavalheiro@usp.br)

---

## 🎯 Objetivos da Aula
Neste tutorial, você aprenderá a realizar o passo final de um pipeline clássico de RNA-seq bulk:
1. Importar contagens brutas individuais geradas pelo **HTSeq-count**.
2. Criar e configurar o objeto de análise estatística no **DESeq2**.
3. Realizar controle de qualidade pós-normalização (PCA e Heatmap de Distâncias).
4. Identificar Genes Diferencialmente Expressos (DEGs) com significância estatística.
5. Visualizar resultados através de Volcano plots, MA plots e Heatmaps.
6. Dar contexto biológico aos DEGs por meio de **Análise de Enriquecimento Funcional** (ORA e GSEA).

---

## 📂 1. Preparação do Ambiente e Estrutura

Antes de abrir o RStudio, certifique-se de que a estrutura do seu projeto está correta. A estrutura recomendada para garantir a reprodutibilidade é:

```
2026-rnaseq-ufrj/
├── README.md
├── TUTORIAL.md                # Este guia
├── scripts/
│   └── analise_diferencial_htseq.R       # Seu script R
├── data/
│   ├── counts_htseq/          # Arquivos .counts individuais das amostras
│   └── metadata_htseq.csv     # Metadados das amostras
└── results/                   # Será criada automaticamente pelo R
    ├── figures/
    ├── tables/
    └── objects/
```

### O Dataset do Nosso Estudo
Utilizaremos dados públicos reais adaptados do estudo **GSE147507**, que investiga a infecção por **SARS-CoV-2** em células epiteliais pulmonares humanas primárias (NHBE) comparadas a células controle (Mock):
*   **Grupo Controle (Mock):** 3 réplicas biológicas (`Series1_NHBE_Mock_1`, `_2`, `_3`)
*   **Grupo Tratamento (Infectado):** 3 réplicas biológicas (`Series1_NHBE_SARS-CoV-2_1`, `_2`, `_3`)

---

## 💻 2. Passo a Passo do Script R

Abra o RStudio, carregue o projeto e crie/abra o script `scripts/analise_diferencial_htseq.R`. Vamos construir a análise etapa por etapa.

---

### Passo 0: Configuração de Diretórios e Pacotes

Primeiro, definimos os caminhos relativos (boas práticas!) e carregamos os pacotes necessários.

```r
# Definir caminhos relativos ao projeto
dir_data    <- "data/counts_htseq"
dir_results <- "results"
dir_figs    <- file.path(dir_results, "figures")
dir_tables  <- file.path(dir_results, "tables")
dir_objects <- file.path(dir_results, "objects")

# Criar os diretórios se não existirem
sapply(c(dir_figs, dir_tables, dir_objects), function(d) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
})

# Carregar pacotes
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)       # Anotação específica para Homo sapiens
library(AnnotationDbi)
library(dplyr)
library(pathview)
```

---

### Passo 1: Importação e Metadados do HTSeq-count

O **HTSeq-count** gera arquivos de contagem individuais por amostra. No R, importamos a tabela de metadados `metadata_htseq.csv` que aponta para cada um desses arquivos.

```r
# Ler a tabela de metadados das amostras
sample_table <- read.csv(
    file.path("data", "metadata_htseq.csv"),
    header = TRUE
)

# Definir a condição como um fator e explicitar que "controle" é o nível de referência
sample_table$condition <- factor(
    sample_table$condition,
    levels = c("controle", "tratamento")
)

# Visualizar a tabela estruturada
print(sample_table)
```

> **⚠️ IMPORTANTE:** Em análises de expressão diferencial, o nível de controle (referência) **deve** ser o primeiro nível do fator. Caso contrário, a direção do log2 Fold Change (LFC) ficará invertida!

---

### Passo 2: Construir o Objeto DESeqDataSet e Filtragem

Diferente do `featureCounts`, onde usamos `DESeqDataSetFromMatrix()`, para dados do HTSeq usamos a função dedicada **`DESeqDataSetFromHTSeqCount()`**:

```r
# Criar o objeto a partir dos arquivos individuais e metadados
dds <- DESeqDataSetFromHTSeqCount(
    sampleTable = sample_table,
    directory   = dir_data,
    design      = ~ condition
)

# Pré-filtragem: Manter apenas genes com pelo menos 10 reads no total de todas as amostras
# Isso reduz o ruído de genes com expressão quase nula e otimiza o tempo de processamento.
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]

cat("Genes mantidos após a filtragem:", nrow(dds), "\n")
```

---

### Passo 3: Executar a Análise de Expressão Diferencial

A função mágica `DESeq()` realiza as quatro etapas principais automaticamente:
1. **Estimação dos size factors:** Corrige a profundidade de sequenciamento das bibliotecas.
2. **Estimação de dispersão:** Mede a variabilidade natural de cada gene entre as réplicas biológicas.
3. **Ajuste do GLM:** Ajusta o modelo linear Binomial Negativo.
4. **Teste de Wald:** Compara estatisticamente a expressão entre tratamento e controle.

```r
# Executar o pipeline estatístico
dds <- DESeq(dds)

# Verificar os Size Factors gerados
print(sizeFactors(dds))

# Gerar e salvar o gráfico de estimativa de dispersão
png(file.path(dir_figs, "01_dispersao.png"), width = 800, height = 600)
plotDispEsts(dds, main = "Estimativa de Dispersão")
dev.off()
```

---

### Passo 4: Extração dos Resultados e Shrinkage

Extraímos os resultados estatísticos aplicando a correção de múltiplos testes (FDR) e o *LFC Shrinkage* (reduz a dispersão artificial de genes de baixa contagem).

```r
# Extrair resultados brutos (alpha = 0.05 define FDR desejada)
res <- results(dds, contrast = c("condition", "tratamento", "controle"), alpha = 0.05)

# Aplicar o encolhimento (Shrinkage) nativo para robustez do log2FC
res_shrunk <- lfcShrink(dds,
    contrast = c("condition", "tratamento", "controle"),
    type = "normal"
)

# Resumo estatístico
summary(res, alpha = 0.05)
```

---

### Passo 5: Controle de Qualidade Pós-Normalização (QC)

Usamos a transformação **VST (Variance Stabilizing Transformation)** para estabilizar a variância antes de rodar algoritmos de agrupamento (PCA/Heatmaps).

```r
# Transformação VST
vsd <- vst(dds, blind = TRUE)

# 1. Gráfico de PCA (Principal Component Analysis)
# Réplicas biológicas do mesmo grupo devem se agrupar próximas uma das outras!
pca_plot <- plotPCA(vsd, intgroup = "condition") +
    theme_minimal(base_size = 14) +
    scale_color_manual(values = c("controle" = "#2196F3", "tratamento" = "#F44336")) +
    labs(title = "PCA — Amostras por Condição", color = "Condição")
ggsave(file.path(dir_figs, "02_PCA.png"), pca_plot, width = 8, height = 6, dpi = 300)

# 2. Heatmap de Distância Euclidiana entre Amostras
sample_dists <- dist(t(assay(vsd)))
dist_matrix  <- as.matrix(sample_dists)
colors_heatmap <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

png(file.path(dir_figs, "03_heatmap_distancia.png"), width = 800, height = 700)
pheatmap(
    dist_matrix,
    clustering_distance_rows = sample_dists,
    clustering_distance_cols = sample_dists,
    col = colors_heatmap,
    main = "Distância Euclidiana entre Amostras (VST)"
)
dev.off()
```

---

### Passo 6: Visualizações dos Genes DE (MA & Volcano Plots)

O **MA Plot** e o **Volcano Plot** ajudam a entender a distribuição global das alterações de expressão gênica.

```r
# 1. MA Plot (Expressão Média vs Magnitude de Mudança)
png(file.path(dir_figs, "04_MA_plot.png"), width = 900, height = 600)
plotMA(res_shrunk, ylim = c(-5, 5), main = "MA Plot (com normal shrinkage)", colSig = "#E53935", alpha = 0.05)
abline(h = c(-1, 1), col = "gray50", lty = 2)
dev.off()

# 2. Volcano Plot (Magnitude vs Signiﬁcância Estatística)
volcano <- EnhancedVolcano(
    res,
    lab = rownames(res),
    x = "log2FoldChange",
    y = "pvalue",
    title = "Volcano Plot — Tratamento vs. Controle",
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 2.0,
    labSize = 3.5,
    col = c("grey70", "#4CAF50", "#2196F3", "#F44336"),
    legendPosition = "right",
    drawConnectors = TRUE
)
ggsave(file.path(dir_figs, "05_volcano_plot.png"), volcano, width = 10, height = 8, dpi = 300)
```

---

### Passo 7: Heatmap dos Top 50 Genes Expressos Diferencialmente

Criamos um Heatmap z-score para examinar os padrões de expressão individual dos 50 genes mais significativos.

```r
# Selecionar os top 50 genes ordenados por padj
res_ordered <- res[order(res$padj), ]
top_n <- 50
top_genes <- head(rownames(res_ordered[!is.na(res_ordered$padj), ]), top_n)

# Centralizar os dados de expressão (VST) por Z-score
mat <- assay(vsd)[top_genes, ]
mat_scaled <- t(scale(t(mat)))

# Gerar e salvar o Heatmap
png(file.path(dir_figs, "06_heatmap_top_genes.png"), width = 900, height = 1200, res = 150)
pheatmap(
    mat_scaled,
    clustering_method = "complete",
    show_rownames = TRUE,
    fontsize_row = 7,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    main = "Top 50 Genes Diferencialmente Expressos (Z-score)"
)
dev.off()
```

---

### Passo 8: Exportar Tabelas

Salvamos os resultados em arquivos CSV legíveis no Excel para compartilhamento ou futuras triagens rápidas.

```r
res_df <- as.data.frame(res_ordered)
res_df$gene_id <- rownames(res_df)

# Escrever a tabela completa
write.csv(res_df, file = file.path(dir_tables, "DE_resultados_completos.csv"), row.names = FALSE)

# Filtrar e exportar apenas os DEGs significativos (padj < 0.05 & |log2FC| > 1)
sig_genes <- res_df %>% filter(padj < 0.05, abs(log2FoldChange) > 1)
write.csv(sig_genes, file = file.path(dir_tables, "DE_genes_significativos.csv"), row.names = FALSE)
```

---

### Passo 9: Análise de Enriquecimento Funcional (ORA)

Agora, passamos da lista fria de genes para o entendimento biológico usando o pacote **`clusterProfiler`**.

```r
# Converter IDs oficiais de Gene Symbol para Entrez ID (essencial para KEGG/GO)
gene_entrez <- bitr(
    sig_genes$gene_id,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
)

# 1. ORA — Gene Ontology (Biological Process)
ora_go <- enrichGO(
    gene         = gene_entrez$ENTREZID,
    OrgDb        = org.Hs.eg.db,
    ont          = "BP",           # Biological Process
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    readable     = TRUE            # Retorna nomes de genes como símbolos legíveis
)

cat("Termos GO enriquecidos:", nrow(ora_go), "\n")

# Salvar gráficos de GO se houver enriquecimento
if (nrow(ora_go) > 0) {
    p_go_dot <- dotplot(ora_go, showCategory = 20, title = "GO — Biological Process (ORA)")
    ggsave(file.path(dir_figs, "07_GO_BP_dotplot.png"), p_go_dot, width = 10, height = 8, dpi = 300)
    
    # Gráfico de Rede Gene-Conceito (Conectividade Biológica)
    p_cnet <- cnetplot(ora_go, showCategory = 5, foldChange = NULL) + ggtitle("Rede Gene–Conceito (GO BP)")
    ggsave(file.path(dir_figs, "09_GO_BP_cnetplot.png"), p_cnet, width = 12, height = 10, dpi = 300)
}

# 2. ORA — KEGG Pathways
ora_kegg <- enrichKEGG(
    gene         = gene_entrez$ENTREZID,
    organism     = "hsa",          # 'hsa' para Homo sapiens
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1
)

cat("Vias KEGG enriquecidas:", nrow(ora_kegg), "\n")

# Salvar gráficos de KEGG se houver enriquecimento
if (nrow(ora_kegg) > 0) {
    p_kegg_dot <- dotplot(ora_kegg, showCategory = 20, title = "KEGG Pathways (ORA)")
    ggsave(file.path(dir_figs, "10_KEGG_dotplot.png"), p_kegg_dot, width = 10, height = 8, dpi = 300)
    
    # Salvar tabela com os resultados do KEGG
    write.csv(as.data.frame(ora_kegg), file = file.path(dir_tables, "KEGG_ORA_resultados.csv"), row.names = FALSE)
}
```

---

## 📝 3. Perguntas Práticas para Discussão

Durante ou após a execução da prática, tente responder às seguintes questões:

1. **Amostras PCA:** No gráfico de PCA gerado (`02_PCA.png`), os grupos controle e tratamento se separaram bem na primeira dimensão (PC1)? O que essa separação significa em termos biológicos?
2. **Número de DEGs:** Quantos genes foram regulados positivamente (Up-regulated) e negativamente (Down-regulated)? Existe alguma assimetria óbvia na resposta celular frente à infecção de SARS-CoV-2?
3. **Biologia da Resposta Humana:** Quais termos biológicos do Gene Ontology (GO) obtiveram maior enriquecimento no dotplot? O resultado faz sentido considerando que o tratamento foi uma infecção viral pulmonar? (Dica: procure por respostas imunológicas, antivirais e de interferon).
4. **Vias KEGG:** Quais vias metabólicas ou de sinalização do KEGG obtiveram maior enriquecimento? Há alguma relação direta com a patologia ou mecanismos de infecção pelo SARS-CoV-2 (por exemplo, vias de infecção viral, resposta inflamatória ou sinalização de citocinas)?

---
**🔬 Bom trabalho e excelente análise!**
