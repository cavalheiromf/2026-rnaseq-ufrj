# RNA-seq Aula Prática — UFRJ 2026

## Análise de Expressão Diferencial de Genes

Aula prática de bioinformática focada na análise de dados de RNA-seq, iniciando a partir da quantificação de contagens (dados do HTSeq-count) até a identificação de genes diferencialmente expressos e enriquecimento funcional.

**Instituição:** Universidade Federal do Rio de Janeiro (UFRJ)  
**Ano:** 2026  
**Instrutor(es):** Mariana Feitosa Cavalheiro (mariana.cavalheiro@usp.br)

---

## Sobre este Repositório (Seu Template de Trabalho)

Este repositório foi estruturado como um **ponto de partida limpo (starter kit)** para você. Para garantir o seu aprendizado ativo e evitar "respostas prontas" desde o início:
* As pastas de scripts (`scripts/`) e saídas geradas (`results/`) estão adicionadas ao `.gitignore` e não são versionadas.
* O repositório já vem com o conjunto de **dados de contagem** na pasta `data/` prontos para uso.
* O arquivo principal de orientação é o [TUTORIAL.md](file:///home/mfcaval/github/2026-rnaseq-ufrj/TUTORIAL.md), que serve como uma apostila/guia passo a passo para você construir seu próprio pipeline R.

---

## Estrutura do Projeto

```
2026-rnaseq-ufrj/
├── README.md                  # Este arquivo explicativo
├── TUTORIAL.md                # Apostila e guia passo a passo da aula prática
├── .gitignore                 # Arquivo de configuração de ignorados (ignora scripts/ e results/)
├── slides/
│   └── aula_expressao_diferencial.html   # Apresentação de slides para abrir no navegador
└── data/                      # Pasta com dados preparados para a prática
    ├── counts_htseq/          # Arquivos individuais de contagem (.counts) do HTSeq-count (GSE147507)
    └── metadata_htseq.csv     # Tabela de metadados das amostras
```

> **Nota:** Quando você executar o script R, as seguintes pastas serão criadas de forma automatizada no seu ambiente local:
> * `scripts/` (onde você criará/salvará o seu script de análise `analise_diferencial_htseq.R`)
> * `results/` (contendo as subpastas `figures/`, `tables/` e `objects/` para você salvar as saídas)

---

## Pré-requisitos

### Software

Como as etapas de controle de qualidade, alinhamento e contagem já foram previamente processadas, você precisará apenas do R e RStudio instalados no seu computador local:

| Ferramenta       | Versão mínima | Finalidade                         |
|-----------------|--------------|--------------------------------------|
| R               | ≥ 4.3        | Análise estatística                  |
| RStudio         | ≥ 2024.04    | IDE para R (recomendado)             |

### Instalação de Pacotes R (Preparação do Computador)

Execute os comandos abaixo no console do RStudio para instalar as bibliotecas necessárias antes de iniciar a aula:

```r
# 1. Instalar Gerenciador do Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# 2. Instalar Pacotes do Bioconductor
BiocManager::install(c(
    "DESeq2",
    "clusterProfiler",
    "org.Hs.eg.db",      # Biblioteca de anotação humana
    "enrichplot",
    "pathview",
    "AnnotationDbi",
    "EnhancedVolcano"
))

# 3. Instalar Pacotes do CRAN
install.packages(c(
    "ggplot2",
    "pheatmap",
    "RColorBrewer",
    "ggrepel",
    "dplyr",
    "readr"
))
```

---

## Como usar este repositório durante a aula

1. **Clone ou baixe** este repositório em sua máquina.
2. Certifique-se de que os pacotes R listados acima estão instalados.
3. Abra o **RStudio** e configure o diretório de trabalho na raiz da pasta `2026-rnaseq-ufrj/`.
4. Crie uma pasta chamada `scripts/` e inicie um script R chamado `analise_diferencial_htseq.R`.
5. Abra e siga as instruções detalhadas contidas no guia [TUTORIAL.md](file:///home/mfcaval/github/2026-rnaseq-ufrj/TUTORIAL.md).
6. Os resultados (figuras e tabelas) serão automaticamente salvos na pasta local `results/` à medida que você executar o script.

---

## Pipeline Geral (Contexto)

As etapas prévias de bioinformática (já processadas e fornecidas como ponto de partida):

```
FASTQ brutos → Controle de Qualidade (FastQC) → Trimming (Trimommatic) → Mapeamento (STAR) → Quantificação (htseq-count) → [Início desta Aula (DESeq2 no R)]
```

---

## Contato

Para dúvidas sobre o material, entre em contato com mariana.cavalheiro@usp.br.

---

## Licença

Este material é disponibilizado sob licença [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/).
