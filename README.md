# RNA-seq Workshop — UFRJ 2026

## Análise de Expressão Diferencial de Genes

Workshop prático de bioinformática focado na análise de dados de RNA-seq, desde o controle de qualidade até a identificação de genes diferencialmente expressos e enriquecimento funcional.

**Instituição:** Universidade Federal do Rio de Janeiro (UFRJ)  
**Ano:** 2026  
**Instrutor(es):** Mariana Feitosa Cavalheiro

---

# RNA-seq Workshop — UFRJ 2026

## Análise de Expressão Diferencial de Genes

Workshop prático de bioinformática focado na análise de dados de RNA-seq, desde o controle de qualidade até a identificação de genes diferencialmente expressos e enriquecimento funcional.

**Instituição:** Universidade Federal do Rio de Janeiro (UFRJ)  
**Ano:** 2026  
**Instrutor(es):** Mariana Feitosa Cavalheiro (mariana.cavalheiro@usp.br)

---

## Sobre este Repositório (Template do Aluno)

Este repositório foi estruturado como um **ponto de partida limpo (starter kit)** para os estudantes. Para garantir o aprendizado ativo e evitar "respostas prontas" desde o início:
* As pastas de scripts (`scripts/`) e saídas geradas (`results/`) estão adicionadas ao `.gitignore` e não são versionadas.
* O repositório já vem com o conjunto de **dados de contagem simulados/reais** na pasta `data/` prontos para uso.
* O arquivo principal de orientação é o [TUTORIAL.md](file:///home/mfcaval/github/2026-rnaseq-ufrj/TUTORIAL.md), que serve como uma apostila/guia passo a passo para os alunos construírem seu próprio pipeline R.

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

> **Nota:** Quando o script R for executado pelos alunos, as seguintes pastas serão criadas de forma automatizada no ambiente local:
> * `scripts/` (onde o aluno criará/salvará o script da análise `analise_diferencial_htseq.R`)
> * `results/` (contendo as subpastas `figures/`, `tables/` e `objects/` para salvar saídas)

---

## Pré-requisitos

### Software

| Ferramenta       | Versão mínima | Finalidade                         |
|-----------------|--------------|--------------------------------------|
| R               | ≥ 4.3        | Análise estatística                  |
| RStudio         | ≥ 2024.04    | IDE para R (recomendado)             |
| FastQC          | ≥ 0.12       | Controle de qualidade das reads      |
| Trim Galore     | ≥ 0.6        | Trimming de adaptadores              |
| STAR            | ≥ 2.7        | Mapeamento splice-aware              |
| HTSeq           | ≥ 2.0        | htseq-count (quantificação)          |
| samtools        | ≥ 1.17       | Manipulação de BAM                   |

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

## Pipeline Completo (Referência)

As etapas anteriores à expressão diferencial (já realizadas no servidor de bioinformática):

```
FASTQ brutos → FastQC/MultiQC → Trim Galore → STAR (mapeamento) → htseq-count → [esta aula]
```

---

## Contato

Para dúvidas sobre o material, entre em contato com mariana.cavalheiro@usp.br.

---

## Licença

Este material é disponibilizado sob licença [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/).

