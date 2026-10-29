# Manual de Implementação e Execução do PROMETA

**Estrutura operacional baseada no ciclo CRISP-DM**
**Versão 1.0 — documento vivo**

Este manual orienta a implementação e a execução do PROMETA. O CRISP-DM organiza o trabalho em seis fases iterativas, dos objetivos científicos à disponibilização dos resultados. Cada fase descreve seu propósito, procedimentos, ferramentas, produtos esperados e verificações para avançar.

O caso inicial do PROMETA é a análise de dados públicos de RNA-seq de *Phthorimaea absoluta* exposta a abamectina, benzoato de emamectina e espinosade. A estrutura poderá ser adaptada a outros problemas biológicos, após revisar os dados, a referência, o desenho experimental e os métodos adequados ao novo contexto.

## Visão geral do ciclo

| Fase | Aplicação no PROMETA | Produto para avançar |
|---|---|---|
| 1. Compreensão do problema científico | Definir a pergunta, o escopo biológico e o uso dos resultados | Pergunta e critérios de sucesso aprovados |
| 2. Compreensão dos dados | Catalogar SRA, metadados, condições e réplicas; explorar qualidade e comparabilidade | Inventário e relatório exploratório das amostras |
| 3. Preparação dos dados | Processar reads, controlar qualidade, quantificar e preparar matrizes | Quantificações e metadados validados |
| 4. Modelagem | Executar DESeq2 por inseticida e calcular recorrência direcional entre grupos | Resultados diferenciais e núcleo de genes candidatos |
| 5. Avaliação | Verificar validade técnica, estatística e biológica | Resultados aprovados para interpretação e comunicação |
| 6. Implantação e comunicação | Organizar ambiente, scripts, parâmetros e entregáveis para repetição | Execução documentada e pacote de resultados rastreável |

O ciclo é iterativo. Se uma verificação falhar, retorne à fase responsável e atualize os produtos afetados.

## 1. Compreensão do problema científico

### 1.1 Contexto e escopo

A análise integra dados públicos de RNA-seq para investigar respostas moleculares de *P. absoluta* a três inseticidas. A execução inicial contempla os grupos de abamectina, benzoato de emamectina e espinosade. Cada grupo será comparado com seus controles disponíveis, e os resultados serão sintetizados para localizar genes com alterações significativas e direção consistente entre os tratamentos.

### 1.2 Pergunta operacional

**Quais genes apresentam expressão diferencial entre larvas tratadas e controles em cada grupo de inseticida, e quais mantêm a mesma direção de alteração nos três grupos?**

### 1.3 Objetivos de dados

- Processar as bibliotecas selecionadas por um fluxo documentado e repetível.
- Estimar expressão diferencial entre tratamento e controle separadamente por grupo de inseticida.
- Aplicar critérios estatísticos e de magnitude definidos antes da interpretação.
- Sintetizar a frequência e a direção das alterações para priorizar genes candidatos.

### 1.4 Critérios de sucesso

- Cada biblioteca está associada ao run, condição, inseticida e experimento corretos.
- Os controles e réplicas biológicas são adequados para cada comparação.
- Transcriptoma, índice Salmon e mapa `tx2gene` são compatíveis.
- As três comparações previstas produzem resultados revisáveis ou têm sua indisponibilidade justificada.
- Scripts, versões, parâmetros e critérios permitem reproduzir a análise.

### 1.5 Restrições e riscos

Os dados públicos podem diferir em população, protocolo, desenho e disponibilidade de controles. O campo `study_id` do metadado atual contém os nomes dos inseticidas, não os identificadores dos bioprojetos. A estrutura experimental de cada grupo precisa ser conferida antes da execução. A exigência de duas amostras por condição no script é apenas um filtro computacional mínimo; não garante poder estatístico adequado. Genes priorizados são candidatos computacionais e requerem validação independente para conclusões biológicas mais fortes.

### 1.6 Decisões a fechar

1. Confirmar se cada grupo de inseticida corresponde a um único bioprojeto e experimento.
2. Confirmar independência das réplicas e comparabilidade dos controles.
3. Definir a regra para grupos ou genes ausentes em alguma comparação.
4. Manter os critérios atuais de `padj ≤ 0,05` e `|log2FoldChange| ≥ 1` ou justificar alterações.
5. Fixar a fórmula, o nome e a interpretação do score de recorrência.

**Produto da fase:** pergunta, escopo, metadados necessários e critérios de sucesso documentados.

**Avançar quando:** as unidades amostrais, os grupos de comparação e o significado de “núcleo compartilhado” estiverem definidos.

## 2. Compreensão dos dados

### 2.1 Fontes e composição atual

O conjunto inicial contém 26 bibliotecas de RNA-seq do NCBI SRA, associadas aos bioprojetos PRJNA1359366, PRJNA749726 e PRJNA1270399. As bibliotecas correspondem a larvas de *P. absoluta* em condições controle ou expostas a abamectina, benzoato de emamectina ou espinosade.

O arquivo `metadata_geral.csv` contém uma linha por run e, atualmente, as colunas `sample_id`, `condition` e `study_id`. Os valores de `condition` são `control` e `treatment`; `study_id` identifica hoje o grupo de inseticida. Antes da análise, verifique se é preciso adicionar colunas como `bioproject_id`, `experiment_id`, `population_id` e `replicate_id`. Preserve compatibilidade com os scripts ou atualize-os junto com o esquema de metadados.

### 2.2 Exploração e auditoria

1. Conferir se cada `sample_id` existe no SRA e se a descrição do registro coincide com os metadados.
2. Confirmar espécie, estágio, tratamento, controle, layout da biblioteca e número de réplicas.
3. Contabilizar amostras por inseticida e condição.
4. Gerar e revisar relatórios de qualidade das reads e quantificações.
5. Explorar distribuição, correlação e agrupamento das amostras. Use transformação adequada para visualização; não use dados transformados como entrada direta do DESeq2.
6. Investigar outliers e variação associada a experimento, população ou lote antes de decidir exclusões ou ajustes.

PCA e correlação são ferramentas exploratórias. Elas ajudam a identificar padrões, mas não corrigem automaticamente efeitos de lote e não substituem o modelo estatístico.

### 2.3 Produtos da fase

- Inventário com identificadores SRA, bioprojeto, inseticida, condição e layout.
- Metadados revisados e contagem de amostras por comparação.
- Relatórios FastQC/Trim Galore e registro das amostras com falha ou exclusão.
- Descrição dos limites de comparabilidade entre experimentos.

**Avançar quando:** metadados, controles, réplicas e correspondência entre run e amostra estiverem auditados.

## 3. Preparação dos dados

### 3.1 Ambiente computacional

O ambiente está definido em `environment.yml` e foi preparado para WSL com Conda, Mamba ou Micromamba. O script `setup_wsl.sh` configura o ambiente `prometa-rnaseq`; `README_WSL.md` contém as instruções detalhadas. A execução ocorre em `~/projects/PROMETA`, onde devem ficar os arquivos grandes e temporários.

### 3.2 Pipeline de reads e quantificação

O fluxo configurado em `prometa_pipeline.py` executa:

1. `prefetch` e `fasterq-dump`, do SRA Toolkit, para baixar e converter runs em FASTQ.
2. Trim Galore! com FastQC para trimming e relatórios de qualidade.
3. Salmon para quantificar transcritos a partir das reads preparadas. O script usa validação de mapeamento e opções de correção de viés de sequência e GC.
4. Registro de logs, status, relatórios e arquivos `quant.sf` por run.

Os resultados são organizados em `sra_bruto/`, `raw_data/`, `quant/`, `logs/` e `qc/`. Esses dados intermediários e brutos não devem ser versionados no Git.

### 3.3 Referência e importação para análise por gene

Antes de quantificar, confirme que o transcriptoma, o índice Salmon e `reference/tx2gene.csv` correspondem à mesma espécie e anotação. O script de instalação pode reutilizar uma referência existente de outro projeto; confira os arquivos e os cabeçalhos antes de iniciar a análise.

O script `analise_diferencial.R` importa os arquivos `quant.sf` por meio do tximport, usando `countsFromAbundance = "lengthScaledTPM"`, `ignoreTxVersion = TRUE` e o mapeamento `TXNAME`–`GENEID`. Essa etapa prepara os dados no nível de gene para a análise diferencial.

### 3.4 Produtos e verificações

- FASTQ processados e relatórios de qualidade.
- Quantificações Salmon vinculadas aos metadados corretos.
- Transcriptoma, índice e `tx2gene.csv` verificados quanto à compatibilidade.
- Versões e parâmetros anotados.
- Nenhum arquivo SRA, FASTQ, índice grande ou segredo incluído no repositório.

**Avançar quando:** todas as bibliotecas incluídas tiverem quantificação válida, metadados coerentes e referência confirmada.

## 4. Modelagem

### 4.1 Expressão diferencial com DESeq2

O script `analise_diferencial.R` percorre os grupos registrados em `study_id` e compara `treatment` com `control` dentro de cada grupo. O desenho atual é `~ condition`. Genes cuja soma das contagens seja inferior a 10 são removidos antes do ajuste. O script requer pelo menos duas amostras em cada condição, executa DESeq2 e salva resultados com log2 fold change, valor-p ajustado e demais estatísticas em `differential_expression/`.

A análise separada evita combinar diretamente os três inseticidas como se fossem um único tratamento. Se um grupo reunir experimentos distintos ou tiver estrutura de lote, atualize os metadados e o desenho antes de estimar os efeitos.

### 4.2 Critérios de classificação

O notebook aplica os critérios `padj ≤ 0,05` e `|log2FoldChange| ≥ 1`. Entre os genes que atendem aos dois critérios, o sinal positivo indica indução e o sinal negativo indica repressão. Esses limiares são a configuração atual do repositório; registre qualquer alteração antes de interpretar os resultados.

### 4.3 Síntese de recorrência e definição do core

Para cada gene, o notebook conta em quantos grupos de inseticida ocorre alteração significativa em cada direção:

- `On_Score_UP`: número de grupos com indução significativa;
- `On_Score_DOWN`: número de grupos com repressão significativa;
- `Bono_Net_Score`: `On_Score_UP - On_Score_DOWN`;
- `On_Ratio_UP` e `On_Ratio_DOWN`: proporção de grupos com cada direção.

O núcleo compartilhado é separado por direção. Com três comparações válidas, `On_Score_UP = 3` define candidatos do núcleo induzido e `On_Score_DOWN = 3` define candidatos do núcleo reprimido. Se uma comparação não for válida, não classifique o gene como comum aos três grupos. A forma de tratar genes ausentes ou não testáveis deve ser registrada explicitamente.

### 4.4 Harmonização dos scripts

O notebook calcula o score a partir dos resultados diferenciais do DESeq2. O script `analise_meta_bono.py` calcula uma rota alternativa diretamente das médias de TPM, com pseudocontagem e limiar de log2-ratio, sem filtrar por `padj`. Essas rotas não são equivalentes. Use o caminho DESeq2 → critérios de classificação → score de recorrência como fluxo de referência. Atualize o script Python para corresponder a esse fluxo ou não o use na execução canônica.

### 4.5 Anotação funcional

A anotação funcional prevista usa DIAMOND blastx contra Swiss-Prot. Registre a versão da ferramenta e do banco, os parâmetros de busca e a regra para selecionar correspondências. Mantenha identificadores originais dos genes para permitir relacionar as anotações às saídas do DESeq2 e do score.

**Produtos da fase:** resultados diferenciais por inseticida, tabela de direção e recorrência por gene, listas do core induzido/reprimido e anotações funcionais.

## 5. Avaliação

### 5.1 Verificação técnica

- Todas as amostras esperadas foram processadas ou tiveram sua exclusão justificada.
- Quantificações, metadados e identificadores coincidem.
- Scripts, parâmetros e versões estão registrados.
- Uma nova execução com as mesmas entradas gera saídas equivalentes.

### 5.2 Verificação estatística

- Confirme o número de réplicas biológicas independentes por grupo.
- Revise os contrastes tratamento-controle e o fator de desenho.
- Verifique genes filtrados, valores ausentes e valores ajustados.
- Recalcule os scores a partir dos resultados diferenciais e confirme os denominadores usados nos On-Ratios.
- Examine a estabilidade do core diante de escolhas justificadas de limiares.

### 5.3 Verificação biológica e limitações

Compare os genes priorizados com a anotação e com evidências confiáveis sobre resposta a xenobióticos, detoxificação e estresse oxidativo. Diferencie recorrência estatística de confirmação funcional. A análise computacional prioriza candidatos; não comprova causalidade nem substitui validação experimental.

**Avançar para comunicação quando:** os resultados forem reproduzíveis, os critérios do core estiverem aplicados de forma consistente e as limitações dos dados estiverem documentadas.

## 6. Implantação e comunicação

### 6.1 Execução operacional

1. No WSL, acesse a pasta do projeto e execute `setup_wsl.sh` conforme `README_WSL.md`.
2. Confira metadados, referência Salmon e `tx2gene.csv`.
3. Execute `prometa_pipeline.py` para obter, preparar e quantificar as bibliotecas.
4. Execute `analise_diferencial.R` para obter as comparações tratamento-controle por inseticida.
5. Execute a etapa do notebook que calcula On-Scores, On-Ratios e o core compartilhado.
6. Revise os relatórios e salve as versões, parâmetros, logs e decisões da execução.

### 6.2 Pacote de resultados

Organize os resultados em tabelas e relatórios legíveis, acompanhados dos identificadores SRA, metadados, versões, parâmetros e logs relevantes. Mantenha FASTQ, arquivos SRA, índices e intermediários fora do Git. O repositório deve conter o código, o ambiente, a documentação e apenas os metadados apropriados para compartilhamento.

### 6.3 Reuso do fluxo

Para aplicar o PROMETA a outra espécie ou pergunta biológica, reinicie o ciclo a partir da compreensão do problema. Reavalie os metadados, a referência, o desenho, as etapas de preparação, o modelo estatístico e os critérios de validação. O funcionamento neste caso agrícola não valida automaticamente o fluxo em outros sistemas, como dados clínicos ou estudos de doenças.

## Glossário operacional

| Termo | Uso neste manual |
|---|---|
| BioProject | Registro que agrupa dados de um projeto no NCBI |
| Réplica biológica | Unidade biológica independente representada por uma amostra |
| DESeq2 | Pacote R usado para estimar expressão diferencial a partir de dados de RNA-seq |
| `padj` | Valor-p ajustado para testes múltiplos |
| `log2FoldChange` | Magnitude e direção da diferença de expressão em escala logarítmica base 2 |
| Core compartilhado | Genes que atendem aos critérios de significância e direção nos três grupos de inseticidas |
| Candidato | Gene priorizado para investigação posterior, ainda sem validação experimental implícita |
