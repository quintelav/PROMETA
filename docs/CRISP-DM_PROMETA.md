# Estrutura CRISP-DM do PROMETA

**Versão inicial — documento vivo**
Este documento organiza o PROMETA segundo o modelo CRISP-DM usado no material iHealth. Cada fase registra o propósito no projeto, o que já está descrito ou implementado, os produtos esperados e as perguntas que precisam ser resolvidas. O estado será atualizado à medida que a análise avançar.

O CRISP-DM é a estrutura de organização do trabalho, não uma ferramenta estatística nem um algoritmo. O estudo de *Phthorimaea absoluta* é o caso de demonstração do PROMETA; a possibilidade de adaptar o fluxo a outros sistemas biológicos é uma perspectiva que exige nova avaliação e validação.

## Visão geral das fases

| Fase | Situação atual no PROMETA | Produto principal |
|---|---|---|
| 1. Compreensão do problema científico | Pergunta e caso de uso definidos em versão inicial; critérios de sucesso precisam ser acordados | Pergunta científica, objetivos, escopo e critérios de sucesso |
| 2. Compreensão dos dados | 26 bibliotecas e três grupos de inseticidas registrados; auditoria exploratória precisa ser executada e documentada | Inventário de amostras, metadados e relatório de qualidade inicial |
| 3. Preparação dos dados | Pipeline SRA Toolkit → Trim Galore/FastQC → Salmon e importação com tximport estão configurados no repositório | Leituras processadas, quantificações, matriz por gene e registro dos parâmetros |
| 4. Modelagem | DESeq2 em R e síntese de recorrência no notebook estão implementados; o script Python alternativo precisa ser harmonizado | Resultados diferenciais por inseticida e tabela de consistência entre grupos |
| 5. Avaliação | Há resultados preliminares no resumo simples; precisam ser reconciliados com a execução final do DESeq2 | Avaliação estatística, biológica e das limitações do estudo |
| 6. Implantação e comunicação | Ambiente, scripts e documentação inicial estão no repositório; publicação dos dados brutos não faz parte do fluxo | Análise reproduzível, relatórios e documentação para reuso |

As fases são iterativas. Uma falha de qualidade pode levar à revisão da seleção ou preparação dos dados; resultados instáveis podem exigir rever o desenho, os critérios ou as comparações.

## 1. Compreensão do problema científico

### 1.1 Contexto

Dados públicos de RNA-seq permitem estudar respostas moleculares sem produzir novo sequenciamento, mas podem vir de experimentos com metadados, controles e condições diferentes. O PROMETA organiza as etapas computacionais necessárias para processar esses dados de forma rastreável e avaliar quais sinais se repetem entre grupos de tratamento.

O primeiro caso de uso investiga a resposta de larvas de *P. absoluta* à exposição a abamectina, benzoato de emamectina e espinosade. O pipeline de bioinformática é o foco técnico do artigo; a pergunta sobre a resposta do inseto demonstra o fluxo em uma situação biológica concreta.

### 1.2 Partes interessadas

| Parte interessada | Necessidade ou interesse |
|---|---|
| Pesquisador e orientador | Entender as decisões analíticas, revisar a interpretação biológica e reproduzir a análise |
| Pesquisadores de bioinformática | Acessar scripts, parâmetros, versões e critérios de processamento |
| Pesquisa e desenvolvimento agrícola | Receber resultados moleculares rastreáveis que possam orientar estudos posteriores |
| Futuros usuários do framework | Adaptar o fluxo a outra espécie ou pergunta, após revisar dados, referência e desenho analítico |

A aplicação em outros setores ou em dados clínicos não é demonstrada pelo estudo atual.

### 1.3 Objetivo científico traduzido em objetivo de dados

- **Objetivo científico:** investigar quais genes apresentam alteração de expressão em larvas de *P. absoluta* expostas a cada inseticida e quais apresentam direção consistente entre os três grupos.
- **Objetivo de dados:** processar as bibliotecas, estimar expressão diferencial tratamento-controle em cada grupo e sintetizar a recorrência dos genes que atendem aos critérios estatísticos e de magnitude definidos.
- **Critério operacional inicial de candidato:** `padj ≤ 0,05` e `|log2FoldChange| ≥ 1` em uma comparação. Para integrar os três grupos, a direção da alteração será contada separadamente como indução ou repressão.

Esses critérios são os que aparecem atualmente no notebook e devem ser justificados e congelados antes da análise final.

### 1.4 Pergunta de análise de dados

**Quais genes apresentam expressão diferencial entre larvas de *P. absoluta* tratadas e controles em cada grupo de inseticida, e quais mantêm uma direção de alteração consistente nos três grupos?**

A resposta computacional prioriza candidatos para investigação; não demonstra que um gene cause resistência nem que já seja um biomarcador validado.

### 1.5 Restrições e riscos

- Os dados provêm de estudos públicos distintos, com possível heterogeneidade de desenho, população, protocolo e controle.
- O metadado atual usa `study_id` para identificar o inseticida. Ele não registra separadamente BioProject, experimento, população ou lote; essa distinção precisa ser verificada antes de interpretar o desenho.
- O DESeq2 exige réplicas biológicas adequadas. O script requer pelo menos duas amostras por condição, mas esse limite computacional não garante poder estatístico suficiente.
- A referência do Salmon, o transcriptoma e o mapa `tx2gene` precisam corresponder à espécie e à mesma anotação.
- O efeito de lote não é removido automaticamente pelo CRISP-DM, pela PCA ou pelo score de recorrência.
- A análise não inclui validação experimental dos candidatos.
- A fórmula e a denominação finais do Net XR-Score precisam ser harmonizadas com a implementação usada.

### 1.6 Perguntas em aberto

1. Cada inseticida corresponde a um único BioProject/desenho experimental ou algum grupo combina estudos diferentes?
2. As amostras tratadas e controles são comparáveis dentro de cada grupo? Há réplicas biológicas independentes?
3. Qual fonte, versão e anotação correspondem ao índice Salmon e ao arquivo `tx2gene.csv`?
4. O critério de *core* será exigir consistência nos três inseticidas ou permitir recorrência em dois de três?
5. Como serão tratados genes ausentes ou não testáveis em alguma comparação?
6. Qual referência metodológica e fórmula serão citadas para o score de consistência?

## 2. Compreensão dos dados

### 2.1 Fontes e composição

O conjunto descrito no resumo simples inclui 26 bibliotecas públicas de RNA-seq do NCBI SRA, associadas aos bioprojetos PRJNA1359366, PRJNA749726 e PRJNA1270399. As amostras são de larvas de *P. absoluta* expostas a abamectina, benzoato de emamectina ou espinosade, além dos controles disponíveis.

O arquivo `metadata_geral.csv` contém uma linha por run e as colunas `sample_id`, `condition` e `study_id`. No arquivo atual, `study_id` assume os nomes dos inseticidas; recomenda-se renomear essa variável para `compound` ou acrescentar campos separados para BioProject e experimento, caso os dados permitam.

### 2.2 O que já está registrado

- A seleção totaliza 26 runs.
- O fluxo diferencia condições `control` e `treatment`.
- A análise está organizada em três grupos: emamectina, espinosade e abamectina.
- O layout single-end ou paired-end é detectado durante o processamento.

### 2.3 Análise exploratória a realizar

Antes da análise diferencial, deve ser produzido um relatório exploratório com:

- contagem de amostras por inseticida e condição;
- associação entre SRR, BioProject, tratamento, controle e réplica;
- layout, número de reads e métricas de qualidade antes e depois do trimming;
- distribuição das abundâncias e presença de amostras com valores ausentes ou inconsistentes;
- correlação e agrupamento das amostras, usando transformação apropriada apenas para visualização exploratória;
- detecção e investigação de possíveis outliers, sem remover amostras apenas por não agruparem como esperado.

A PCA pode ajudar a visualizar variação entre amostras, mas não substitui uma estratégia de correção de lote nem deve ser usada como teste de expressão diferencial.

### 2.4 Entregáveis e perguntas de exploração

**Entregáveis:** catálogo auditável de amostras, tabela de metadados revisada, relatório de controle de qualidade e registro de inclusão/exclusão de runs.

**Perguntas:** os grupos têm controles e réplicas suficientes? Os runs de uma mesma condição pertencem a um desenho comparável? A variação acompanha tratamento, BioProject, população ou características técnicas da biblioteca?

## 3. Preparação dos dados

### 3.1 Fluxo técnico configurado

O repositório contém o seguinte fluxo:

1. **SRA Toolkit** (`prefetch` e `fasterq-dump`) para obter as reads públicas e convertê-las em FASTQ.
2. **Trim Galore!** com **FastQC** para trimming e geração de relatórios de qualidade.
3. **Salmon** para quantificar transcritos contra um índice de referência, com validação de mapeamento e opções de viés habilitadas no script.
4. **tximport** em R para importar os arquivos `quant.sf`, agregar transcritos a genes por `tx2gene` e gerar a entrada do DESeq2.
5. Organização dos metadados e dos resultados por run e por grupo de inseticida.

### 3.2 O que precisa ser confirmado

- Correspondência biológica e de anotação entre transcriptoma, índice do Salmon e `tx2gene.csv`.
- Parâmetros e versões efetivamente usados em cada etapa.
- Qualidade e independência das réplicas biológicas.
- Tratamento de bibliotecas com baixa qualidade, reads insuficientes ou falhas de quantificação.
- Se BioProject, população ou lote precisam ser incorporados aos metadados e ao desenho estatístico.
- Que dados brutos e referências grandes permanecerão fora do Git; o repositório deve conter scripts, metadados permitidos, parâmetros e instruções, não FASTQ, SRA ou índices pesados.

### 3.3 Entregáveis

FASTQ preparados; quantificações do Salmon; metadados processados; matriz importada por gene; relatórios de QC; e registro de ferramentas, versões e parâmetros. Os entregáveis computacionais só serão considerados reproduzíveis quando forem ligados a uma versão identificável do pipeline e da referência.

## 4. Modelagem

### 4.1 Expressão diferencial com DESeq2

O script `analise_diferencial.R` realiza comparações separadas por grupo de inseticida. Ele importa quantificações do Salmon com tximport usando `countsFromAbundance = "lengthScaledTPM"`, filtra genes cuja soma das contagens seja inferior a 10, define o desenho `~ condition` e estima a diferença entre tratamento e controle com DESeq2. O script exige pelo menos duas amostras em cada condição e solicita o contraste tratamento versus controle com `alpha = 0,05`.

A análise por grupo evita tratar diretamente os três inseticidas como se fossem um único tratamento. Se um grupo reunir mais de um experimento ou lote, o desenho deverá ser revisto para representar essa estrutura, desde que existam controles e amostras suficientes para estimar os efeitos.

### 4.2 Score de recorrência e núcleo compartilhado

O notebook integra os arquivos diferenciais por `gene_id`. Em cada comparação, classifica como induzidos ou reprimidos os genes que satisfazem `padj ≤ 0,05` e `|log2FoldChange| ≥ 1`. O On-Score de cada direção conta em quantas comparações o gene satisfaz esses critérios; o Net Score é a diferença entre On-Score de indução e repressão; e o On-Ratio representa a proporção de comparações com a direção correspondente.

Como o objetivo é encontrar o núcleo de genes compartilhado nos três grupos, a definição operacional inicial será: On-Score de indução igual a 3 para o núcleo induzido e On-Score de repressão igual a 3 para o núcleo reprimido. Se uma comparação não puder ser realizada, o gene não será considerado parte de um núcleo de três grupos; a regra para comparações não disponíveis deverá ser informada.

### 4.3 Anotação funcional

O resumo simples descreve anotação funcional com DIAMOND blastx contra Swiss-Prot. Para tornar a etapa reproduzível, a execução final deverá registrar a versão do DIAMOND, a versão/data do banco, os parâmetros de busca e a regra de seleção das correspondências. O script de anotação e seus parâmetros também precisam ser incorporados à documentação do repositório, se essa análise for parte do artigo.

### 4.4 Harmonização necessária

O notebook calcula o score a partir das estatísticas do DESeq2. Já `analise_meta_bono.py` calcula uma versão alternativa diretamente das médias de TPM, com pseudocontagem e limiar de log2-ratio, sem usar `padj`. Essas rotas não são equivalentes. Para o artigo, a rota de referência será DESeq2 → classificação por significância e magnitude → score de recorrência. O script alternativo deve ser atualizado ou explicitamente excluído dos resultados apresentados.

## 5. Avaliação

### 5.1 Avaliação técnica

Verificar se todas as bibliotecas selecionadas foram processadas, se os arquivos de quantificação correspondem aos metadados, se o índice e a anotação são compatíveis e se os scripts produzem os mesmos resultados a partir dos mesmos arquivos e parâmetros. Registrar falhas e decisões de exclusão.

### 5.2 Avaliação estatística e biológica

- Revisar a adequação de controles, réplicas e desenho por grupo.
- Conferir os genes diferenciais e suas direções em cada inseticida.
- Verificar a recorrência do núcleo compartilhado e a sensibilidade dos resultados aos limiares definidos.
- Anotar e interpretar os candidatos com fontes científicas confiáveis.
- Distinguir consistência computacional de validação experimental.

O resumo simples relata candidatos associados a citocromos P450, peroxidases, carboxilesterases e proteínas cuticulares. Esses achados devem ser recalculados pela rota final baseada em DESeq2 antes de serem tratados como resultados definitivos.

### 5.3 Limitações e decisão de avanço

A heterogeneidade dos dados, a disponibilidade de controles e o número de réplicas limitam a inferência. O fluxo poderá avançar para interpretação apenas quando o desenho por grupo estiver documentado, a referência for confirmada, os resultados diferenciais forem reproduzíveis e a regra do score estiver fechada. A validação experimental fica fora do escopo computacional atual e será apresentada como etapa futura.

## 6. Implantação e comunicação

### 6.1 Como o PROMETA será utilizado neste estudo

O pipeline é executado no ambiente WSL descrito em `README_WSL.md`, com dependências definidas em `environment.yml`. O notebook coordena o processamento, chama o script R de expressão diferencial e calcula a síntese de recorrência. Os arquivos produzidos incluem quantificações, resultados diferenciais por grupo e tabelas consolidadas de score.

### 6.2 Entrega reprodutível

A divulgação do PROMETA deve incluir scripts, metadados compatíveis com compartilhamento, versões das ferramentas, parâmetros, identificadores de acesso das bibliotecas e instruções de execução. Os arquivos FASTQ/SRA, índices de referência e demais arquivos volumosos ficam fora do Git; os dados permanecem recuperáveis por seus identificadores públicos, quando disponíveis.

### 6.3 Reuso e alcance

O caso agrícola demonstra o fluxo em um sistema biológico concreto. Sua adaptação a outras pragas, culturas ou perguntas de bioinformática exigirá revisar a pergunta, os metadados, a referência, o desenho experimental, os modelos e os critérios de validação. O PROMETA não será descrito como validado para saúde humana, câncer ou outras doenças sem aplicações independentes nesses contextos.


## Glossário

| Termo | Significado no PROMETA |
|---|---|
| RNA-seq | Sequenciamento usado para estimar a abundância de transcritos em amostras biológicas |
| BioProject | Registro do NCBI que agrupa dados e amostras associados a um projeto de pesquisa |
| Réplica biológica | Amostra independente que representa uma unidade biológica no desenho experimental |
| Expressão diferencial | Diferença de expressão estimada entre tratamento e controle para cada gene |
| `padj` | Valor de significância ajustado para testes múltiplos |
| `log2FoldChange` | Tamanho e direção da diferença de expressão em escala logarítmica base 2 |
| Net XR-Score / On-Score | Métrica de recorrência e direção da regulação entre comparações; fórmula e nome final devem ser harmonizados |
| Núcleo compartilhado (*core*) | Genes que atendem aos critérios de significância e direção consistente nos três grupos definidos |
| Candidato a biomarcador | Gene priorizado computacionalmente para investigação; não implica validação experimental |
