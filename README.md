# PROMETA

**PROMETA** é um projeto de análise reprodutível de dados biológicos. Seu primeiro caso de uso é uma metanálise de dados públicos de RNA-seq de *Phthorimaea absoluta* exposta a inseticidas, para priorizar transcritos candidatos a biomarcadores. O CRISP-DM organiza o trabalho da pergunta científica à comunicação dos resultados.

O repositório reúne o fluxo computacional atualmente em desenvolvimento. A aplicação a outras espécies e a expansão para uma plataforma genérica são possibilidades futuras, que exigem adaptação e validação em cada novo conjunto de dados.

## Fluxo CRISP-DM

1. **Compreensão do problema científico** — delimitar a pergunta, os tratamentos e o que será considerado um candidato a biomarcador.
2. **Compreensão dos dados** — catalogar bibliotecas e metadados; verificar controles, réplicas e comparabilidade entre estudos.
3. **Preparação dos dados** — obter leituras do SRA, avaliar qualidade, realizar trimming, quantificar transcritos e organizar metadados e matrizes.
4. **Modelagem** — estimar expressão diferencial por estudo e resumir direção e recorrência dos sinais entre tratamentos com métricas como o Net XR-Score.
5. **Avaliação** — revisar qualidade, consistência, limitações e plausibilidade biológica dos resultados.
6. **Implantação e comunicação** — registrar scripts, parâmetros e versões; gerar resultados rastreáveis e relatórios para reprodução e revisão.

O ciclo é iterativo: limitações dos dados ou dos resultados podem levar à revisão das etapas anteriores.

## Ferramentas de análise

| Etapa | Ferramentas atualmente configuradas | Uso no projeto |
|---|---|---|
| Obtenção e preparação | SRA Toolkit, Trim Galore, FastQC e `pigz` | Obter runs públicos, inspecionar qualidade e preparar leituras |
| Quantificação | Salmon 1.10.3 | Estimar abundância de transcritos usando um índice de referência |
| Análise estatística | R, `tximport`, DESeq2, `readr` e `dplyr` | Importar quantificações e estimar expressão diferencial por estudo |
| Métricas e exploração | Python 3.12, Pandas, NumPy e JupyterLab | Organizar dados, calcular métricas de recorrência/direção e inspecionar resultados |
| Ambiente | Conda, Mamba ou Micromamba no WSL | Criar o ambiente descrito em `environment.yml` |

Nextflow, nf-core/rnaseq, contêineres, machine learning, WGCNA e dashboards aparecem em documentos de planejamento, mas não fazem parte do fluxo configurado neste repositório.

## Arquivos principais

- `metadata_geral.csv` — metadados com as colunas `sample_id`, `condition` e `study_id`.
- `prometa_pipeline.py` — automatiza obtenção do SRA, controle de qualidade, trimming e quantificação com Salmon.
- `analise_diferencial.R` — importa resultados do Salmon com `tximport` e executa DESeq2 por estudo.
- `analise_meta_bono.py` — calcula métricas de recorrência e direção a partir das abundâncias e metadados.
- `projeto_paralelo_vitor_quintela.ipynb` — notebook para executar e inspecionar etapas da análise.
- `environment.yml` — dependências do ambiente computacional.
- `setup_wsl.sh` — prepara o projeto e o ambiente no WSL.
- `README_WSL.md` — instruções técnicas detalhadas para configuração e execução no WSL.

## Configuração rápida no WSL

Com Conda, Mamba ou Micromamba instalado no Ubuntu/WSL, execute:

```bash
cd /mnt/c/Users/quint/OneDrive/Documentos/ChatGPT/PROMETA
bash setup_wsl.sh
cd ~/projects/PROMETA
```

Abra a pasta no VS Code a partir do WSL e selecione o kernel **Python (PROMETA RNA-seq)** no notebook. A sequência geral é: pipeline de processamento e quantificação, análise diferencial e resumo das métricas. Consulte [`README_WSL.md`](README_WSL.md) para detalhes e arquivos esperados.

> **Referência biológica:** o índice Salmon, o transcriptoma e o mapa `tx2gene.csv` precisam corresponder à espécie e à mesma anotação. Antes da análise, confirme que a referência configurada é adequada a *P. absoluta*; o instalador pode reutilizar uma referência existente de outro projeto.

## Reprodutibilidade e limites

Cada linha do CSV representa um run/amostra. Confirme que os runs correspondem a réplicas biológicas independentes e que cada estudo tem controles compatíveis antes de interpretar DESeq2 ou integrar resultados. Diferenças entre estudos não são corrigidas automaticamente apenas por aplicar PCA ou calcular métricas de recorrência.

Os transcritos priorizados são **candidatos computacionais**. A análise não demonstra causalidade nem valida biomarcadores experimentalmente. Dados brutos e arquivos intermediários são volumosos e devem ser armazenados fora do Git quando necessário; versionem-se os scripts, metadados apropriados e parâmetros, sem incluir dados restritos ou credenciais.
