# Ambiente PROMETA no WSL

O ambiente `prometa-rnaseq` reúne Python/Jupyter, SRA Toolkit, Trim Galore/FastQC,
Salmon e R com tximport/DESeq2. O CSV mantém as colunas existentes
`sample_id,condition,study_id`.

## Preparação

1. No Ubuntu/WSL, instale Miniforge, Mambaforge ou micromamba caso ainda não tenha
   `conda`, `mamba` ou `micromamba` disponível.
2. Execute o instalador a partir do WSL:

   ```bash
   bash /mnt/c/Users/quint/OneDrive/Documentos/ChatGPT/PROMETA/setup_wsl.sh
   ```

   Ele copia o código e o CSV para `~/projects/PROMETA`, cria/atualiza o ambiente
   definido em `environment.yml` e registra o kernel Jupyter. Se o caminho Windows
   for diferente, defina `PROMETA_WINDOWS_PROJECT` antes de executar.

3. Os dados grandes devem ficar no filesystem Linux do WSL, dentro de
   `~/projects/PROMETA`, para evitar a penalidade de I/O em `/mnt/c`.
4. O instalador reutiliza `~/projects/tuta_rnaseq/reference/salmon_index` e
   `tuta_transcripts.fna`, se existirem, criando um link simbólico e gerando
   `reference/tx2gene.csv` (colunas `TXNAME,GENEID`) pelos IDs dos cabeçalhos.
   Se não encontrar essa referência, disponibilize um índice Salmon e um mapa
   TXNAME/GENEID correspondentes à mesma anotação.
5. Abra a pasta no VS Code a partir do terminal WSL:

   ```bash
   cd ~/projects/PROMETA
   code .
   ```

   Instale/ative a extensão **WSL** (Microsoft) no VS Code para abrir a janela no
   Ubuntu. No notebook, escolha o kernel **Python (PROMETA RNA-seq)**.

## Fluxo de execução

1. Rode a célula de pipeline no notebook. Os logs ficam em `logs/`, relatórios de
   QC em `qc/`, status em `run_status.csv` e quantificações em `quant/`.
2. Rode a célula de análise diferencial depois que as quantificações estiverem
   prontas. O script R analisa cada `study_id` como tratamento versus controle.
3. Rode a célula de resumo On-Score/On-Ratio após a análise diferencial.

O pipeline mantém FASTQ e arquivos SRA por padrão. Para removê-los após quantificar,
altere `KEEP_INTERMEDIATES = False` em `prometa_pipeline.py`.

> Cada linha do CSV é tratada como uma amostra independente. Confirme que os SRRs
> são réplicas biológicas distintas antes de interpretar os resultados do DESeq2.
